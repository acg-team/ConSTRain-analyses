#!/usr/bin/env python3
import json

from cyvcf2 import VCF
import numpy as np
import pandas as pd

def dfs_from_vcf(filename: str, samples: list = None, vcf_format: str = "ConSTRain") -> pd.DataFrame:
    df_repeats = {
        "str_id": [],
        "chr": [],
        "start": [],
        "end": [],
        "unit": [],
        "period": [],
        "ref": [],
    }
    df_samples = {
        "sample": [],
        "str_id": [],
        "copy_number": [],
        "frequencies": [],
        "genotype": [],
    }

    try:
        parse_info_fields, parse_format_fields = VCF_FORMATS[vcf_format]
    except KeyError as e:
        return KeyError(f"Unsupported VCF format specified, use one of {list(VCF_FORMATS.keys())}")
    
    vcf = VCF(filename, samples=samples)
    for variant in vcf:           
        # df_repeats["chr"].append(variant.CHROM)
        # df_repeats["start"].append(variant.POS)
        # df_repeats["str_id"].append(f"{variant.CHROM}_{variant.POS}")

        df_repeats = parse_info_fields(df_repeats, variant)
        for sample_idx, sample in enumerate(vcf.samples):
            df_samples["sample"].append(sample)
            df_samples["str_id"].append(df_repeats["str_id"][-1]) # we should have just appended this
            df_samples = parse_format_fields(df_samples, variant, sample_idx)
            
    return pd.DataFrame(df_repeats), pd.DataFrame(df_samples)

def parse_constrain_info(df_repeats: pd.DataFrame, variant) -> pd.DataFrame:    
    df_repeats["chr"].append(variant.CHROM)
    df_repeats["start"].append(variant.POS)
    df_repeats["str_id"].append(f"{variant.CHROM}_{variant.POS}")
    df_repeats["end"].append(variant.INFO.get("END"))
    df_repeats["unit"].append(variant.INFO.get("RU"))
    df_repeats["period"].append(variant.INFO.get("PERIOD"))
    df_repeats["ref"].append(int(variant.INFO.get("REF")))

    return df_repeats

def parse_constrain_format(df_samples: pd.DataFrame, variant, sample_idx) -> pd.DataFrame:
    copy_number = variant.format("CN")[sample_idx]
    df_samples["copy_number"].append(copy_number[0])
    
    try:
        frequencies = variant.format("FREQS")[sample_idx]
        freq_dict = dict()
        for i in frequencies.split("|"):
            i = i.split(",")
            freq_dict[int(i[0])] = int(i[1])
        df_samples["frequencies"].append(freq_dict)            
    except (TypeError, IndexError):
        df_samples["frequencies"].append(np.nan)

    try:
        genotypes = variant.format("REPCN")[sample_idx]
        if genotypes == ".":
            raise ValueError
        genotypes = [int(i) for i in genotypes.split(",")]
        df_samples["genotype"].append(genotypes)
    except (TypeError, ValueError):
        df_samples["genotype"].append(np.nan)
    
    return df_samples

def parse_gangstr_info(df_repeats: pd.DataFrame, variant) -> pd.DataFrame:
    # GangSTR INFO fields are consistent with ConSTRain
    return parse_constrain_info(df_repeats, variant)

def parse_gangstr_format(df_samples: pd.DataFrame, variant, sample_idx) -> pd.DataFrame:
    df_samples["copy_number"].append(np.nan)
    try:
        frequencies = variant.format("ENCLREADS")[sample_idx]
        freq_dict = dict()
        for i in frequencies.split("|"):
            i = i.split(",")
            freq_dict[int(i[0])] = int(i[1])
        df_samples["frequencies"].append(freq_dict)            
    except (TypeError, IndexError):
        df_samples["frequencies"].append(np.nan)
    try:
        genotypes = variant.format("REPCN")[sample_idx]        
        if genotypes.size == 0: 
            raise TypeError
        # sometimes get int32 minimum value here (-2147483648), check for this
        if (genotypes < 0).sum() > 0:
            raise ValueError
        df_samples["genotype"].append(genotypes)
        # print(df_samples["genotype"])
    except (TypeError, ValueError):
        df_samples["genotype"].append(np.nan)

    return df_samples 

def parse_hipstr_info(df_repeats: pd.DataFrame, variant) -> pd.DataFrame:
    df_repeats["chr"].append(variant.CHROM)
    df_repeats["start"].append(variant.POS)
    df_repeats["str_id"].append(f"{variant.CHROM}_{variant.POS}")

    ref_allele = variant.REF
    period = variant.INFO.get("PERIOD")

    # Since HipSTR can extend the repeat region, we need to
    # do some work to find the original repeat unit
    offset = variant.INFO.get("START") - variant.POS
    unit = ref_allele[offset:offset+period]

    end = variant.POS + len(ref_allele) - 1

    df_repeats["end"].append(end)    
    df_repeats["unit"].append(unit)
    df_repeats["period"].append(period)
    df_repeats["ref"].append(len(ref_allele))

    return df_repeats

def parse_hipstr_format(df_samples: pd.DataFrame, variant, sample_idx) -> pd.DataFrame:
    # "Copy number"
    df_samples["copy_number"].append(np.nan)
    # "Allele frequencies {n basepairs: freq}"
    allele_lengths = [len(variant.REF)]
    for alt in variant.ALT:
        allele_lengths.append(len(alt))
    try:
        frequencies = variant.format("ALLREADS")[sample_idx]
        freq_dict = dict()
        for i in frequencies.split(";"):
            i = i.split("|")
            freq_dict[allele_lengths[0] + int(i[0])] = int(i[1])
        df_samples["frequencies"].append(freq_dict)            
    except (TypeError, IndexError):
        df_samples["frequencies"].append(np.nan)

    # "Genotype [n basepairs A, n basepairs B]"
    try:
        genotypes = [allele_lengths[i] for i in variant.genotypes[sample_idx][:-1]]
        df_samples["genotype"].append(genotypes)
    except (TypeError, ValueError):
        df_samples["genotype"].append(np.nan)
    
    return df_samples

def parse_eh_info(df_repeats: pd.DataFrame, variant) -> pd.DataFrame:
    df_repeats["chr"].append(variant.CHROM)
    df_repeats["start"].append(variant.POS + 1)
    df_repeats["str_id"].append(f"{variant.CHROM}_{variant.POS + 1}")
    df_repeats["end"].append(variant.INFO.get("END"))
    df_repeats["unit"].append(variant.INFO.get("RU"))
    df_repeats["period"].append(len(variant.INFO.get("RU")))
    df_repeats["ref"].append(int(variant.INFO.get("REF")))

    return df_repeats

def parse_eh_format(df_samples: pd.DataFrame, variant, sample_idx) -> pd.DataFrame:
    # "copy_number": [x],
    # "frequencies": [x],
    # "genotype": [x],
    df_samples["copy_number"].append(np.nan)

    allele_lengths = [int(variant.INFO.get("REF"))]
    for alt in variant.ALT:
        if alt == ".":
            continue
        allele_len = int(alt.replace("<STR", "").replace(">", ""))
        allele_lengths.append(allele_len)

    try:
        coverage = round(float(variant.format("LC")[sample_idx]))
        df_samples["frequencies"].append(coverage)

        # frequencies = variant.format("ADSP")[sample_idx]
        # freq_dict = dict()
        # for i, freq in enumerate(frequencies.split("/")):
        #     freq_dict[allele_lengths[i]] = int(freq)
        # df_samples["frequencies"].append(freq_dict)            
    except (ValueError, TypeError, IndexError):
        df_samples["frequencies"].append(np.nan)

    try:
        genotypes = [allele_lengths[i] for i in variant.genotypes[sample_idx][:-1]]
        df_samples["genotype"].append(genotypes)
    except (TypeError, ValueError):
        df_samples["genotype"].append(np.nan)

    return df_samples


VCF_FORMATS = {
    "ConSTRain": (parse_constrain_info, parse_constrain_format),
    "GangSTR": (parse_gangstr_info, parse_gangstr_format),
    "HipSTR": (parse_hipstr_info, parse_hipstr_format),
    "ExpansionHunter": (parse_eh_info, parse_eh_format),
}

def dfs_from_EH_json(filename: str) -> pd.DataFrame:
    df_repeats = {
        "str_id": [],
        "chr": [],
        "start": [],
        "end": [],
        "unit": [],
        "period": [],
        "ref": [],
    }

    df_samples = {
        "sample": [],
        "str_id": [],
        "copy_number": [],
        "frequencies": [],
        "genotype": [],
    }

    data = None
    with open(filename, 'r') as f:
        data = json.load(f)

    sample = data["SampleParameters"]["SampleId"]
    for record in data["LocusResults"]:
        str_id = record.keys()[0]
        df_repeats["str_id"].append(str_id)

        chrom, coords = record[str_id]["Variants"][str_id]["ReferenceRegion"].split(":")
        start, end = coords.split("-")
        start, end = int(start), int(end)
        df_repeats["chr"].append(chrom)
        df_repeats["start"].append(start)
        df_repeats["end"].append(end)


        unit = record[str_id]["Variants"][str_id]["ReferenceUnit"]
        df_repeats["unit"].append(unit)
        df_repeats["period"].append(len(unit))
        df_repeats["ref"].append(int(record[str_id]["AlleleCount"]))

        df_samples["sample"].append(sample)
        df_samples["str_id"].append(str_id)
        df_samples["copy_number"].append(np.nan)
        freq_dict = dict()
        for i in eval(record[str_id]["Variants"][str_id]["CountsOfSpanningReads"]):
            freq_dict[i[0]] = i[1]
        df_samples["frequencies"].append(freq_dict)
        genotype = record[str_id]["Variants"][str_id]["CountsOfSpanningReads"].split("/")
        df_samples["genotype"].append(genotype)

def load_haplotype_vars(filename: str, min_mapq=20) -> pd.DataFrame:
    colnames = ["chr", "start", "end", "qdepth", "mapq", "REF", "ALT", "qname", "qstart", "qend", "qorientation"]
    
    df_haplotype = (
        pd.read_csv(filename, comment="R", sep="\t", names=colnames)
            .query(f"qdepth == 1 and mapq >= {min_mapq}")
            .reset_index(drop=True)[["chr", "start", "end", "REF", "ALT"]]
    )

    return df_haplotype

def load_mreps_repeats(filename: str, seqname: str):    
    # Line of dashes shows up right before and right after mreps reports repeats
    dashes = " ---------------------------------------------------------------------------------------------\n"
    df_repeats = {
        "str_id": [],
        "chr": [],
        "start": [],
        "end": [],
        "unit": [],
        "period": [],
        "ref": [],
    }

    with open(filename, 'r') as f:
        parse_lines = False
        for line in f:
            if line == dashes:
                parse_lines = False if parse_lines else True
                continue
            if parse_lines:
                line_split = line.strip().split("\t")
                
                start, end = [int(i) for i in line_split[0].replace(" :   ", "").split("  ->  ")]
                unit = line_split[-1].split(" ")[0] # assume resolution parameter was not set for mreps, i.e., all repeats are perfect                
                period = int(line_split[2].replace(" <", "").replace("> ", ""))
                ref = int(float(line_split[3].replace(" [", "").replace("] ", ""))) # exp will never be negative so int() == math.floor()

                df_repeats["str_id"].append(f"{seqname}_{start}")
                df_repeats["chr"].append(seqname)                
                df_repeats["start"].append(start)
                df_repeats["end"].append(end)
                df_repeats["unit"].append(unit)
                df_repeats["period"].append(period)
                df_repeats["ref"].append(ref)

    df_repeats = pd.DataFrame(df_repeats)
    return df_repeats

if __name__ == "__main__":
    df_repeats = load_mreps_repeats("data/repeats/chr1_mreps.txt", "chr1")
    print(df_repeats.head())
    print(df_repeats.shape)
