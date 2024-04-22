#!/usr/bin/env python3
import numpy as np
import pandas as pd

from cyvcf2 import VCF

def dfs_from_vcf(filename: str, samples: list = None, format: str = None) -> pd.DataFrame:    
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
    
    vcf = VCF(filename, samples=samples)
    for variant in vcf:
        str_id = f"{variant.CHROM}_{variant.POS}"
        df_repeats["str_id"].append(str_id)
        df_repeats["chr"].append(variant.CHROM)
        df_repeats["start"].append(variant.POS)
        df_repeats["end"].append(variant.INFO.get("END"))
        df_repeats["unit"].append(variant.INFO.get("RU"))
        df_repeats["period"].append(variant.INFO.get("PERIOD"))
        df_repeats["ref"].append(int(variant.INFO.get("REF")))

        for sample_idx, sample in enumerate(vcf.samples):
            df_samples["sample"].append(sample)
            df_samples["str_id"].append(str_id)
            try:
                copy_number = variant.format("CN")[sample_idx]
                df_samples["copy_number"].append(copy_number[0])
            except KeyError:
                df_samples["copy_number"].append(np.nan)
            try:
                if format == "GangSTR":
                    frequencies = variant.format("ENCLREADS")[sample_idx]
                else:
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
                if format == "GangSTR":
                    if genotypes.size == 0:
                        raise TypeError
                else:
                    if genotypes == ".":
                        raise TypeError
                    genotypes = [int(i) for i in genotypes.split(",")]
                df_samples["genotype"].append(genotypes)
            except TypeError:
                df_samples["genotype"].append(np.nan)
            
    return pd.DataFrame(df_repeats), pd.DataFrame(df_samples)

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
