#!/usr/bin/env python3
import numpy as np
import pandas as pd

from cyvcf2 import VCF

def dfs_from_vcf(file: str, samples: list) -> pd.DataFrame:    
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
    
    vcf = VCF(file, samples=samples)
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
            except TypeError:
                df_samples["copy_number"].append(np.nan)
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