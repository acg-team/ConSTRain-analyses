#!/usr/bin/env python3
import numpy as np
import pandas as pd

def range_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start) + 1)

def expand_allele_lengths(lengths: dict):
    expanded = []
    for k, v in lengths.items():
        expanded.extend([k] * int(v))
    return np.array(sorted(expanded))
    
def str_len_from_haplotype(locus: dict, df_haplotype: pd.DataFrame) -> int:
    """
    IMPORTANT: df_haplotype must be pre-filtered to only contain indel variants 
    from the same chromosome as locus! This is assumed here for performance reasons.
    """
    str_region_length = locus["end"] - locus["start"] + 1

    indels = df_haplotype.query(f"start <= {locus['end']} and {locus['start'] - 1} <= end")
    for indel in indels.to_dict(orient="records"):
        if indel["REF"] == "-":
            # This is an insertion            
            str_region_length += len(indel["ALT"])
        elif indel["ALT"] == "-":
            # This is a deletion
            overlap = range_overlap(locus["start"], locus["end"], indel["start"] + 1, indel["end"])
            str_region_length -= overlap
        else: 
            raise ValueError(f"Non-indel variant encountered:{indel['chr']}\t{indel['start']}")
            
    return str_region_length
