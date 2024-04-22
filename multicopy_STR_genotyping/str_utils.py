#!/usr/bin/env python3
from itertools import combinations

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

def possible_chromosomes(genotype: np.array) -> set:
    genotype.sort()
    if genotype.shape[0] == 0:
        # when the CN is 0, we return two placeholder alleles -1
        return {(-1,), (-1,)}
    if genotype.shape[0] == 1:
        # when the CN is 1, we add in a placeholder allele -1        
        return {(-1,), tuple(genotype)}
        
    indexes = np.arange(0, genotype.shape[0])
    possibilities = set()
    for i in range(1, int(indexes.shape[0] / 2) + 1):
        idx_combos = combinations(indexes, i)
        for c in idx_combos:
            possibilities.add(tuple(genotype[c, ]))
            possibilities.add(tuple(np.delete(genotype, c)))

    return possibilities

def child_gt_possible(child: np.array, parent_1: np.array, parent_2: np.array) -> bool:        
    p1_chroms = possible_chromosomes(parent_1)
    p2_chroms = possible_chromosomes(parent_2)

    # Iterate over all possible combinations of two chromosomes based on estimated
    # genotype of child. For each combination, see if it could have originated 
    # from the genotypes estimated for the parents
    # Unfortunately, we repeat basically all of possible_chromosomes() here, as we have to
    # iterate over possible pairs of complementary child chromosomes that could occur together
    # based on the child genotype
    child.sort()
    if child.shape[0] == 0:
        # when the CN is 0, we add two placeholder alleles -1
        child = np.array([-1, -1])
    if child.shape[0] == 1:
        # when the CN is 1, we add in a placeholder allele -1
        child = np.array([-1, *child])

    indexes = np.arange(0, child.shape[0])
    for i in range(1, int(indexes.shape[0] / 2) + 1):
        idx_combos = combinations(indexes, i)
        for c in idx_combos:
            chrom1 = tuple(child[c, ])
            chrom2 = tuple(np.delete(child, c)) # complement of chrom1
            if chrom1 in p1_chroms and chrom2 in p2_chroms:
                return True
            if chrom1 in p2_chroms and chrom2 in p1_chroms:
                return True
    return False
