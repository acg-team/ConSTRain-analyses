#!/usr/bin/env python3
"""Remove STR loci that contain 'N' characters in their flanking regions (ExpansionHunter does not allow for this)."""
import pandas as pd
import pysam

def main():    
    names = ["chr", "start", "end", "period", "unit"]
    df_strs = pd.read_csv("../../data/hg38_ver13_0boe_mononucleotides.bed", names=names, sep="\t")

    fasta = pysam.FastaFile("../../data/reference/hg38.fa.gz")
    length_map = dict()
    for contig, length in zip(fasta.references, fasta.lengths):
        length_map[contig] = length

    flanksize = 1000
    mask = []
    nucleotides = ['A', 'C', 'G', 'T']
    blacklist = [
        ("chr14", 16142377),
        ("chr14", 16154135),
        ("chr14", 16158054),        
    ]
    for idx, record in df_strs.iterrows():
        # hard coding some weird regions without N or n in flank, that EH still gets
        # mad about somehow (all lower case region?)
        # skip = False
        # for bl in blacklist:
        #     if bl[0] == record["chr"] and bl[1] == record["start"]:
        #         skip = True
        # if skip:
        #     continue
        chrom = record["chr"]
        start = max(0, record["start"] - flanksize)
        end = min(length_map[chrom], record["end"] + flanksize)

        region = fasta.fetch(chrom, start, end)
        characters = set(region)
        if "N" in characters or "n" in characters:
             # just remove everything that has an 'N'/'n' in it, even though EH allows for up to 5 'N's
            continue
        elif not any([c in nucleotides for c in characters]):
            # I suspect that EH also gets mad if the region consists entirely of lower case characters
            continue

        mask.append(idx)

    df_strs_filt = df_strs.loc[mask, :].reset_index(drop=True)
    df_strs_filt.to_csv("../../data/hg38_ver13_0boe_mononucleotides_EH.bed", index=False, header=None, sep="\t")




if __name__ == "__main__":
    main()
