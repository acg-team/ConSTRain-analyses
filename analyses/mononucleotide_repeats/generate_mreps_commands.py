#!/usr/bin/env python3
from Bio import SeqIO

TARGET_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

def find_flanks(sequence):
    start, end = 0, 0
    for idx, nucleotide in enumerate(str(sequence.seq)):
        if nucleotide != "N":
            start = idx
            break
    for idx, nucleotide in enumerate(str(sequence.seq[::-1])):
        if nucleotide != "N":
            end = len(sequence) - idx
            break
    return start, end


def main():
    parser = SeqIO.parse("../../data/reference/hg38.fa", "fasta")
    print("#!/usr/bin/env bash")
    print("set -eou pipefail")

    for record in parser:
        if not record.id in TARGET_CHROMS:
            continue
        start, end = find_flanks(record)
        print(f"echo 'annotating STRs on sequence: {record.id}'")
        print(f"mreps -exp 9.0 -minperiod 1 -maxperiod 1 -from {start} -to {end} -fasta ../../data/reference/hg38_chromosomes/{record.id}.fasta > ../../data/repeats/{record.id}_mreps.txt")

if __name__ ==  "__main__":
    main()
