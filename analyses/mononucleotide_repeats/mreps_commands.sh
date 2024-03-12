#!/usr/bin/env bash
set -eou pipefail
echo 'annotating STRs on sequence: chr1'
mreps -exp 9.0 -minperiod 1 -maxperiod 1 -from 10000 -to 248946422 -fasta ../../data/reference/hg38_chromosomes/chr1.fasta > ../../data/repeats/chr1_mreps.txt
echo 'annotating STRs on sequence: chr10'
mreps -exp 9.0 -minperiod 1 -maxperiod 1 -from 10000 -to 133787422 -fasta ../../data/reference/hg38_chromosomes/chr10.fasta > ../../data/repeats/chr10_mreps.txt
echo 'annotating STRs on sequence: chr11'
mreps -exp 9.0 -minperiod 1 -maxperiod 1 -from 60000 -to 135076622 -fasta ../../data/reference/hg38_chromosomes/chr11.fasta > ../../data/repeats/chr11_mreps.txt
echo 'annotating STRs on sequence: chr12'
mreps -exp 9.0 -minperiod 1 -maxperiod 1 -from 10000 -to 133265309 -fasta ../../data/reference/hg38_chromosomes/chr12.fasta > ../../data/repeats/chr12_mreps.txt
echo 'annotating STRs on sequence: chr13'
mreps -exp 9.0 -minperiod 1 -maxperiod 1 -from 16000000 -to 114354328 -fasta ../../data/reference/hg38_chromosomes/chr13.fasta > ../../data/repeats/chr13_mreps.txt
echo 'annotating STRs on sequence: chr14'
mreps -exp 9.0 -minperiod 1 -maxperiod 1 -from 16000000 -to 106883718 -fasta ../../data/reference/hg38_chromosomes/chr14.fasta > ../../data/repeats/chr14_mreps.txt
echo 'annotating STRs on sequence: chr15'
mreps -exp 9.0 -minperiod 1 -maxperiod 1 -from 17000000 -to 101981189 -fasta ../../data/reference/hg38_chromosomes/chr15.fasta > ../../data/repeats/chr15_mreps.txt
echo 'annotating STRs on sequence: chr16'
mreps -exp 9.0 -minperiod 1 -maxperiod 1 -from 10000 -to 90228345 -fasta ../../data/reference/hg38_chromosomes/chr16.fasta > ../../data/repeats/chr16_mreps.txt
echo 'annotating STRs on sequence: chr17'
mreps -exp 9.0 -minperiod 1 -maxperiod 1 -from 60000 -to 83247441 -fasta ../../data/reference/hg38_chromosomes/chr17.fasta > ../../data/repeats/chr17_mreps.txt
echo 'annotating STRs on sequence: chr18'
mreps -exp 9.0 -minperiod 1 -maxperiod 1 -from 10000 -to 80263285 -fasta ../../data/reference/hg38_chromosomes/chr18.fasta > ../../data/repeats/chr18_mreps.txt
echo 'annotating STRs on sequence: chr19'
mreps -exp 9.0 -minperiod 1 -maxperiod 1 -from 60000 -to 58607616 -fasta ../../data/reference/hg38_chromosomes/chr19.fasta > ../../data/repeats/chr19_mreps.txt
echo 'annotating STRs on sequence: chr2'
mreps -exp 9.0 -minperiod 1 -maxperiod 1 -from 10000 -to 242183529 -fasta ../../data/reference/hg38_chromosomes/chr2.fasta > ../../data/repeats/chr2_mreps.txt
echo 'annotating STRs on sequence: chr20'
mreps -exp 9.0 -minperiod 1 -maxperiod 1 -from 60000 -to 64334167 -fasta ../../data/reference/hg38_chromosomes/chr20.fasta > ../../data/repeats/chr20_mreps.txt
echo 'annotating STRs on sequence: chr21'
mreps -exp 9.0 -minperiod 1 -maxperiod 1 -from 5010000 -to 46699983 -fasta ../../data/reference/hg38_chromosomes/chr21.fasta > ../../data/repeats/chr21_mreps.txt
echo 'annotating STRs on sequence: chr22'
mreps -exp 9.0 -minperiod 1 -maxperiod 1 -from 10510000 -to 50808468 -fasta ../../data/reference/hg38_chromosomes/chr22.fasta > ../../data/repeats/chr22_mreps.txt
echo 'annotating STRs on sequence: chr3'
mreps -exp 9.0 -minperiod 1 -maxperiod 1 -from 10000 -to 198235559 -fasta ../../data/reference/hg38_chromosomes/chr3.fasta > ../../data/repeats/chr3_mreps.txt
echo 'annotating STRs on sequence: chr4'
mreps -exp 9.0 -minperiod 1 -maxperiod 1 -from 10000 -to 190204555 -fasta ../../data/reference/hg38_chromosomes/chr4.fasta > ../../data/repeats/chr4_mreps.txt
echo 'annotating STRs on sequence: chr5'
mreps -exp 9.0 -minperiod 1 -maxperiod 1 -from 10000 -to 181478259 -fasta ../../data/reference/hg38_chromosomes/chr5.fasta > ../../data/repeats/chr5_mreps.txt
echo 'annotating STRs on sequence: chr6'
mreps -exp 9.0 -minperiod 1 -maxperiod 1 -from 60000 -to 170745979 -fasta ../../data/reference/hg38_chromosomes/chr6.fasta > ../../data/repeats/chr6_mreps.txt
echo 'annotating STRs on sequence: chr7'
mreps -exp 9.0 -minperiod 1 -maxperiod 1 -from 10000 -to 159335973 -fasta ../../data/reference/hg38_chromosomes/chr7.fasta > ../../data/repeats/chr7_mreps.txt
echo 'annotating STRs on sequence: chr8'
mreps -exp 9.0 -minperiod 1 -maxperiod 1 -from 60000 -to 145078636 -fasta ../../data/reference/hg38_chromosomes/chr8.fasta > ../../data/repeats/chr8_mreps.txt
echo 'annotating STRs on sequence: chr9'
mreps -exp 9.0 -minperiod 1 -maxperiod 1 -from 10000 -to 138334717 -fasta ../../data/reference/hg38_chromosomes/chr9.fasta > ../../data/repeats/chr9_mreps.txt
echo 'annotating STRs on sequence: chrX'
mreps -exp 9.0 -minperiod 1 -maxperiod 1 -from 10000 -to 156030895 -fasta ../../data/reference/hg38_chromosomes/chrX.fasta > ../../data/repeats/chrX_mreps.txt
echo 'annotating STRs on sequence: chrY'
mreps -exp 9.0 -minperiod 1 -maxperiod 1 -from 10000 -to 57217415 -fasta ../../data/reference/hg38_chromosomes/chrY.fasta > ../../data/repeats/chrY_mreps.txt
