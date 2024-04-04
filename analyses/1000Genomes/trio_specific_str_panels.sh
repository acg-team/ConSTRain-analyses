#!/usr/bin/env bash
set -euo pipefail

DATADIR="../../data/1000Genomes/cnvs_per_trio/"
PANEL="../../data/hg38_ver13_0boe_mononucleotides.bed"

for TRIODIR in $(ls ${DATADIR});
do
    bedtools intersect \
        -wa -f 1 \
        -a ${PANEL} \
        -b ${DATADIR}/${TRIODIR}/trio${TRIODIR}_cnv_regions.bed > \
        ${DATADIR}/${TRIODIR}/trio${TRIODIR}_strs_0boe.bed
done