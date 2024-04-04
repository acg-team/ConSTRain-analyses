#!/usr/bin/env bash
set -euo pipefail

DATADIR="../../data/1000Genomes/cnvs_per_trio/"

for TRIODIR in $(ls ${DATADIR});
do
    sort -k1,1 -k2,2n ${DATADIR}/${TRIODIR}/*.bed | bedtools merge > ${DATADIR}/${TRIODIR}/trio${TRIODIR}_cnv_regions.bed;
done