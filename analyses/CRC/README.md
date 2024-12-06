# Analysing STRs in colorectal cancer (CRC) sequencing data

Data used here are part of the EGAD50000000411 dataset in the human genome-phenome archive ([https://ega-archive.org/datasets/EGAD50000000411](https://ega-archive.org/datasets/EGAD50000000411)).

The dataset is described in this [pre-print](https://www.biorxiv.org/content/10.1101/2024.02.26.582054v1.full).
The authors of the pre-print ran ConSTRain on their aligned WGS data and provided us with four VCF files: LMO (bulk MSI organoid before any cloning), 01-0 (4n clone), 05-0 (2n clone), and 07-0 (4n clone).
It looks like ConSTRain was run with CNV information for each of these samples.

## Determining normalised depth filtering parameters

The script ["depth_overview.py"](../../bin/depth_overview.py) was used the generate plots of normalised sequencing depth from hte VCF files.
E.g.:

```bash
python3 depth_overview.py \
    -v CRC0282-07-0_constrain.vcf.gz \
    -o CRC0282-07-0_constrain.png \
    -a 0.05
```

Based on these plots, the minimal and maximal normalised depth values were determined for each plot, which are used as filter parameters in <jupyter notebook here>.