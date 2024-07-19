# Running ExpansionHunter

Our STR panel had to be filtered for it to be usable by ExpansionHunter, as EH does not allow for loci that contain more than 5 'N' characters in their flanking regions (which is 1000 bp by default). Interestingly, ExpansionHunter also claimed to find too many 'N's for STR loci where the flank regions were fully lower-case, so these were removed as well (not sure if this is intended behaviour for EH...). This filtering step is implemented in [filter_flanking_N.py](filter_flanking_N.py), and reduced the number of loci from 1733646 to 1659608, so 74'038 loci had to be removed. Subsequently, [convert_str_panel_to_eh.py](convert_str_panel_to_eh.py) was used to convert the filtered STR panel to a json file in the format that ExpansionHunter expects (see documentation [here](https://github.com/Illumina/ExpansionHunter/blob/master/docs/04_VariantCatalogFiles.md)).

```bash
# File contains 1659608 entries total
split -l 829804 hg38_ver13_0boe_mononucleotides_EH.bed
```

```bash
python3 analyses/ExpansionHunter/convert_str_panel_to_eh.py -i data/tmp/hg38_ver13_EH_part1.bed -o data/tmp/hg38_ver13_EH_part1.json
python3 analyses/ExpansionHunter/convert_str_panel_to_eh.py -i data/tmp/hg38_ver13_EH_part2.bed -o data/tmp/hg38_ver13_EH_part2.json
```

ExpansionHunter was run using to following command line arguments:
```bash
# Tried to run with --analysis-mode streaming, but ran out of memory even with 180GB RAM on cluster
ExpansionHunter \
    --threads 32 \
    --analysis-mode seeking \
    --reads HG002.GRCh38.2x250.cram \
    --sex female \
    --reference hg38.fa \
    --variant-catalog hg38_ver13_mononucleotides_EH.json \
    --output-prefix HG002.GRCh38.2x250_EH
```