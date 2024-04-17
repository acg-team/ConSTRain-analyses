# Analysis of STRs in CNVs in 1000 Genomes triads
## Overview of samples of interest
Pedigree information for the 3202 samples included in the latest phase of the 1000 Genomes project is provided here: [https://www.internationalgenome.org/data-portal/data-collection/30x-grch38](https://www.internationalgenome.org/data-portal/data-collection/30x-grch38).
```bash
curl -L http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/1kGP.3202_samples.pedigree_info.txt > 1kGP.3202_samples.pedigree_info.txt
```

This file gets filtered and reformatted to end up with an overview of all the triads in the 1000 Genomes project in the [2024-04-02_reformat_1000G_pedigree.ipynb](2024-04-02_reformat_1000G_pedigree.ipynb) notebook. There are 602 triads, which corresponds to 1000 Genomes information online. Final format is like this:

|    |   triadID | sampleID   | role   | pgx_id          | sex   |
|---:|----------:|:-----------|:-------|:----------------|:------|
|  0 |         0 | HG00403    | father | onekgbs-HG00403 | M     |
|  1 |         0 | HG00404    | mother | onekgbs-HG00404 | F     |
|  2 |         0 | HG00405    | child  | onekgbs-HG00405 | F     |
|  3 |         1 | HG00406    | father | onekgbs-HG00406 | M     |
|  4 |         1 | HG00407    | mother | onekgbs-HG00407 | F     |

## Getting copy number calls from Progenetix
For this analysis, we copy the Progenetix database from the main server to the local computer. We connect to it using the `brew` installation of `MongoDB` and use `MongoDB Compass` to interact with the database and explore the schema. 1000 Genomes CN calls are in the `variants` collection, and can be accessed using the 1000 Genomes sample IDs, prefixed with `onekgbs-`.

We first generate a query usable by `mongoexport` from our previously generated pedigree file.
```bash
cut -d "," -f 4 1000_genomes_triad_pedigrees.csv | \ # select 4th column with Progenetix IDs
    tail -n +2 | \  # skip header line
    sort | \    # make IDs unique
    uniq | \
    awk '{print "\x22"$0"\x22"}' | \    # add double quotation marks
    tr "\n" "," | \ # replace linebreaks with commas
    pbcopy  # send to clipboard (macOS specific I guess)
```

We copy the resulting biosample IDs into the query in [pgx_mongodb_query.sh](pgx_mongodb_query.sh).

 **IMPORTANT:** the list of IDs cannot end with a comma, otherwise the mongoexport will not work. I manually removed the last comma from the list of IDs with a text editor.

 We can then export all variants belonging to out samples of interest using mongoexport with the following command:

```bash
mongoexport \
    -c variants \
    --type csv \
    --fields biosample_id,info.cn_count,location.chromosome,location.start,location.end \
    --query "`analyses/1000Genomes/pgx_mongodb_query.sh`" \     
    "mongodb://127.0.0.1:27017/progenetix?directConnection=true&serverSelectionTimeoutMS=2000&appName=mongosh+2.0.2" > \
    data/1000Genomes/2024-04-02_trio_cn_variants.csv
```

In total, this exported 1104925 records (every record represents one copy number variant). Next, we need to split this output file into files specific for each trio. First, we make a directory for each trio:

```bash
for i in $(seq 0 601); do mkdir $i; done
```

Then, we parse the `mongoexport` output file, generate bed files of CNV regions for each individual, and write them to the corresponding trio identifier directory that we just created. This happens in the notebook [2024-04-02_trio_specific_cnv_files.ipynb](2024-04-02_trio_specific_cnv_files.ipynb).

## Generating trio-specific STR panels to use for genotyping
For this analysis, we will investigate STRs that are located in CNVs. Since we have trios, we will analyze the STRs in each member of the trio that are affected by a CNV in *any* member of that trio. This may allow us to quantify Mendelian consistency of STR calls.

First, we merge the Progenetix bed files of each trio to get a single bed file per trio that covers all regions affected by a CNV in that trio. We use `bedtools merge` for this, see the bash script [merge_cnv_beds.sh](merge_cnv_beds.sh).

Then, we need to make a filtered STR panel for each trio that only contains loci located in a CNV in any member of that trio. We do this using `bedtools intersect` in the script [trio_specific_str_panels.sh](trio_specific_str_panels.sh).

## Calling STR lengths using ConSTRain

