# Assesing genotyper accuracy using HG002 assembly

Idea: compare estimated STR lengths to those in the HG002 assembly to determine how accurate our tool is. The HG002 assembly is based on multiple orthogonal sequencing approaches and is expected to be of very high quality, it will thus function as a gold standard to compare our STR calls against.

## Determining HG002 variants relative to GRCh38 from haplotypes

#### downloading HG002 haplotypes
```bash
aws s3 cp --no-sign-request s3://human-pangenomics/working/HPRC_PLUS/HG002/assemblies/year1_f1_assembly_v2_genbank/HG002.maternal.f1_assembly_v2_genbank.fa.gz ./

aws s3 cp --no-sign-request s3://human-pangenomics/working/HPRC_PLUS/HG002/assemblies/year1_f1_assembly_v2_genbank/HG002.paternal.f1_assembly_v2_genbank.fa.gz ./
```

#### mapping HG002 haplotypes to GRCh38 using minimap2

Explanation on the [minimap2 GitHub page](https://github.com/lh3/minimap2/tree/master/misc#asmvar)
```bash
# Maternal haplotype
minimap2 -cx asm5 -t 8 --cs hg38.fa.gz HG002.maternal.f1_assembly_v2_genbank.fa.gz > HG002.maternal.f1_assembly_v2_genbank.paf

sort -k6,6 -k8,8n HG002.maternal.f1_assembly_v2_genbank.paf > HG002.maternal.f1_assembly_v2_genbank_sort.paf
k8 paftools.js call HG002.maternal.f1_assembly_v2_genbank_sort.paf > HG002.maternal.f1_assembly_v2_genbank_sort_var.txt

# Paternal haplotype
minimap2 -cx asm5 -t 8 --cs hg38.fa.gz HG002.paternal.f1_assembly_v2_genbank.fa.gz > HG002.paternal.f1_assembly_v2_genbank.paf

sort -k6,6 -k8,8n HG002.paternal.f1_assembly_v2_genbank.paf > HG002.paternal.f1_assembly_v2_genbank_sort.paf
k8 paftools.js call HG002.paternal.f1_assembly_v2_genbank_sort.paf > HG002.paternal.f1_assembly_v2_genbank_sort_var.txt
```

## Calling STR variants in HG002 from Illumina short read sequencing

#### Downloading 2x250 Illumina reads of HG002 aligned to GRCh32
```bash
# 100x depth, filesize ~122GB
curl -L -O https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/novoalign_bams/HG002.GRCh38.2x250.bam
```

#### Downsampling to 30x and 10x
```bash
samtools view -@ 15 -T hg38.fa.gz -Cs 42.1 HG002.GRCh38.2x250.bam > HG002.GRCh38.2x250_depth10x.cram

# For 30X, could not convert directly to CRAM -> samtools: cram/cram_io.c:3180: cram_ref_decr_locked: Assertion `r->ref_id[id]->count == 0' failed.
samtools view -@ 15 -bs 42.3 HG002.GRCh38.2x250.bam > HG002.GRCh38.2x250_depth30x.bam
```

#### convert to CRAM
```bash
samtools view \
    -@ 15 \
    -T hg38.fa.gz \
    -C HG002.GRCh38.2x250.bam \
    chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY > \
    HG002.GRCh38.2x250.cram

samtools index -@ 15 HG002.GRCh38.2x250.cram
```

## STR length calling
### CN-guided STR genotyping
#### Get and modify GangSTR STR annotation of HG38
```bash
curl -L https://s3.amazonaws.com/gangstr/hg38/genomewide/hg38_ver13.bed.gz | \
    gunzip | \
    awk 'BEGIN {OFS = "\t" } {$2 -= 1} {print $0}' > \
    hg38_ver13.bed
```
#### Run genotyper
```bash
cn-guided-str-genotying \
    --reads-per-allele 0 \
    --threads 16 \
    --repeats hg38_ver13_0boe_mononucleotides.bed \ # modified from https://s3.amazonaws.com/gangstr/hg38/genomewide/hg38_ver13.bed.gz
    --ploidy h_sapiens_male.json \
    --sample HG002.GRCh38.2x250_depth10x \
    --alignment HG002.GRCh38.2x250_depth10x.cram \
    --reference hg38.fa.gz > \
    HG002.GRCh38.2x250_depth10x.vcf
```

### GangSTR

#### Run GangSTR
```bash
GangSTR \
    --bam HG002.GRCh38.2x250.cram \
    --ref hg38.fa \
    --regions hg38_ver13_1bce_mononucleotides.tsv \ # modified from https://s3.amazonaws.com/gangstr/hg38/genomewide/hg38_ver13.bed.gz
    --out HG002.GRCh38.2x250 \
    --bam-samps HG002.GRCh38.2x250 \
    --samp-sex M
```

## Determining which STR loci are covered by HG002 haplotypes
To properly determine the accuracy of our STR genotyping approach, we need to make sure to restrict our analysis
to loci that are covered by both the HG002 haplotypes *and* the Illumina short reads. 

We first extract the STR loci that are covered by the aligned HG002 genotypes using bedtools:
```bash
grep "^R" HG002.maternal.f1_assembly_v2_genbank_sort_var.txt | cut -f 2-4 > HG002_maternal_regions.txt
bedtools intersect -u -a hg38_ver13_0boe_mononucleotides.bed -b HG002_maternal_regions.txt > hg38_ver13_0boe_mononucleotides_mat.bed

grep "^R" HG002.paternal.f1_assembly_v2_genbank_sort_var.txt | cut -f 2-4 > HG002_paternal_regions.txt
bedtools intersect -u -a hg38_ver13_0boe_mononucleotides.bed -b HG002_paternal_regions.txt > hg38_ver13_0boe_mononucleotides_pat.bed

# switching -a and -b files gives exactly the same result (checked with openssl sha256)
bedtools intersect -u -a hg38_ver13_in_maternal.bed -b hg38_ver13_in_paternal.bed > hg38_ver13_union.bed

grep "^chrX" hg38_ver13_0boe_mononucleotides_mat.bed >> hg38_ver13_0boe_mononucleotides_union.bed # 81268 loci
grep "^chrY" hg38_ver13_0boe_mononucleotides_pat.bed >> hg38_ver13_0boe_mononucleotides_union.bed # 10549 loci
```
1695865 out of 1733646 STR loci (97.82%) are covered by the HG002 haplotypes 
