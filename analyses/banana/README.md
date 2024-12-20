# Assessing ConSTRain's ability to analyse STRs in a polyploid banana (*Musa acuminata* Dwarf Cavendish)
To investigate whether ConSTRain is able to call STRs in short read sequencing data of a polyploid organism, we downloaded short sequencing reads of a banana sample.
The sample we analysed was taken from a *Musa acuminata* Dwarf Cavendish banana, which is an important commercial cultivar.
Sequencing reads were obtained from [Busche et al. (2020)](https://academic.oup.com/g3journal/article/10/1/37/6020287), and mapped to the DH-Pahang v4 reference genome ([Liu et al. 2023](https://www.nature.com/articles/s41597-023-02546-9)).

## Downloading the reference genome
The DH-Pahang v4 reference published by [Liu et al. 2023](https://www.nature.com/articles/s41597-023-02546-9) was downloaded from [Genoscope](https://www.genoscope.cns.fr/externe/plants/index.html).
```bash
curl -O https://www.genoscope.cns.fr/externe/plants/data/Musa_acuminata_pahang_v4.fasta
samtools faidx Musa_acuminata_pahang_v4.fasta
```

## Annotating STRs in the reference genome
We first split out chromosomes into individual files, and removed the putative mitochondrial sequence.
Taken from [stack overflow](https://stackoverflow.com/questions/21476033/splitting-a-multiple-fasta-file-into-separate-files-keeping-their-original-names)

```bash
awk 'BEGIN{RS=">";FS="\n"} NR>1{fnme=$1".fasta"; print ">" $0 > fnme; close(fnme);}' Musa_acuminata_pahang_v4.fasta
rm putative_mitochondrion.fasta
```

Then, we run mreps on each chromsome individually.

```bash
for CHR in {01..11};
do
    mreps -allowsmall -exp 3.0 -minperiod 1 -maxperiod 6 -fasta chr${CHR}.fasta > chr${CHR}_mreps.out;
done
```

By running mreps like this, we get MANY repeats.
We need to do some filtering to cut this set of repeats down to the more relevant STRs with allele length thresholds.
Also, mreps replaces 'N's in the reference sequence with a random nucleotide, so there may be spurious hits.
Thus, we also checked for each repeat whether the sequence reported by mreps exactly matches the reference genome in that position.
These steps are implemented in [combine_filter_mreps.ipynb](combine_filter_mreps.ipynb).
There are 183345 STRs left after combining and filtering.
These loci were written to a file to be used by ConSTRain later on.

## Downloading sequencing reads
[Busche et al. (2020)](https://academic.oup.com/g3journal/article/10/1/37/6020287) released their sequencing reads to the ENA under study accession [PRJEB33317](https://www.ebi.ac.uk/ena/browser/view/PRJEB33317).
We download them using a bash script generated by the site:
```bash
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR341/001/ERR3413471/ERR3413471_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR341/003/ERR3413473/ERR3413473_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR341/004/ERR3412984/ERR3412984_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR341/002/ERR3413472/ERR3413472_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR341/003/ERR3412983/ERR3412983_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR341/003/ERR3412983/ERR3412983_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR341/001/ERR3413471/ERR3413471_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR341/003/ERR3413473/ERR3413473_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR341/004/ERR3413474/ERR3413474_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR341/004/ERR3413474/ERR3413474_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR341/004/ERR3412984/ERR3412984_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR341/002/ERR3413472/ERR3413472_1.fastq.gz
```
The runs ERR3412983 and ERR3412984 are runs from a HiSeq1500 experiment.
the runs ERR3413471-4 are from a NextSeq500 experiment of the same sample.
The runs were merged (SAMEA5752290 is the sample accession in the ENA).
```bash
cat ERR341298*_1* > SAMEA5752290_HiSeq1500_1.fastq.gz
cat ERR341298*_2* > SAMEA5752290_HiSeq1500_2.fastq.gz

cat ERR341347*_1* > SAMEA5752290_NextSeq500_1.fastq.gz
cat ERR341347*_2* > SAMEA5752290_NextSeq500_2.fastq.gz
```

## Mapping reads to the reference genome
We used minimap2 to map the short sequencing reads to DH-Pahang v4 reference.
We remove improperly paired and duplicate alignments, and require at least Q30.
```bash
minimap2 \
    -ax sr -t 8 \
    Musa_acuminata_pahang_v4.fasta \
    SAMEA5752290_HiSeq1500_1.fastq.gz \
    SAMEA5752290_HiSeq1500_2.fastq.gz | \
    samtools view -@ 8 -F 2048 -f 2 -q 30 -b | \
    samtools sort -@ 8 > \
    SAMEA5752290_HiSeq1500.bam

samtools index SAMEA5752290_HiSeq1500.bam
```

Alignments for the NextSeq500 and HiSeq1500 experiments were merged to create a higher-depth alignment as well.
```bash
samtools cat -@ 8 SAMEA5752290_HiSeq1500.bam SAMEA5752290_NextSeq500.bam | \
    samtools sort -@ 8 > \
    SAMEA5752290_merged.bam

samtools index SAMEA5752290_merged.bam
```

## Genotyping STRs (with and without CNA information)
We then run ConSTRain on the resulting alignments.
Initially, without CNA information.

```bash
RUST_LOG=debug ConSTRain alignment \
    --repeats dh_pahang_v4_strs.bed \
    --alignment SAMEA5752290_NextSeq500.bam \
    --karyotype m_acuminata.json \
    --sample SAMEA5752290 \
    --threads 16 > \
    SAMEA5752290_NextSeq500.vcf
```

Then, we create updated VCFs by running ConSTRain on the previously generated VCFs but now with CNA info.

```bash
RUST_LOG=debug ConSTRain vcf \
    -k m_acuminata.json \
    -s SAMEA5752290 \
    -v SAMEA5752290_merged.vcf \
    --cnvs chr02_duplication.bed > \
    SAMEA5752290_merged_chr2_dup.vcf
```