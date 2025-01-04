# Simulating copy number variants (and more) to demonstrate ConSTRain
## Obtaining assemblies to simulate short reads from: first attempt
The [human pangenome reference consortium](https://humanpangenome.org/) released the full assemblies of the samples currently incorporated into the pangenome. Links two download the assemblies are indexed here: [https://github.com/human-pangenomics/HPP_Year1_Assemblies/blob/main/assembly_index/Year1_assemblies_v2_genbank.index](https://github.com/human-pangenomics/HPP_Year1_Assemblies/blob/main/assembly_index/Year1_assemblies_v2_genbank.index)


I'm planning to simulate trisomy 21, for which I will need assemblies from two individuals. Furthermore, I want to simulate a copy number variant (CNV), for which I'll need assemblies of just a single individual. Thus, I randomly sample three individuals from the assembly index file (excluding reference genomes and cell lines).
```python
import pandas as pd

df = pd.read_csv("https://raw.githubusercontent.com/human-pangenomics/HPP_Year1_Assemblies/main/assembly_index/Year1_assemblies_v2_genbank.index", sep="\t")
skip = ["HG002", "HG005", "CHM13_v1.1", "GRCh38_no_alt_analysis_set"]

df = df[~df["sample"].isin(skip)]

print(df["sample"].sample(n=3, random_state=1).values)
```

These are the sample ids that are returned (luckily, there are both male and female samples):
* HG00735 (female)
* HG00673 (male)
* HG01109 (male)

From the assembly index file, we can get the urls to download the maternal and paternal haplotypes of each individual.
```bash
# E.g., for HG00735
aws s3 cp --no-sign-request s3://human-pangenomics/working/HPRC/HG00735/assemblies/year1_f1_assembly_v2_genbank/HG00735.paternal.f1_assembly_v2_genbank.fa.gz ./
aws s3 cp --no-sign-request s3://human-pangenomics/working/HPRC/HG00735/assemblies/year1_f1_assembly_v2_genbank/HG00735.maternal.f1_assembly_v2_genbank.fa.gz ./

# For convenience, unzip and recompress with bgzip
gunzip HG00735.*
bgzip HG00735.paternal.f1_assembly_v2_genbank.fa
bgzip HG00735.maternal.f1_assembly_v2_genbank.fa
```

Upon closer inspection, these assemblies do not seem very practical to use for our purposes. The scaffolds contained in the haplotypes are not assembled into chromosomes. This means we would need to do work to construct chromosomes from these scaffolds, and then simulate reads from them. On the other hand, we could look for haplotypes that are assembled into chromosomes so we don't have to deal with this. One idea might be to use the GRCh38 and HG002 (and CHM13?) assemblies.

## Obtaining assemblies to simulate short reads from: second attempt
Extract chromosome 21 from GRCh38.
```bash
samtools faidx GRCh38_full_analysis_set_plus_decoy_hla.fa "chr21" | bgzip > GRCh38_chr21.fa.gz
samtools faidx GRCh38_chr21.fa.gz
```

All good so far, but HG002 has the same issues as the other assemblies... I guess we need te generate artificial chromosomes 21 from the STR variants called when we aligned the HG002 haplotypes to GRCh38 to determine ConSTRain's accuracy.

## Obtaining assemblies to simulate short reads from: third attempt
In the end, we just generate a new reference sequence ourselves that only differes from GRCh38 at STR loci that are not the same as the reference in the HG002 assemblies. This is done based on the STR variant calls from the assembly to reference alignment we generated previously (see [STR_lengths_from_haplotypes.ipynb](../HG002/STR_lengths_from_haplotypes.ipynb)). We generate a chr21 representation for the HG002 maternal and paternal haplotypes separately and write them to fasta file. See [reference_from_STR_genotypes.ipynb](reference_from_STR_genotypes.ipynb) for details.

We also bgzip and index the new fasta files.
```bash
bgzip HG002_chr21_maternal.fa
samtools faidx HG002_chr21_maternal.fa.gz
```

## Simulating short sequencing reads
We are not interested in generating realistic sequencing errors here, we just want to see how different variant callers handle CNVs/ aneuploidies. Therefore, we will use [wgsim](https://github.com/lh3/wgsim) to generate perfect reads (paired, 150bp) from the chr21 reference sequences we obtained. We will simulate 30X sequencing, so we will simulate reads to a depth of 15x for each haplotype. Wgsim takes as command-line argument `-N`, which is the number of reads to generate. To find the number of reads that will give us the desired sequencing depth, we use the following relationship ([Sims *et al.*, 2014](https://www.nature.com/articles/nrg3642)):
$$
\begin{equation}
    \begin{aligned}
        c &= \frac{LN}{G},
    \end{aligned}
\end{equation}
$$
where $c$ is the sequencing depth, $L$ is the read length, $N$ is the number of reads, and $G$ is the haploid genome length. We know that $G=46709983$ for chr21 in GRCh38. We plug this into the equation and solve for $N$:
$$
\begin{equation}  
    \begin{aligned}
        \frac{150N}{46709983} &= 15 \\[5pt]
        150N &= 15 * 46709983 \\[5pt]
        N &= \frac{15 * 46709983}{150} = 0.1 * 46709983 = 4670998.3
    \end{aligned}
\end{equation}
$$
Huh, because of how we chose our $L$ and $c$, we can just multiply $L$ by $0.1$ to get $N$. Finally, we need to divide this value by 2 before we give it to wgsim to account for the fact we're simulating paired reads. This means the `-N` value we give to wgsim will be 2335499, 2335790, and 2335625 for GRCh38, HG002 maternal, and HG002 paternal, respectively.

We simulate the reads... (example shown for GRCh38 chr21)
```bash
# wgsim takes around 41s, unzipped output files are 778M each, 144M each when gzipped
wgsim -e 0 -r 0 -R 0 -X 0 -S 42 -1 150 -2 150 -N 2335499 GRCh38_chr21.fa.gz GRCh38_chr21_1.fq GRCh38_chr21_2.fq
gzip GRCh38_chr21_*.fq
```

...and then map them to the GRCh38 reference representation of chr21
```bash
# whole pipeline takes 1m12s
minimap2 -ax sr haplotypes/GRCh38_chr21.fa.gz reads/GRCh38_chr21_1.fq.gz reads/GRCh38_chr21_2.fq.gz | \
    samtools view -b | \
    samtools sort -@ 3 > \
    alignments/GRCh38_chr21_GRCh38.bam

samtools index alignments/GRCh38_chr21_GRCh38.bam
```

And indeed, `samtools coverage` tells us that the average read depth in the alignments are 15, 14.9978, and 14.9981 for the GRCh38, maternal, and paternal haplotype reads. Finally, we merge all files into one big alignment to analyze later.
```bash
samtools cat GRCh38_chr21_GRCh38.bam HG002_chr21_maternal.bam HG002_chr21_paternal.bam | samtools sort -@ 3 > simulated_trisomy_21.bam

samtools index simulated_trisomy_21.bam
```
## Genotyping 
<!-- chr21_13941477, chr21_15583008, chr21_15678835 could be used as example loci -->

```bash
# Takes 42.8s to run on single core
ConSTRain \
    --repeats hg38_ver13_0boe_mononucleotides_chr21.bed \
    --karyotype h_sapiens_male_tri21.json \
    --sample simulated_trisomy_21 \
    --alignment simulated_trisomy_21.bam \
    --reference GRCh38_chr21.fa.gz > \
    HG002.GRCh38.2x250_depth10x.vcf
```

```bash
# Takes 7m44s to run (roughly 11x of ConSTRain on single core)
GangSTR \
    --bam simulated_trisomy_21.bam \
    --ref GRCh38_chr21.fa.gz \
    --regions hg38_ver13_1bce_mononucleotides_chr21.tsv \
    --out GangSTR_simulated_trisomy_21 \
    --bam-samps simulated_trisomy_21 \
    --samp-sex M
```

```bash
HipSTR \
   --bams simulated_trisomy_21.bam \
   --bam-samps tri21 \
   --bam-libs tri21 \
   --fasta GRCh38_chr21.fa.gz \
   --regions hg38_ver13_mononucleotides_hipstr_chr21.tsv \
   --def-stutter-model \
   --min-reads 0 \
   --str-vcf HipSTR_simulated_trisomy_21.vcf.gz
```