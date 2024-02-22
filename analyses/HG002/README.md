# Assesing genotyper accuracy using HG002 assembly

Idea: compare estimated STR lengths to those in the HG002 assembly to determine how accurate our tool is. The HG002 assembly is based on multiple orthogonal sequencing approaches and is expected to be of very high quality, it will thus function as a gold standard to compare our STR calls against.

**downloading HG002 haplotypes**
```
aws s3 cp --no-sign-request s3://human-pangenomics/working/HPRC_PLUS/HG002/assemblies/year1_f1_assembly_v2_genbank/HG002.maternal.f1_assembly_v2_genbank.fa.gz ./

aws s3 cp --no-sign-request s3://human-pangenomics/working/HPRC_PLUS/HG002/assemblies/year1_f1_assembly_v2_genbank/HG002.paternal.f1_assembly_v2_genbank.fa.gz ./
```