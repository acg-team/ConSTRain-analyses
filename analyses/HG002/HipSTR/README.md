# Running HipSTR

Our STR panel had to be slightly modified for it to be usable by HipSTR (documentation [here](https://github.com/HipSTR-Tool/HipSTR)). This was done using [convert_str_panel_to_hipstr.py](convert_str_panel_to_hipstr.py). Since HipSTR does not allow for repeat regions with longer unit sizes, we remove all repeats where the period was longer than 6. This reduced the panel from 1733646 to 1713165 loci.

*Note: For installing on the HPC, I had to modify the Makefiles of both HipSTR and the version of htslib that is included in the HipSTR repo by adding -L/path/to/some/lib and -idirafter/path/to/some/include at the appropriate locations to correctly include and link against bzlib2. The resulting binary then also needs to be run with LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/path/to/some/lib*

HipSTR was run using to following command line arguments:
```bash
HipSTR \
    --bams HG002.GRCh38.2x250.cram \
    --fasta data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa \
    --regions data/hipstr_test_regions.tsv \
    --str-vcf hipstr_test.vcf.gz \
    --def-stutter-model \
    --min-reads 0
```