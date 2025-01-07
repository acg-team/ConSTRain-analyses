# Annotating mononucleotide repeats in GRCh38
The GangSTR STR reference panel does not contain mononucleotide STRs.
Since we are also interested in investigating this type of STR, we extend the GangSTR panel by detecting mononucleotide STRs ourselves and merging the two panels.


## Splitting the GRCh38 fasta into individual chromosomes

Taken from [stack overflow](https://stackoverflow.com/questions/21476033/splitting-a-multiple-fasta-file-into-separate-files-keeping-their-original-names)
```bash
awk 'BEGIN{RS=">";FS="\n"} NR>1{fnme=$1".fasta"; print ">" $0 > fnme; close(fnme);}' hg38.fasta
```

## Removing non-standard chromosomes
```bash
find ./ -name "*_*.fasta" -exec rm {} \;

rm chrM.fasta
```

## Detecting mononucleotide repeats using mreps

[Paper](https://academic.oup.com/nar/article/31/13/3672/2904239?login=true)

[GitHub](https://github.com/gregorykucherov/mreps)

**ATTENTION:** mreps was modified before doing this to allow for a larger fraction of the input sequence to consist of 'N' characters. Specifically, line 149 of `defs.h` was changed: `#define LIMIT_N_PROPORTION 0.05` -> `#define LIMIT_N_PROPORTION 0.75`.

```bash
python3 generate_mreps_commands.py > mreps_commands.sh

bash mreps_commands.sh
```

## Check for spurious STRs introduced by MREPS
MREPS randomly replaces 'N' characters in the input sequence with A, C, G, or T.
This has a chance to introduce artificial repeats that do not exist in the reference genome. 
This is even more of a concern since we are interested in mononucleotide repeats. 
To address this, we need to check for every mreps mononucleotide whether the annotated sequences fully matches the specified region in the reference genome.

This is implemented in [this notebook](./include_mreps_repeats.ipynb).


