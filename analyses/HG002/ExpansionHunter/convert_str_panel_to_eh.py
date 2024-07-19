#!/usr/bin/env python3
import argparse
import json

import pandas as pd

def parse_cli():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i", "--input", type=str, required=True, help="Path to bed file to convert to EH json"
    )
    parser.add_argument(
        "-o", "--output", type=str, required=True, help="Path EH json will be written to"
    )

    return parser.parse_args()

def main():
    args = parse_cli()

    names = ["chr", "start", "end", "period", "unit"]
    df_str_panel = pd.read_csv(args.input, sep="\t", names=names)
    
    df_str_panel = df_str_panel.assign(
        locus_id = lambda x: [f"{chrom}_{start}" for chrom, start in zip(x["chr"], x["start"])],
        locus_structure = lambda x: [f"({unit})*" for unit in x["unit"]],
        reference_region = lambda x: [f"{chrom}:{start}-{end}" for chrom, start, end in zip(x["chr"], x["start"], x["end"])]
    )
    
    json_out = []
    for _, row in df_str_panel.iterrows():
        record = {
            "LocusId": row["locus_id"],
            "LocusStructure": row["locus_structure"],
            "ReferenceRegion": row["reference_region"],
            "VariantType": "Repeat"
        }
        json_out.append(record)
    

    with open(args.output, 'w') as o:
        o.write(json.dumps(json_out, indent=4, sort_keys=True))

if __name__ == "__main__":
    main()
