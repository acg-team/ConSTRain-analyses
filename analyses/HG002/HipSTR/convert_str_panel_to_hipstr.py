#!/usr/bin/env python3
import pandas as pd


def main():
    names = ["chr", "start", "end", "period", "unit"]
    df_str_panel = pd.read_csv("../../data/hg38_ver13_0boe_mononucleotides.bed", sep="\t", names=names)
    
    # Determine number of units in repeat region, increment start postion, remove repeats with unit size >6
    df_str_panel = (
        df_str_panel
            .query("period <= 6")
            .assign(n_units = lambda x: (x["end"] - x["start"]) / x["period"]))

    df_str_panel["start"] = df_str_panel["start"] + 1

    # rearrange columns to match HipSTR spec
    df_str_panel = df_str_panel[["chr", "start", "end", "period", "n_units",]]
    
    df_str_panel.to_csv("../../data/hg38_ver13_mononucleotides_hipstr.tsv", index=False, header=False, sep="\t")

if __name__ == "__main__":
    main()
