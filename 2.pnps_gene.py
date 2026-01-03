#!/usr/bin/env python3
"""
Aggregate normalized pN/pS components to the gene level
across samples and species.
"""

import pandas as pd
import glob
import os
import re

INPUT_DIR = "./pnps_normalized"
OUT_FILE = "pnps_gene_level.tsv"


def parse_meta(filename):
    base = os.path.basename(filename)
    m = re.match(r"([A-Z]+)_(.+?)_normalized\.tsv", base)
    return m.group(1), m.group(2)


dfs = []

for f in glob.glob(os.path.join(INPUT_DIR, "*_normalized.tsv")):
    species, sample = parse_meta(f)
    df = pd.read_csv(f, sep="\t")

    if df.empty:
        continue

    g = (
        df.groupby(["gene", "variant_type"])
          .agg(weighted_value=("weighted_value", "sum"),
               syn_sites=("syn_sites", "sum"),
               nonsyn_sites=("nonsyn_sites", "sum"))
          .reset_index()
    )

    piv = g.pivot(index="gene", columns="variant_type")
    piv.columns = ["_".join(c) for c in piv.columns]
    piv = piv.fillna(0)

    piv["pN"] = piv["weighted_value_nonsynonymous"] / piv["nonsyn_sites_nonsynonymous"].replace(0, pd.NA)
    piv["pS"] = piv["weighted_value_synonymous"] / piv["syn_sites_synonymous"].replace(0, pd.NA)

    piv["species"] = species
    piv["sample"] = sample

    dfs.append(piv.reset_index())

pd.concat(dfs, ignore_index=True).to_csv(OUT_FILE, sep="\t", index=False)
