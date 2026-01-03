#!/usr/bin/env python3
"""
Summarize gene-level pN/pS values at the COG level.
"""

import pandas as pd

PNPS_FILE = "pnps_gene_level.tsv"
COG_MAP_FILE = "gene_to_cog.tsv"
OUT_FILE = "pnps_cog_level.tsv"

pnps = pd.read_csv(PNPS_FILE, sep="\t")
cog = pd.read_csv(COG_MAP_FILE, sep="\t")

df = pnps.merge(cog, on=["gene", "species"], how="left")
df["COG"] = df["COG"].fillna("Unclassified")

agg = (
    df.groupby(["species", "sample", "COG"], as_index=False)
      .agg(pN=("pN", "sum"),
           pS=("pS", "sum"),
           n_genes=("gene", "count"))
)

agg.to_csv(OUT_FILE, index=False)
