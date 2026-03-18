#!/usr/bin/env python3
"""
Compute AF-weighted disruptive, in-frame, and noncoding indel rates per gene
for ALL samples and emit one long-format table.

Long table columns:
sample, species, gene, gene_length_nt, disruptive_rate, inframe_rate, noncoding_rate

Only true indels included (len(REF) != len(ALT)).
Designed to mirror pnps.py logic but without codon weighting.

Author: Y. Kwak (2025)
"""

import csv
import gzip
import os
from glob import glob
from collections import defaultdict

# =============================================================================
# Config paths
# =============================================================================
VCF_DIR = "./"
CODON_TABLES = {
}

OUT_FILE = "g.indel/indel_rates_all_samples.tsv"

# =============================================================================
# Indel effect categories
# =============================================================================
DISRUPTIVE = {
    "frameshift_variant",
    "disruptive_inframe_insertion",
    "disruptive_inframe_deletion",
    "stop_gained", "stop_lost", "start_lost"
}

INFRAME = {
    "inframe_insertion", "inframe_deletion",
    "conservative_inframe_insertion",
    "conservative_inframe_deletion"
}

NONCODING = {
    "intergenic_region",
    "intragenic_variant",
    "feature_elongation"
}

# =============================================================================
# Load gene lengths using codon table
# =============================================================================
def load_gene_lengths(codon_table_path):
    """
    Returns dictionary: {gene : coding_length_nt}
    Based on max codon_position × 3.
    """
    gene_max_pos = defaultdict(int)

    with open(codon_table_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            try:
                gene = row["gene"]
                pos = int(row["codon_position"])
                if pos > gene_max_pos[gene]:
                    gene_max_pos[gene] = pos
            except:
                continue

    return {g: pos * 3 for g, pos in gene_max_pos.items()}


# =============================================================================
# Process a single VCF → returns gene-level AF sums
# =============================================================================
def process_vcf(vcf_file, gene_lengths):
    """
    Reads a normalized SnpEff-annotated VCF and returns:
    dicts: disruptive_sum[gene], inframe_sum[gene], noncoding_sum[gene]
    """
    open_fn = gzip.open if vcf_file.endswith(".gz") else open

    disruptive_sum = defaultdict(float)
    inframe_sum = defaultdict(float)
    noncoding_sum = defaultdict(float)

    with open_fn(vcf_file, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue

            parts = line.strip().split("\t")
            if len(parts) < 10:
                continue

            chrom, pos, _, ref, alt, _, _, info, fmt, sample_field = parts

            # Extract AF using FORMAT
            try:
                fmt_keys = fmt.split(":")
                vals = sample_field.split(":")
                fmt_map = dict(zip(fmt_keys, vals))
                af = float(fmt_map["AF"].split(",")[0])
            except:
                continue

            # Strictly identify indels
            if len(ref) == len(alt):
                continue

            # Extract ANN field
            anns = [x for x in info.split(";") if x.startswith("ANN=")]
            if not anns:
                continue

            ann_entries = anns[0].replace("ANN=", "").split(",")

            for ann in ann_entries:
                fields = ann.split("|")
                if len(fields) < 11:
                    continue

                allele = fields[0]
                if allele != alt:
                    continue

                effects = fields[1].split("&")
                gene = fields[3]

                if gene not in gene_lengths:
                    continue

                # Classify indel
                if any(e in DISRUPTIVE for e in effects):
                    disruptive_sum[gene] += af
                elif any(e in INFRAME for e in effects):
                    inframe_sum[gene] += af
                elif any(e in NONCODING for e in effects):
                    noncoding_sum[gene] += af

    return disruptive_sum, inframe_sum, noncoding_sum


# =============================================================================
# MAIN: Build long-format table
# =============================================================================
if __name__ == "__main__":
    rows = []

    with open(OUT_FILE, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow([
            "sample", "species", "gene", "gene_length_nt",
            "disruptive_rate", "inframe_rate", "noncoding_rate"
        ])

        for species in ["1.NAS-ALF", "2.SUL-ALF"]:
            tag = "NAS" if "NAS" in species else "SUL"
            gene_lengths = load_gene_lengths(CODON_TABLES[tag])

            vcf_path = f"{VCF_DIR}/{species}/6.annotated/2.pop/4.normalized/"
            vcfs = sorted(glob(os.path.join(vcf_path, "*.norm.vcf.gz")))

            if not vcfs:
                print(f"⚠ No VCFs found for {species} in {vcf_path}")
                continue

            for vcf in vcfs:
                sample = os.path.basename(vcf).replace(".norm.vcf.gz", "")

                # Parse VCF to get AF sums
                disruptive_sum, inframe_sum, noncoding_sum = process_vcf(vcf, gene_lengths)

                # Write rows for ALL genes (even genes without indels → rate = 0)
                for gene, L in gene_lengths.items():
                    dis = disruptive_sum[gene] / L if L > 0 else 0
                    inf = inframe_sum[gene] / L if L > 0 else 0
                    non = noncoding_sum[gene] / L if L > 0 else 0

                    writer.writerow([sample, tag, gene, L, dis, inf, non])

    print(f"✔ Complete long-format table written → {OUT_FILE}")
