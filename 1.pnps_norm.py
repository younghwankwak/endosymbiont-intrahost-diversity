#!/usr/bin/env python3
"""
Normalize SNP variants for pN/pS calculation using allele frequency (AF)
and codon mutational opportunity.

Input:
- SnpEff-annotated VCF files
- Codon table TSV files (gene, codon_position, ref_codon)

Output:
- Per-sample normalized variant tables

Author: Y. Kwak
"""

import csv
import gzip
import os
from glob import glob
from collections import defaultdict
from Bio.Data import CodonTable

# =========================
# User-configurable paths
# =========================
VCF_DIR = "./vcf"
CODON_TABLES = {
    "NAS": "./codon_tables/NAS_codon_table.tsv",
    "SUL": "./codon_tables/SUL_codon_table.tsv",
}
OUT_DIR = "./pnps_normalized"
os.makedirs(OUT_DIR, exist_ok=True)

# =========================
# Variant category definitions
# =========================
SYN_VARIANTS = {"synonymous_variant"}
NONSYN_VARIANTS = {"missense_variant", "stop_gained"}

# =========================
# Genetic code
# =========================
genetic_code = CodonTable.unambiguous_dna_by_id[11].forward_table
stop_codons = CodonTable.unambiguous_dna_by_id[11].stop_codons
valid_codons = set(genetic_code.keys()) | set(stop_codons)


def codon_syn_nonsyn_counts(ref_codon):
    """Count possible synonymous and nonsynonymous single-base substitutions."""
    if ref_codon not in valid_codons:
        return 0, 0

    syn, nonsyn = 0, 0
    ref_aa = genetic_code.get(ref_codon, "*")

    for i in range(3):
        for base in "ACGT":
            if base == ref_codon[i]:
                continue
            alt = ref_codon[:i] + base + ref_codon[i + 1 :]
            if alt not in valid_codons:
                continue
            alt_aa = genetic_code.get(alt, "*")
            if alt_aa == ref_aa:
                syn += 1
            else:
                nonsyn += 1

    return syn, nonsyn


def load_codon_table(path):
    """Load codon table into gene → codon_position → ref_codon."""
    lookup = defaultdict(dict)
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            try:
                lookup[row["gene"]][int(row["codon_position"])] = row["ref_codon"].upper()
            except Exception:
                continue
    return lookup


def process_vcf(vcf_file, codon_table, output_path):
    """Parse one annotated VCF and write normalized pN/pS components."""
    rows = []
    opener = gzip.open if vcf_file.endswith(".gz") else open

    with opener(vcf_file, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue

            parts = line.strip().split("\t")
            if len(parts) < 10:
                continue

            info = parts[7]
            sample_field = parts[9]

            try:
                af = float(sample_field.split(":")[2])
            except Exception:
                continue

            ann = next((x for x in info.split(";") if x.startswith("ANN=")), None)
            if ann is None:
                continue

            for entry in ann.replace("ANN=", "").split(","):
                fields = entry.split("|")
                if len(fields) < 11:
                    continue

                effects = fields[1].split("&")
                gene = fields[3]
                aa_change = fields[10]

                if not aa_change.startswith("p."):
                    continue

                try:
                    codon_pos = int("".join(filter(str.isdigit, aa_change)))
                except ValueError:
                    continue

                ref_codon = codon_table.get(gene, {}).get(codon_pos)
                if not ref_codon:
                    continue

                syn_sites, nonsyn_sites = codon_syn_nonsyn_counts(ref_codon)

                for effect in effects:
                    if effect in SYN_VARIANTS and syn_sites:
                        weight = 3 / syn_sites
                        vtype = "synonymous"
                    elif effect in NONSYN_VARIANTS and nonsyn_sites:
                        weight = 3 / nonsyn_sites
                        vtype = "nonsynonymous"
                    else:
                        continue

                    rows.append({
                        "gene": gene,
                        "codon_pos": codon_pos,
                        "ref_codon": ref_codon,
                        "variant_type": vtype,
                        "AF": af,
                        "syn_sites": syn_sites,
                        "nonsyn_sites": nonsyn_sites,
                        "weight": round(weight, 6),
                        "weighted_value": round(af * weight, 6),
                        "effect": effect
                    })

    if not rows:
        return

    with open(output_path, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    for species, codon_path in CODON_TABLES.items():
        codon_table = load_codon_table(codon_path)
        vcfs = glob(os.path.join(VCF_DIR, species, "*.vcf.gz"))

        for vcf in vcfs:
            sample = os.path.basename(vcf).replace(".vcf.gz", "")
            out = os.path.join(OUT_DIR, f"{species}_{sample}_normalized.tsv")
            process_vcf(vcf, codon_table, out)
