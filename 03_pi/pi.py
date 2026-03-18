#!/usr/bin/env python3

import csv
from collections import defaultdict

# === CONFIGURATION ===
input_file = "b.afs/af.csv"

genome_lengths = {
    "CP006059.2": 113350,
    "CP006060.2": 190804,
    "NC_034781.1": 16863  # Mitochondria
}

# For both sliding and non-overlapping
chrom_window_steps = {
    "CP006059.2": (1000, 200),
    "CP006060.2": (1000, 200),
    "NC_034781.1": (100, 20)
}

include_chroms = set(genome_lengths.keys())

# === Step 1: Load allele frequency data ===
print("[*] Loading allele frequency data...")

site_alleles = defaultdict(list)

with open(input_file, "r", encoding="utf-8-sig") as f:
    reader = csv.reader(f, delimiter=",")
    header = [col.strip() for col in next(reader)]

    species_idx = header.index("Species")
    chrom_idx = header.index("CHROM")
    pos_idx = header.index("POS")
    dp_idx = header.index("DP")
    af_idx = header.index("AF")
    line_idx = header.index("Line")
    sex_idx = header.index("Sex")

    for row in reader:
        chrom = row[chrom_idx]
        if chrom not in include_chroms:
            continue

        af = float(row[af_idx])
        if af < 0.01:
            continue

        line = row[line_idx]
        species = row[species_idx]
        pos = int(row[pos_idx])
        dp = int(row[dp_idx])
        sex = row[sex_idx]

        key = (line, species, sex, chrom)
        site_alleles[key].append((pos, af, dp))

# === Step 2: Compute per-site π ===
print("[*] Calculating per-site π...")

site_data = []

with open("pi_site_1kb.tsv", "w", newline="") as f_site:
    writer = csv.writer(f_site, delimiter="\t")
    writer.writerow(["Line", "Species", "Sex", "CHROM", "POS", "N_ALLELES", "PI", "MEAN_DP"])

    for (line, species, sex, chrom), records in sorted(site_alleles.items()):
        for pos, af, dp in sorted(records):
            mean_dp = dp
            ref_af = 1.0 - af
            all_afs = [ref_af, af]
            pi = 1.0 - sum(p ** 2 for p in all_afs)
            site_data.append((line, species, sex, chrom, pos, pi, mean_dp))
            writer.writerow([line, species, sex, chrom, pos, 2, round(pi, 6), round(mean_dp, 2)])

# === Step 3: Non-overlapping window π (with fill and variant count)
print("[*] Computing non-overlapping π...")

window_data = defaultdict(lambda: defaultdict(lambda: {"sum_pi": 0.0, "n_variants": 0}))

for (line, species, sex, chrom, pos, pi, _) in site_data:
    win = chrom_window_steps[chrom][0]
    start = (pos // win) * win
    key = (line, species, sex, chrom)
    window_data[key][start]["sum_pi"] += pi
    window_data[key][start]["n_variants"] += 1

with open("pi_1kb.tsv", "w", newline="") as f_out:
    writer = csv.writer(f_out, delimiter="\t")
    writer.writerow(["Line", "Species", "Sex", "CHROM", "START", "END", "WINDOW_SIZE",
                     "N_VARIANT", "TOTAL_PI", "PI_per_base"])

    for key in sorted(window_data.keys()):
        line, species, sex, chrom = key
        genome_size = genome_lengths[chrom]
        win = chrom_window_steps[chrom][0]

        for start in range(0, genome_size, win):
            data = window_data[key].get(start, {"sum_pi": 0.0, "n_variants": 0})
            total_pi = data["sum_pi"]
            n_var = data["n_variants"]
            pi_per_base = total_pi / win
            writer.writerow([
                line, species, sex, chrom, start, start + win, win,
                n_var, round(total_pi, 6), round(pi_per_base, 8)
            ])

# === Step 4: Sliding window π ===
print("[*] Computing sliding window π...")

with open("pi_1kb_sliding.tsv", "w", newline="") as f_out:
    writer = csv.writer(f_out, delimiter="\t")
    writer.writerow(["Line", "Species", "Sex", "CHROM", "START", "END", "WINDOW_SIZE",
                     "STEP_SIZE", "N_SITES", "TOTAL_PI", "PI_per_base"])

    for (line, species, sex, chrom), records in sorted(site_alleles.items()):
        records.sort()
        genome_size = genome_lengths[chrom]
        win, step = chrom_window_steps[chrom]

        for start in range(0, genome_size - win + 1, step):
            end = start + win
            pi_sum = 0
            count = 0
            for pos, af, _ in records:
                if start <= pos < end:
                    ref_af = 1.0 - af
                    pi = 1.0 - (ref_af ** 2 + af ** 2)
                    pi_sum += pi
                    count += 1
            pi_per_base = pi_sum / win
            writer.writerow([line, species, sex, chrom, start, end, win, step, count,
                             round(pi_sum, 6), round(pi_per_base, 8)])

# === Step 5: Genome-wide π ===
print("[*] Calculating genome-wide π...")

genome_pi = defaultdict(lambda: {"sum_pi": 0, "variant_sites": 0})

for (line, species, sex, chrom, pos, pi, _) in site_data:
    key = (line, species, sex, chrom)
    genome_pi[key]["sum_pi"] += pi
    genome_pi[key]["variant_sites"] += 1

with open("pi_gw.tsv", "w", newline="") as f_out:
    writer = csv.writer(f_out, delimiter="\t")
    writer.writerow(["Line", "Species", "Sex", "CHROM", "N_SITES", "N_VARIANT", "TOTAL_PI", "PI_per_site"])

    for (line, species, sex, chrom), values in sorted(genome_pi.items()):
        total_sites = genome_lengths[chrom]
        pi_per_site = values["sum_pi"] / total_sites if total_sites > 0 else 0
        writer.writerow([
            line, species, sex, chrom,
            total_sites, values["variant_sites"],
            round(values["sum_pi"], 6), round(pi_per_site, 8)
        ])

print("✅ π calculated at site, window (both types), and genome-wide levels.")
