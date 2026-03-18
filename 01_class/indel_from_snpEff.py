#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import glob
import csv
from pathlib import Path
from collections import defaultdict

# ---- where to look ----
# Adjust ROOT if you run the script from a different working dir
ROOT = "./"
PATTERNS = [
    f"{ROOT}/1.NAS-ALF/6.annotated/2.pop/1.snpEff/*_snpEff_summary.html",
    f"{ROOT}/2.SUL-ALF/6.annotated/2.pop/1.snpEff/*_snpEff_summary.html",
    f"{ROOT}/3.MITO-ALF/6.annotated/2.pop/1.snpEff/*_snpEff_summary.html",
]

# ---- helpers ----
SPECIES_BY_DIR = {
    "1.NAS-ALF": "NAS",
    "2.SUL-ALF": "SUL",
    "3.MITO-ALF": "MITO",
}

# host ID: everything after the first underscore up to "_snpEff"
# e.g., NAS_42C_24H_5_snpEff_summary.html -> 42C_24H_5
HOST_ID_RE = re.compile(r"^[A-Za-z]+[_-]([^_]+(?:_[^_]+)*)_snpEff", re.I)

def get_species_from_path(path: Path) -> str:
    for key, lab in SPECIES_BY_DIR.items():
        if key in str(path):
            return lab
    # fallback: guess from basename prefix
    base = path.name
    if re.match(r"^NAS[_-]", base, re.I): return "NAS"
    if re.match(r"^SUL[_-]", base, re.I): return "SUL"
    if re.match(r"^MITO[_-]", base, re.I): return "MITO"
    return "UNK"

def host_id_from_basename(base: str) -> str:
    m = HOST_ID_RE.search(base)
    return m.group(1) if m else base

def parse_indel_histogram(html_text: str):
    """
    Parse the InDels section that looks like:

        <table class="histo">
          <tr><th>Min</th><td>0</td></tr>
          <tr><th>Max</th><td>50</td></tr>
          <tr><th>Mean</th><td>3.389</td></tr>
          <tr><th>Median</th><td>1</td></tr>
          <tr><th>Standard deviation</th><td>11.642</td></tr>
          <tr><th>Values</th><td>0,1,50</td></tr>
          <tr><th>Count</th><td>6,11,1</td></tr>

    Returns:
      summary: dict with min, max, mean, median, sd
      hist: list of (length:int, count:int)
    """
    # Narrow to the block after the InDels anchor (robust to whitespace)
    block = re.search(
        r"<a name=\"indels\".*?Insertions and deletions length:.*?<table class=\"histo\">(.*?)</table>",
        html_text, flags=re.I | re.S
    )
    if not block:
        return None, []

    tbl = block.group(1)

    def grab(label):
        m = re.search(rf"<th[^>]*>\s*{label}\s*</th>\s*<td>(.*?)</td>", tbl, flags=re.I | re.S)
        return m.group(1).strip() if m else None

    # Summary stats
    txt_min = grab("Min")
    txt_max = grab("Max")
    txt_mean = grab("Mean")
    txt_median = grab("Median")
    txt_sd = grab("Standard deviation")
    txt_vals = grab("Values")
    txt_cnts = grab("Count")

    summary = {
        "indel_min": float(txt_min) if txt_min is not None else None,
        "indel_max": float(txt_max) if txt_max is not None else None,
        "indel_mean": float(txt_mean) if txt_mean is not None else None,
        "indel_median": float(txt_median) if txt_median is not None else None,
        "indel_sd": float(txt_sd) if txt_sd is not None else None,
    }

    hist = []
    if txt_vals and txt_cnts:
        # e.g., Values: "0,1,50"  Count: "6,11,1"
        vals = [v.strip() for v in txt_vals.split(",") if v.strip() != ""]
        cnts = [c.strip() for c in txt_cnts.split(",") if c.strip() != ""]
        if len(vals) == len(cnts):
            for v, c in zip(vals, cnts):
                # SnpEff reports lengths as integers in this table
                try:
                    hist.append((int(float(v)), int(float(c))))
                except ValueError:
                    # skip any odd token
                    pass

    return summary, hist

def parse_variant_type_counts(html_text: str):
    """
    From the 'Number variants by type' table, retrieve counts for INS and DEL (and total InDels).
    """
    # capture the 'Number variants by type' table block
    block = re.search(
        r"<a name=\"changesByType\".*?<table.*?>(.*?)</table>",
        html_text, flags=re.I | re.S
    )
    counts = defaultdict(int)
    if not block:
        return counts

    rows = re.findall(r"<tr>\s*<td>\s*<b>\s*([A-Z]+)\s*</b>\s*</td>\s*<td.*?>\s*([0-9]+)\s*</td>\s*</tr>",
                      block.group(1), flags=re.I | re.S)
    for typ, num in rows:
        counts[typ.upper()] = int(num)

    # Combined indels if available
    counts["INDEL_TOTAL"] = counts.get("INS", 0) + counts.get("DEL", 0)
    return counts

# ---- main: collect files and parse ----
files = []
for pat in PATTERNS:
    files.extend(glob.glob(pat))

if not files:
    raise SystemExit("No SnpEff HTML files found. Check ROOT/PATTERNS.")

# outputs
summary_rows = []
hist_rows = []

for f in sorted(files):
    html = Path(f).read_text(errors="ignore")
    base = Path(f).name
    species = get_species_from_path(Path(f))
    host_id = host_id_from_basename(base)

    # per-sample variant type counts
    vt_counts = parse_variant_type_counts(html)

    # indel histogram & stats
    summary, hist = parse_indel_histogram(html)
    if summary is None:
        # still keep a row noting missing histogram (rare)
        summary = {"indel_min": None, "indel_max": None, "indel_mean": None, "indel_median": None, "indel_sd": None}
        hist = []

    # compute some handy derived metrics
    total_indels = vt_counts.get("INDEL_TOTAL", sum(c for _, c in hist))
    n_len1 = next((c for L, c in hist if L == 1), 0)
    frac_len1 = (n_len1 / total_indels) if total_indels else None

    # one wide summary row per sample
    summary_rows.append({
        "Sample": base.replace("_snpEff_summary.html", ""),
        "Species": species,
        "HostID": host_id,
        "INS_count": vt_counts.get("INS", None),
        "DEL_count": vt_counts.get("DEL", None),
        "INDEL_total": total_indels if total_indels else None,
        "Indel_min": summary["indel_min"],
        "Indel_max": summary["indel_max"],
        "Indel_mean": summary["indel_mean"],
        "Indel_median": summary["indel_median"],
        "Indel_sd": summary["indel_sd"],
        "Frac_len1_indels": frac_len1,
        "File": f,
    })

    # long-form histogram (one row per length)
    for L, C in hist:
        hist_rows.append({
            "Sample": base.replace("_snpEff_summary.html", ""),
            "Species": species,
            "HostID": host_id,
            "Indel_length_bp": L,
            "Count": C,
            "File": f,
        })

# write outputs
out1 = "indel_length_summary_per_sample.csv"
out2 = "indel_length_histogram_long.csv"

with open(out1, "w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=list(summary_rows[0].keys()))
    writer.writeheader()
    writer.writerows(summary_rows)

with open(out2, "w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=list(hist_rows[0].keys()) if hist_rows else
                            ["Sample","Species","HostID","Indel_length_bp","Count","File"])
    writer.writeheader()
    for r in hist_rows:
        writer.writerow(r)

print(f"Wrote:\n  {out1}\n  {out2}")
print("Columns in summary:",
      ", ".join(summary_rows[0].keys()))
