#!/usr/bin/env python3
import csv
from collections import defaultdict, Counter
from pathlib import Path

IN_CSV  = "af.csv"
OUT_TSV = "af_bins.tsv"

# ---- Configurable AF binning ----
LOW_MAX = 0.25
HIGH_MIN = 0.75
# Everything between (LOW_MAX, HIGH_MIN) is "medium"

def categorize_af(af: float) -> str:
    if af <= LOW_MAX:
        return "low"
    elif af < HIGH_MIN:
        return "medium"
    else:
        return "high"

def safe_float(x, default=None):
    try:
        return float(x)
    except Exception:
        return default

def main():
    if not Path(IN_CSV).exists():
        raise SystemExit(f"Input file not found: {IN_CSV}")

    bin_counts = defaultdict(Counter)

    with open(IN_CSV, newline="") as f:
        rdr = csv.DictReader(f)
        rdr.fieldnames = [h.strip().lstrip("\ufeff") for h in rdr.fieldnames]

        for col in ("Species", "Sample", "id", "AF"):
            if col not in rdr.fieldnames:
                raise SystemExit(f"Missing required column '{col}' in {IN_CSV} (found {rdr.fieldnames})")

        for row in rdr:
            af = safe_float(row["AF"])
            if af is None:
                continue

            af_bin = categorize_af(af)
            key = (row["Species"], row["Sample"], row["id"])
            bin_counts[key][af_bin] += 1

    # Build summary
    summary_rows = []
    for (species, sample, id), ctr in sorted(bin_counts.items()):
        low = ctr.get("low", 0)
        med = ctr.get("medium", 0)
        high = ctr.get("high", 0)
        total = low + med + high

        if total > 0:
            pct_low  = round(low  / total * 100.0, 2)
            pct_med  = round(med  / total * 100.0, 2)
            pct_high = round(high / total * 100.0, 2)
        else:
            pct_low = pct_med = pct_high = 0.0

        summary_rows.append({
            "Species": species,
            "Sample": sample,
            "id": id,
            "% low": f"{pct_low:.2f}",
            "% medium": f"{pct_med:.2f}",
            "% high": f"{pct_high:.2f}",
        })

    Path(OUT_TSV).parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_TSV, "w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["Species", "Sample", "id", "% low", "% medium", "% high"],
            delimiter="\t"
        )
        w.writeheader()
        w.writerows(summary_rows)

    print(f"Summary saved to {OUT_TSV} ✅")

if __name__ == "__main__":
    main()
