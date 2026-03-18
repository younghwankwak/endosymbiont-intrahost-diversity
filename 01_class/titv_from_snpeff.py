# titv_from_snpeff_summary_only_v2.py
import re, glob
from pathlib import Path
import pandas as pd
from bs4 import BeautifulSoup

ROOT = "root"
SUBDIR = "subdir"
PATTERNS = [
    f"{ROOT}/1.NAS-ALF/{SUBDIR}/*_snpEff_summary.html",
    f"{ROOT}/2.SUL-ALF/{SUBDIR}/*_snpEff_summary.html",
    f"{ROOT}/3.MITO-ALF/{SUBDIR}/*_snpEff_summary.html",
]

SPECIES_FROM_DIR = {
    "1.NAS-ALF": "NAS",
    "2.SUL-ALF": "SUL",
    "3.MITO-ALF": "MITO",
}

def infer_species_from_path(p: Path) -> str:
    for key, lab in SPECIES_FROM_DIR.items():
        if key in str(p):
            return lab
    # fallback to basename prefix (NAS_, SUL_, MITO_)
    b = p.name
    if re.match(r"^NAS[_-]", b, re.I): return "NAS"
    if re.match(r"^SUL[_-]", b, re.I): return "SUL"
    if re.match(r"^MITO[_-]", b, re.I): return "MITO"
    return "UNK"

def host_id_from_basename(basename: str) -> str:
    # NAS_25C_24H_1_F2_snpEff_summary.html -> 25C_24H_1_F2
    # Split at first "_", drop species token, cut trailing "_snpEff..."
    core = basename.split("_snpEff", 1)[0]
    parts = core.split("_", 1)
    return parts[1] if len(parts) > 1 else core

def parse_titv_with_bs4(html_path: str):
    """Parse Ts/Tv from the Ts/Tv table; robust fallbacks included."""
    txt = Path(html_path).read_text(errors="ignore")
    soup = BeautifulSoup(txt, "lxml")

    # 1) Try the Ts/Tv section anchored by <a name="tstv">
    titv = Ti = Tv = None
    anchor = soup.find("a", attrs={"name": "tstv"})
    if anchor:
        tbl = anchor.find_next("table")
        if tbl:
            rows = tbl.find_all("tr")
            for r in rows:
                cells = [c.get_text(strip=True) for c in r.find_all(["th", "td"])]
                if not cells:
                    continue
                label = cells[0].lower()
                if "transition" in label and len(cells) > 1:
                    try: Ti = int(cells[1].replace(",", ""))
                    except: pass
                if "transversion" in label and len(cells) > 1:
                    try: Tv = int(cells[1].replace(",", ""))
                    except: pass
                if "ts/tv" in label and len(cells) > 1:
                    try: titv = float(cells[1])
                    except: pass

    # 2) Fallback: global text regex (just in case)
    if titv is None or Ti is None or Tv is None:
        flat = re.sub(r"\s+", " ", txt)
        if Ti is None:
            m = re.search(r"Transitions[^0-9]*([0-9,]+)", flat, re.I)
            if m: Ti = int(m.group(1).replace(",", ""))
        if Tv is None:
            m = re.search(r"Transversions[^0-9]*([0-9,]+)", flat, re.I)
            if m: Tv = int(m.group(1).replace(",", ""))
        if titv is None:
            m = re.search(r"Ti/Tv[^0-9]*([0-9]+(?:\.[0-9]+)?)", flat, re.I)
            if m: titv = float(m.group(1))

    # 3) Compute if needed
    if titv is None and Ti is not None and Tv not in (None, 0):
        titv = Ti / Tv

    # 4) Haldane–Anscombe corrected ratio if counts present
    titv_corr = None
    if Ti is not None and Tv is not None:
        titv_corr = (Ti + 0.5) / (Tv + 0.5) if titv is None else titv

    return Ti, Tv, titv, titv_corr

# ---- crawl files
files = []
for pat in PATTERNS:
    files.extend(glob.glob(pat))
files = sorted(files)
if not files:
    raise SystemExit("No *_snpEff_summary.html found under expected directories.")

rows = []
for f in files:
    p = Path(f)
    species = infer_species_from_path(p)
    host_id = host_id_from_basename(p.name)
    Ti, Tv, TiTv, TiTv_corr = parse_titv_with_bs4(f)
    rows.append({
        "Sample": p.name.replace("_snpEff_summary.html",""),
        "Species": species,
        "HostID": host_id,
        "Ti": Ti, "Tv": Tv, "TiTv": TiTv,
        "TiTv_corrected": TiTv_corr,
        "File": str(p)
    })

df = pd.DataFrame(rows).sort_values(["Species","HostID","Sample"]).reset_index(drop=True)
summ = (
    df.groupby("Species")["TiTv_corrected"]
      .agg(median="median",
           q1=lambda x: x.quantile(0.25),
           q3=lambda x: x.quantile(0.75),
           n="count")
      .reset_index()
)

df.to_csv("titv_per_sample.csv", index=False)
summ.to_csv("titv_species_summary.csv", index=False)

print("Wrote titv_per_sample.csv and titv_species_summary.csv")
print(summ.to_string(index=False))
