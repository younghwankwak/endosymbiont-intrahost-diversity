#!/usr/bin/env bash
set -euo pipefail

module load bcftools

NAS_ROOT=""
SUL_ROOT=""
MITO_ROOT=""   # <- adjust if needed

NAS_SUM="${NAS_ROOT}/nas_af_summary_ref.tsv"
SUL_SUM="${SUL_ROOT}/sul_af_summary_ref.tsv"
MITO_SUM="${MITO_ROOT}/mito_af_summary_ref.tsv"
COMBINED="af_combined_ref.tsv"

# headers
echo -e "Sample\tCHROM\tPOS\tREF\tALT\tDP\tAF" > "$NAS_SUM"
echo -e "Sample\tCHROM\tPOS\tREF\tALT\tDP\tAF" > "$SUL_SUM"
echo -e "Sample\tCHROM\tPOS\tREF\tALT\tDP\tAF" > "$MITO_SUM"
echo -e "Species\tSample\tCHROM\tPOS\tREF\tALT\tDP\tAF" > "$COMBINED"

extract_dir () {
  local species_name="$1"   # "Nasuia" | "Sulcia" | "Mitochondria"
  local root="$2"           # directory containing *_snpEff.vcf.gz
  local summary="$3"        # output summary path for this species

  shopt -s nullglob
  for vcf in "${root}"/*.replaced.vcf.gz; do
    [[ -f "$vcf" ]] || continue

    # index if needed (best effort)
    if [[ ! -s "${vcf}.tbi" && ! -s "${vcf}.csi" ]]; then
      bcftools index -f "$vcf" >/dev/null 2>&1 || true
    fi

    # get samples in this VCF (expect 1; handle >1 safely)
    mapfile -t samples < <(bcftools query -l "$vcf")
    if [[ "${#samples[@]}" -eq 0 ]]; then
      echo "[WARN] no samples in: $vcf"
      continue
    fi

    for sid in "${samples[@]}"; do
      # Extract per-sample DP, AF; single line per site
      bcftools query -s "$sid" -f '%CHROM\t%POS\t%REF\t%ALT\t[%DP\t%AF]\n' "$vcf" \
      | awk -v s="$sid" -v sp="$species_name" -v specsum="$summary" -v comb="$COMBINED" '
          BEGIN{OFS="\t"}
          {
            chrom=$1; pos=$2; ref=$3; alt=$4; dp=$5; af=$6;
            if (dp=="." || dp=="") dp=0;
            if (af=="." || af=="") af="0";
            # split multiallelic ALT and AF; align lengths; coerce "." -> 0
            nA=split(alt, A, /,/);
            nF=split(af,  F,  /,/);
            if (nF < nA) { for (i=nF+1; i<=nA; i++) F[i]=0 }
            # if AF has extra fields, ignore beyond nA
            for (i=1; i<=nA; i++) {
              if (F[i]=="." || F[i]=="") F[i]=0;
              print s, chrom, pos, ref, A[i], dp, F[i] >> specsum;
              print sp, s, chrom, pos, ref, A[i], dp, F[i] >> comb;
            }
          }'
    done

    echo "✓ ${species_name}: $(basename "$vcf")"
  done
}

extract_dir "Nasuia"       "$NAS_ROOT"  "$NAS_SUM"
extract_dir "Sulcia"       "$SUL_ROOT"  "$SUL_SUM"
extract_dir "Mitochondria" "$MITO_ROOT" "$MITO_SUM"

echo "Wrote:"
echo "  - $NAS_SUM"
echo "  - $SUL_SUM"
echo "  - $MITO_SUM"
echo "  - $COMBINED"

