# Intrahost mutational dynamics parallel long-term genome evolution in endosymbionts

**Younghwan Kwak** and **Gordon M. Bennett**

---

## Overview

This repository contains custom analysis scripts and processed variant call files (VCFs) used in downstream analyses. Large datasets, intermediate files, and scripts for standard bioinformatic tools are available upon request.

Kwak Y, Bennett GM. 2026. *Intrahost mutational dynamics parallel long-term genome evolution in endosymbionts.*  
bioRxiv 2026.01.27.701382  
https://www.biorxiv.org/content/10.64898/2026.01.27.701382v1

The scripts are organized by analysis modules and implement downstream processing of variant data, including classification, allele frequency summaries, nucleotide diversity, and pN/pS analyses.

All raw sequencing data, assemblies, and associated metadata are available under:

NCBI BioProject: PRJNA1391824

---

## Repository structure

### `01_class/`
Scripts for variant classification and annotation-based summaries.

- `indel_from_snpEff.py`  
  Extracts and classifies indels from SnpEff-annotated VCF outputs.

- `indel_impact.py`  
  Summarizes predicted functional impacts of indels from annotation results.

- `titv_from_snpeff.py`  
  Calculates transition/transversion (Ti/Tv) patterns from annotated SNP data.

### `02_afs/`
Scripts for allele frequency spectrum analyses.

- `af_bin.py`  
  Bins allele frequency values for downstream summarization and visualization.

- `bcftools_af.sh`  
  Extracts allele frequency information from VCF files using `bcftools`.

### `03_pi/`
Scripts for nucleotide diversity analysis.

- `pi.py`  
  Computes nucleotide diversity (π) and related summary statistics from variant data.

### `04_pnps/`
Scripts for pN/pS analyses at multiple levels.

- `1.pnps_norm.py`  
  Performs normalization steps required for pN/pS calculations.

- `2.pnps_gene.py`  
  Calculates gene-level pN/pS ratios.

- `3.pnps_cog.py`  
  Summarizes pN/pS patterns across COG functional categories.
  
## Contact

For questions, clarification, or access to additional data and analysis workflows, please contact:

Younghwan Kwak  
University of California, Merced  
Email: ykwak@ucmerced.edu