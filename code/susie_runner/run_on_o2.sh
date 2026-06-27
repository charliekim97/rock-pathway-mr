#!/bin/bash
## Run the coloc.susie audit on Harvard O2 — reproduce the ORIGINAL environment.
## Usage:  bash run_on_o2.sh  [INDIR]  [OUTDIR]
set -e

## 1) same toolchain as the original run (adjust module names if yours differ)
module load gcc/14.2.0 R/4.5.2 2>/dev/null || module load R/4.5.2 2>/dev/null || true

## 2) inputs = the frozen *.betas.tsv / *.ld.gz / *.meta (8 loci)
INDIR=${1:-$HOME/rock2_o2/susie_inputs}
OUTDIR=${2:-$HOME/rock2_o2/susie_audit}

## if the inputs are not on O2 anymore, upload them from the repo first, e.g.:
##   scp -r data/processed/susie_inputs  USER@o2.hms.harvard.edu:~/rock2_o2/
echo "INDIR=$INDIR  OUTDIR=$OUTDIR"
ls "$INDIR"/*.betas.tsv >/dev/null 2>&1 || { echo "ERROR: no *.betas.tsv in $INDIR"; exit 1; }

## 3) confirm the CITED versions (manuscript says coloc 5.2.3)
Rscript -e 'cat("coloc",as.character(packageVersion("coloc")),"| susieR",as.character(packageVersion("susieR")),"| R",R.version.string,"\n")'
## If coloc is NOT 5.2.3, pin it (do NOT silently use 6.x and call it robustness):
##   Rscript -e 'remotes::install_version("coloc","5.2.3",repos="https://cloud.r-project.org",upgrade="never")'

## run on a compute node (not the login node):
##   srun -t 0:30:00 --mem=8G -c 2 bash run_on_o2.sh
Rscript audit_susie.R "$INDIR" "$OUTDIR"

echo
echo "DONE. Rebuilt Table S1 -> $OUTDIR/TableS1_rebuilt.tsv"
echo "Send back: TableS1_rebuilt.tsv  +  every <locus>/*_outcome_CS.tsv  +  *_outcome_CS_sensitivity.tsv"
