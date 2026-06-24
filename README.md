# rock-pathway-mr

**cis-eQTL Mendelian randomization + colocalization of the RhoA–ROCK pathway (ROCK2 / ROCK1 / RHOA) against cardiometabolic and fibrosis/liver traits.**

Public-data, lab-independent, fully reproducible. This repository accompanies the manuscript:

> *Human genetics of RhoA–ROCK expression highlights limits of expression-based validation for an intracellular kinase target.* Kim S. & Kim Y.-B. (manuscript; DOI on publication).

---

## Overview

The RhoA–ROCK pathway has been a candidate cardiometabolic/fibrosis drug target for two
decades, and a selective ROCK2 inhibitor (belumosudil) is clinically available. We ask
whether **genetically proxied expression** of ROCK2, ROCK1 and RHOA is causally associated
with human cardiometabolic and liver outcomes, using two-sample **cis-eQTL MR** and
**Bayesian colocalization** (single-variant `coloc.abf` and multi-signal `coloc.susie`),
benchmarked against tissue-appropriate positive controls.

**Headline result —** no robust evidence that constitutive ROCK-pathway expression affects
these traits; the null is interpreted as a limit of *expression-based* validation for an
intracellular, activity-dependent kinase (expression ≠ catalytic activity), not as evidence
against pharmacological ROCK2 inhibition.

| Analysis | Result |
|---|---|
| MR, 3 genes × 10 primary outcomes × 5 tissues | null; no PP4 exceeded 0.40 |
| Well-instrumented ROCK2 / ROCK1 colocalization | all PP4 ≤ 0.07 |
| T2D (European DIAMANTE / Finnish FinnGen / trans-ancestry) | null in all three |
| Multi-signal colocalization (`coloc.susie`, incl. nominal loci) | all SuSiE-evaluable PP.H4 < 0.11 |
| Positive controls | KCNJ11→T2D PP4 = 0.89; liver SORT1→LDL PP4 = 1.00 |
| Rare LoF burden (Genebass) | non-significant (underpowered) |

---

## Repository layout

```
code/
  python_pipeline/   Python: instrument freeze, eQTLGen Z→β/SE conversion, LD clumping,
                     MR (Wald / IVW / LD-corrected), coloc.abf, regional fetch,
                     positive controls, fibrosis axis, and figure generation
                     (fig.py → Fig 1–4,6 ; fig5.py → conceptual schematic)
  O2_R_suite/        R pipeline run on an HPC cluster (config / helpers / driver +
                     task0–task11); task4_susie_coloc.R = conditional colocalization
  susie_runner/      run_susie.R — coloc.susie runner for the 8 loci (+ positive control)
data/
  processed/         small reproducibility-grade inputs committed here
    instruments_frozen.csv             frozen instrument set (gene, SNP, tier, alleles, β, se, eaf)
    eqtl_region_converted.csv          ROCK2/KCNJ11/SORT1/GIPR cis-region eQTL
    rho_eqtl_region_conv.csv,
      rho_region_full_eqtl.txt         ROCK1/RHOA cis-region eQTL
    region_full_eqtl.txt               full ROCK2 cis-region eQTL (with alleles)
    finngen_*_region.tsv,
      nafld_meta_region.tsv            outcome cis-region extracts
    susie_inputs/                      per-locus z/β + 1000G-EUR LD + meta for coloc.susie
  README.md          how to obtain the raw (large) public summary statistics
results/             per-analysis reports (.md) + result tables (.csv/.txt), block0–block7
figures/             Fig 1–6 (PNG + PDF)
REPO_METADATA.md     copy-paste GitHub description / topics / Zenodo metadata
LICENSE              MIT
```

---

## Methods in brief

- **Instruments:** cis-eQTLs (±1 Mb) from eQTLGen whole blood and GTEx (muscle, liver,
  subcutaneous & visceral adipose); eQTLGen Z→β/SE per Zhu et al. 2016; LD clumping in
  1000 Genomes EUR; lead variant primary, near-independent / muscle sets as sensitivity.
- **Outcomes (10 primary):** MRI hepatic fat, fasting insulin, fasting glucose, HbA1c,
  2-h glucose, T2D, BMI (cardiometabolic) + cirrhosis, NAFLD, ICD-10 K74 (fibrosis/liver);
  imaging-NAFLD replication; small NASH GWAS in supplement.
- **MR:** Wald ratio / IVW / LD-corrected IVW; Steiger filtering; cross-tissue concordance.
- **Colocalization:** `coloc.abf` (PP4 threshold 0.75) over full cis regions; conditional
  multi-signal `coloc.susie` (susieR) at nominal + selected loci, LD from 1000G EUR.
- **Controls / robustness:** KCNJ11→T2D, tissue-matched SORT1→LDL; prior- and
  power-sensitivity; multi-ancestry T2D (DIAMANTE EUR + trans-ancestry, FinnGen); Genebass.

---

## Reproduce

**Requirements:** Python 3 (numpy, pandas, pysam, matplotlib); R ≥ 4.2 (coloc ≥ 5.2.3,
susieR, data.table); plink + 1000 Genomes phase-3 EUR for LD.

```bash
# 1. obtain raw public summary stats (not stored here) — see data/README.md
# 2. MR + coloc.abf + figures
python code/python_pipeline/freeze.py        # freeze instruments
python code/python_pipeline/mr_rock2.py      # MR (ROCK2); mr_rho.py for ROCK1/RHOA
python code/python_pipeline/coloc_all.py     # coloc.abf across genes × outcomes
python code/python_pipeline/fig.py           # Figures 1–4, 6
python code/python_pipeline/fig5.py          # Figure 5 (conceptual)

# 3. conditional (multi-signal) colocalization
Rscript code/susie_runner/run_susie.R        # uses data/processed/susie_inputs/
#   full HPC pipeline: code/O2_R_suite/  (edit config.R paths, then run_all.sh)
```
Processed inputs in `data/processed/` let you re-run MR/coloc and `coloc.susie` without
re-downloading the multi-GB raw files.

---

## Data sources (public)

eQTLGen · GTEx / eQTL Catalogue · GWAS Catalog · MAGIC · GLGC · FinnGen R12 ·
DIAMANTE (Mahajan et al. 2022) · Genebass · 1000 Genomes.
Accessions, builds and download links: see `data/README.md` and the manuscript Methods.
Raw summary statistics are **not** redistributed here.

---

## Citation

If you use this code or the processed inputs, please cite the manuscript (above) and this
repository (Zenodo DOI to be added on release). See `REPO_METADATA.md`.

## License

MIT — see [LICENSE](LICENSE).

## Contact

Sangeon (Charles) Kim · Young-Bum Kim — Division of Endocrinology, Diabetes, and
Metabolism, Beth Israel Deaconess Medical Center and Harvard Medical School.
Issues and questions: please use the GitHub issue tracker.
