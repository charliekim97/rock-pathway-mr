# ROCK pathway MR — O2 reinforcement suite (v3)

Defensive analyses to harden the 11-trait expression-level null against reviewer "underpowered" attacks.
**The conclusion does not change.** Every task either removes a stated limitation or pre-empts a reviewer question.

## Quick start
```bash
# 1) put this folder on O2, edit paths in R/config.R (CONFIG list)
# 2) install R deps:  install.packages(c("data.table","coloc","susieR","jsonlite","ggplot2"))
#    plink 1.9 on PATH; 1000G EUR plink bfiles at CONFIG$ld_plink_bfile
# 3) place/download outcome files (table below) into CONFIG$data
bash run_all.sh            # or run R/taskN_*.R individually
```

## Files
| File | Purpose |
|---|---|
| `R/config.R` | **edit paths here** (data, out, LD bfile, outcome files, priors) |
| `R/helpers.R` | Zhu Z→β/se, harmonise, MR (Wald/IVW), `run_coloc` (coloc.abf), outcome readers |
| `R/driver.R` | `analyze_mr`, `analyze_coloc`, plink LD-proxy |
| `R/eqtl_catalogue.R` | GTEx-uniform cis-eQTL via eQTL Catalogue REST (avoids 462GB) |
| `data/instruments_frozen.csv` | **frozen 17-SNP instruments** (EA/OA/β/se/eaf) — do not recompute |
| `data/eqtl_region_converted.csv` | staged eQTL cis regions: ROCK2, SORT1, KCNJ11, GIPR |
| `data/rho_eqtl_region_conv.csv` | staged eQTL cis regions: ROCK1, RHOA |
| `data/rock2_eqtlgen_converted.txt` | for task0 round-trip check |

## Tasks (run order)
| # | Script | Removes / shields | Out |
|---|---|---|---|
| 0 | `task0_sanity.R` | integrity: Z-conversion round-trip + cross-study **direction** (not magnitude) | direction PDF |
| 1 | `task1_diamante_t2d.R` | **T2D coverage** limit (Xue→DIAMANTE 2022 EUR) | MR+coloc csv |
| 1b | `task1b_chen_nafld.R` | NAFLD power (Ghodsian **+ Chen 2023** dual) | dual csv |
| 2 | `task2_sort1_ldl_dense.R` | positive control (SORT1 dense LDL, PP4) | coloc csv |
| 7 | `task7_prior_sensitivity.R` | "prior arbitrary" | p12×PP4 table |
| 9 | `task9_power_mde.R` | "underpowered" → MDE (tone-down legend) | power fig |
| 4 | `task4_susie_coloc.R` | single-causal-variant assumption (SuSiE, 6 loci, RHOA excluded) | susie csv |
| 5 | `task5_regional_table.R` | signal-present vs absent (regional min-p + locuscompare) | table + png |
| 6 | `task6_muscle_coloc.R` | muscle nominal (BMI, cirrhosis) | coloc csv |
| 3 | `task3_tissue_feasibility.R` | tissue scope (liver/adipose, **gated**, Supplementary) | feasibility csv |
| 10/11 | `task10_11_lookups.md` | Genebass pLoF / ROCK2 sQTL (web lookups) | hand tables |

## Outcome downloads (place in CONFIG$data; rsID matching throughout)
| Outcome | Accession / source | Build | URL |
|---|---|---|---|
| DIAMANTE T2D EUR (80,154/853,816) | Mahajan 2022 | confirm | kp4cd.org (node/872) **or** diagram-consortium access form |
| NAFLD Chen 2023 (66,814 imaging) | **GCST90271622** (confirmed; multi-ancestry) | hg38 | GWAS Catalog FTP `GCST90271001-GCST90272000/GCST90271622/harmonised/` |
| NAFLD Ghodsian 2021 (8,434/770,180) | GCST90091033 | hg38 | GWAS Catalog FTP (have it) |
| LDL GLGC (~1.3M) | GCST90239658 | confirm | GWAS Catalog FTP `GCST90239001-GCST90240000/GCST90239658/harmonised/` |
| Cirrhosis (5,545) | FinnGen R12 CIRRHOSIS_BROAD | hg38 | `storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_CIRRHOSIS_BROAD.gz` |
| BMI (806,834) | Pulit 2019 GCST009004 | hg19 | GWAS Catalog FTP (have it) |

## ⚠ Confirm before running (flagged in scripts)
- **DIAMANTE column names** (`task1` `COL` list) — set to the actual file header.
- **DIAMANTE build** (`CONFIG$diamante_build`) from its readme.
- **Chen NAFLD ancestry**: prefer EUR-only file if released; else trans-ancestry → state in MS + LD caution.
- **eQTL Catalogue dataset_id** for GTEx muscle/liver/adipose (`ec_dataset_id` resolves; verify printed id).

## v3 scope (deliberately trimmed — trimmed = more honest/defensible)
- **Task 8 (Egger/MR-PRESSO) DELETED** — cis instruments are 1–3 SNPs → meaningless. Methods one-liner instead (see B9 below).
- **Task 5 permutation DELETED** → regional min-p + locuscompare (same job, no fragile "matched SNP set").
- **Task 3 gated** — feasibility first; analysis only if F>10, 4 outcomes, Supplementary.
- **Task 4** excludes RHOA (single weak instrument).
- **Task 9 / Figure 5 toned down** — "limits of expression MR for intracellular kinases", not "MR fails".

## Manuscript edits (not scripts — do in build_ms.js / docx)
- B9: Methods one-liner — "Because most cis instruments comprised one to three variants and residual LD was present in the broader sensitivity sets, MR-Egger and MR-PRESSO were not treated as informative primary pleiotropy diagnostics; pleiotropy was addressed through colocalization, cross-tissue concordance, conditional (SuSiE) colocalization and regional visualization."
- B7: add **Sveinbjornsson 2022 deCODE** (Nat Genet 54:1652) as Discussion contrast (MTARC1/GPAM protective LoF validated vs ROCK2 null) — citation only, not an outcome.
- B8: NAFLD dual (Chen + Ghodsian) in Methods/Results + outcome table row.
- B6: NASH → Supplementary (main = 10 traits).
- ref16 = Commun Biol 2023;6:1176 ; ref18 = AJRCMB 2024;71(4):430 — fill first authors.

## ★ HARD STOP
Tasks 1–11 + Figure 5 = dry completion line. **No single-cell, no drug-signature projection, no wet.**
Wet (HSC/CRISPRi/phosphoproteomics) = next project (PhD Ch.2). Even if task 10/11 shows signal → decide separately, no pre-emptive expansion.
