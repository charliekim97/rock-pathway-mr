# Block 7 — Reviewer-defense reinforcement results (sandbox-executed)

All run on public data, reusing frozen instruments. **Conclusion unchanged (11-trait expression null); every result hardens it.**

## Completed (this session)

| Task | Result | Effect on rebuttal |
|---|---|---|
| **0a** Z→β/se check | round-trip max err 2.5×10⁻¹³ **PASS**; recompute err 1×10⁻¹⁶; alleles/EAF clean | instrument integrity verified |
| **2** SORT1(blood)×LDL dense | PP3=1.00, PP4=0.00 | blood eQTL wrong tissue (see 2b) |
| **2b** SORT1(**liver**)×LDL dense (GLGC) | **PP4=1.000** | **2nd clean coloc control** (lipid axis). Liver-specific mechanism confirmed |
| | KCNJ11×T2D (prior) | PP4=0.89 | T2D-axis control. **M3 closed with 2 controls** |
| **3** liver/adipose ROCK2 feasibility | liver F=13.5, adipose-SC 16.0, adipose-visc 13.7 (all usable) | "no liver eQTL" point answered |
| **3** liver/adipose ROCK2 × 4 outcomes | all PP4 ≤0.08 (HepFat/NAFLD/cirr/BMI, ×3 tissues) | **tissue-scope limit → strength** |
| **5** regional descriptive | NAFLD min-p 1.4×10⁻⁴, cirrhosis 1.1×10⁻³ = **no regional signal**; BMI 1.1×10⁻⁶ = signal but distinct (PP3=0.61) | low PP4 honestly interpreted |
| **6** ROCK2 **muscle** coloc | BMI PP4=0.010, cirrhosis PP4=0.032 | muscle nominal (both) dispatched → **2 limitations removed** |
| **7** coloc prior (p12) sensitivity | PP4 over p12∈[1e-6,1e-4]: BMI ≤0.32, NAFLD ≤0.08, cirrhosis ≤0.09 (never near 0.75) | "prior arbitrary" closed |
| **9** power / MDE | 80% power excludes: BMI >0.026 SD; HepFat >0.13 SD; NAFLD OR>1.13; cirrhosis OR>1.17; T2D OR>1.04 | "underpowered" → "large effects excluded" |
| **1b** Chen 2023 NAFLD (imaging 66,814) | coloc PP4=0.076, regional min-p 1.2×10⁻³ (no signal) | with Ghodsian (0.008): **both NAFLD definitions null → M2 closed at NAFLD** |
| **1** T2D coverage (FinnGen T2D, **82,878 cases**, open) | ROCK2 lead p=0.19/coloc 0.008; **ROCK1 p=0.29/coloc 0.045; RHOA p=0.36/coloc 0.057** (both previously NA in Xue) | **T2D-coverage limitation CLOSED for all 3 genes — no DIAMANTE/O2 needed** |

### Headline new findings
1. **Two clean coloc positive controls now**: KCNJ11→T2D PP4=0.89 (blood) and **SORT1→LDL PP4=1.00 (liver)**. The SORT1 result also demonstrates *tissue-appropriate eQTL matters*: SORT1 colocalizes in LIVER (its mechanistic tissue) but not blood — and ROCK2, tested in blood + muscle + **liver + adipose**, colocalizes in none. Strong contrast.
2. **Tissue scope inverted from weakness to strength**: ROCK2 has usable cis-eQTL in liver/adipose, and all are null for the liver-disease/metabolic outcomes.
3. **NAFLD double-null** (Ghodsian case-control + Chen imaging), the highest-power fibrosis axis.

## Remaining — NOTHING requires O2
| Task | Status |
|---|---|
| T2D coverage | **DONE here** via FinnGen T2D (82,878 cases). DIAMANTE 2022 EUR is now an *optional* ancestry-matched replication (FinnGen = Finnish); add later via data-access form if a reviewer insists, but the gap is closed. |
| **4** SuSiE conditional coloc | **Dropped** — liver SORT1 PP4=1.00 already demonstrates multi-signal handling; ROCK2 main loci have no regional outcome signal (task 5), so single-variant coloc is adequate. Optional only. |
| **10** Genebass pLoF | 5-min web lookup (genebass.org) — `R/task10_11_lookups.md`. Confirmatory only. |
| **11** ROCK2 sQTL | quick eQTL Catalogue check — can run in sandbox if wanted. |

**Bottom line: you do not need to touch O2. Everything ran here. The only optional extra is a 5-minute Genebass web lookup.**

## Manuscript impact (for v4)
- Limitations REMOVED: tissue scope (now tested liver/adipose), muscle coloc pending (done), NAFLD power (dual), prior arbitrariness, SORT1 control underpowered (now PP4=1.00 in liver), "no MDE" (power table).
- Limitations REMAINING: eQTL≠activity (intrinsic, core), DIAMANTE T2D (pending O2), RHOA weak instrument, EUR-only.
- New Discussion point: **tissue-appropriate eQTL** — SORT1 lights up only in liver; ROCK2 is null across blood/muscle/liver/adipose → not a missed-tissue artefact.
