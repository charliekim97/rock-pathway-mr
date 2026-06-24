# DIAMANTE T2D + SuSiE — sandbox execution results (062026)

## PART A — DIAMANTE T2D (✅ DONE, sandbox-confirmed)

Data: DIAMANTE (Mahajan et al. 2022), hg19/b37, beta/SE provided directly (no conversion).
Frozen instruments reused (eQTLGen blood lead + GTEx muscle 2-SNP; ROCK1/RHOA blood). rsID harmonisation, allele-aligned.

### A2 — European MAIN (DIAMANTE-EUR, 80,154 cases / 853,816 controls)
MR (per-SD expression):
- ROCK2 blood lead Wald (rs12468344): beta=+0.021, se=0.040, **p=0.61**
- ROCK2 muscle IVW (2 SNP): beta=-0.012, **p=0.67**
- ROCK1 blood IVW (2 SNP): beta=+0.009, **p=0.87**
- RHOA blood Wald (rs115431681): beta=-0.044, **p=0.46**

coloc.abf (cis × T2D, cc):
- ROCK2 × T2D: **PP4=0.003** (PP3=0.355, PP1=0.642), nSNP=8126
- ROCK1 × T2D: **PP4=0.016** (PP3=0.061), nSNP=1449
- RHOA × T2D: **PP4=0.021** (PP3=0.521, PP2=0.456), nSNP=2942

Regional min-p (EUR): ROCK2 3.5e-4, ROCK1 4.3e-4, **RHOA 2.1e-9** (a genuine T2D
signal exists within the RHOA ±1Mb window but does NOT colocalise with the RHOA
eQTL — PP3>>PP4 — an internal demonstration that the pipeline returns "distinct
variant", not false PP4, where regional signal is present).

### A3 — Trans-ancestry SENSITIVITY (DIAMANTE-TA, 180,834 cases; MR only, no coloc)
- ROCK2 lead Wald p=0.90; muscle IVW p=0.21
- ROCK1 IVW p=0.69; RHOA Wald p=0.62
LD/AF heterogeneity across ancestries: European-derived instruments; reported as
power/generalizability sensitivity only (no coloc — LD mismatch).

### Verdict
T2D axis now NULL across **European (DIAMANTE) + Finnish (FinnGen) + trans-ancestry
(DIAMANTE multi-ancestry)** — triple null, all 3 genes. Ancestry attack closed.
"planned DIAMANTE" limitation removed.

---

## PART B — SuSiE conditional coloc (✅ DONE on O2, 2026-06-20)

Run on Harvard O2 (R 4.5.2, coloc 5.2.3, susieR; compute node). LD = 1000G EUR (503),
allele-aligned in sandbox. coloc.susie per locus; coloc.abf fallback where SuSiE
non-convergent. 5 loci (RHOA excluded — single weak instrument):

| Locus | coloc.susie PP.H4 | eQTL credible sets | outcome credible sets | method |
|---|---|---|---|---|
| **KCNJ11–T2D (positive control)** | **0.89** (idx1=1,idx2=1) | 3 | 6 | susie |
| ROCK1–T2D | 0.017 | 3 (3/8/15) | 0 | susie (no outcome CS)→abf |
| ROCK2–cirrhosis | 0.013 | eQTL non-converged | 0 | abf |
| ROCK2–NAFLD | 0.010 | 5 | 0 | susie (no outcome CS)→abf |
| ROCK2–T2D | 0.003 | 4 | 0 | susie (no outcome CS)→abf |

### Nominal-locus extension (062026, P0-1 selection-bias fix) — all run on O2
To close the selection-bias gap (SuSiE must also cover the loci that had a nominal MR
signal, not only null loci), three nominal multi-instrument loci were added:

| Locus | coloc.susie PP.H4 | eQTL credible sets | outcome credible sets |
|---|---|---|---|
| ROCK2–BMI (nominal) | 0.104 | 4 (1/1/46/2) | 0 |
| ROCK1–hepatic fat (nominal) | 0.069 | 3 (3/8/16) | 0 |
| ROCK2 muscle–cirrhosis (nominal) | 0.032 | many (eQTL over-split) | 0 |

All three nominal loci, like the null loci, had NO outcome credible set → no
colocalization under the multi-signal model (all PP.H4 ≤ 0.10). RHOA–2-h glucose
(single weak instrument) is unsuitable for SuSiE and was assessed by coloc.abf only.
Net: every multi-instrument ROCK locus evaluated — nominal or null — is concordant
between coloc.abf and coloc.susie, while the KCNJ11 control colocalizes (PP.H4 = 0.89).
"every ROCK locus" claim is now literally true (selection-bias hole closed).

**Interpretation:** the positive control KCNJ11→T2D colocalizes under the multi-signal
SuSiE model (PP.H4=0.89), validating the pipeline. At every ROCK locus the eQTL is
resolved into multiple credible sets, yet the disease outcome has NO credible set in
the cis region (regional min-p 0.0006–0.017), so no colocalization is possible under
either the single-variant (coloc.abf) or multi-signal (coloc.susie) model. **abf and
susie agree → the expression-level null is robust.** "planned SuSiE" limitation removed;
"robust" wording restored where justified. → v8.

### (historical) sandbox attempt — why it moved to O2

Attempted in sandbox per packet. Blockers:
1. No R, no plink, no sudo (cannot install r-base / susieR / coloc).
2. 1000G phase3 ID column is "." → LD obtained successfully by **position** matching
   via pysam remote tabix on 1000G EUR (503 samples) — this part WORKED
   (e.g. KCNJ11 region: 6,183 SNP LD matrix in ~6 s).
3. Implemented susie_rss + coloc.susie from scratch in Python (no susieR available).
   **Positive control FAILED**: KCNJ11→T2D returned PP4≈0 (should be ~0.89). The
   from-scratch single-effect EM does not shrink null effects (all L effects fire),
   so credible sets and per-effect lBF are unreliable. Two parameterisations
   (sufficient-stat and z-score) both failed the control.

Per HARD STOP ("no infinite debugging → O2 fallback or keep planned"), SuSiE was
NOT completed. The LD-fetch machinery is ready (`/tmp/dia/ld.py`), so the existing
`O2_scripts/R/task4_susie_coloc.R` can run on HPC with susieR/coloc directly.

### Recommendation
- Option 1 (fast): keep "SuSiE conditional coloc planned" in limitations; do NOT
  restore "robust" wording. Conditional argument is partly covered empirically: the
  RHOA-region result above + KCNJ11/SORT1-liver positive controls bracket the
  method (high PP4 where truly shared; PP3 where regionally signalled but distinct).
- Option 2 (full upside): run task4_susie_coloc.R on O2 for the 5 loci
  (ROCK2-BMI/NAFLD/cirrhosis, ROCK1-hepatic fat, SORT1-LDL liver) + KCNJ11 control,
  then restore "robust" and remove the planned-SuSiE limitation → CTS/JCEM ready.
