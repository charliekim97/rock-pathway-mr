# coloc.susie audit — run on O2, then rebuild Table S1 from fact

## Why
The current Table S1 reports a `PP.H4` for ROCK loci that have **Outcome CS = 0**.
`coloc.susie()` returns one row per (exposure-signal × outcome-signal) pair, so if the
outcome has **no** credible set, there is no genuine multi-signal PP.H4 — the number is a
`coloc.abf` fallback. Before we change a single sentence in the manuscript, we must look at
the **raw SuSiE objects** and decide which of three cases is true at each locus:

- **A — outcome credible set actually exists.** Then "Outcome CS = 0" was a transcription
  error. Keep `coloc.susie`, fix the CS counts, add the signal-pair PP.H3/PP.H4 to the supplement.
- **B — outcome credible set is genuinely 0 at every reasonable threshold.** Then multi-signal
  colocalization is **not evaluable**; `coloc.abf` + the absence of any regional outcome signal
  is the evidence. (This is what v21 currently assumes — provisionally.)
- **C — components exist but are filtered out by coverage/purity.** Then do **not** call it a
  negative SuSiE test; lower the wording to "no high-confidence credible set; interpreted cautiously."

## Run (Harvard O2)
```bash
# inputs already on O2 from the original run (~/rock2_o2/susie_inputs), else upload:
#   scp -r data/processed/susie_inputs  USER@o2.hms.harvard.edu:~/rock2_o2/
module load gcc/14.2.0 R/4.5.2          # the toolchain you used originally
Rscript -e 'cat("coloc",as.character(packageVersion("coloc")),"\n")'   # must be 5.2.3
srun -t 0:30:00 --mem=8G -c 2 bash run_on_o2.sh
```
Reproduces the original call exactly (`runsusie(D, maxit=100, repeat_until_convergence=FALSE)`,
default L=10, coverage 0.95) and **adds** a coverage=0.50 and a relaxed-purity (min_abs_corr=0.1)
sensitivity so B vs C is decided empirically.

## What it writes (per locus, under `susie_audit/<locus>/`)
`*_susie_objects.rds` (raw S1/S2) · `*_exposure_CS.tsv` / `*_outcome_CS.tsv` (size, coverage,
purity) · `*_outcome_CS_sensitivity.tsv` (CS count at coverage 0.95/0.50 and relaxed purity) ·
`*_exposure_summary.txt` / `*_outcome_summary.txt` · `*_coloc_susie_summary.tsv` + `*_results.tsv`
(only if evaluable) · `*_QC.rds` · `*_row.tsv`.
Combined: **`susie_audit/TableS1_rebuilt.tsv`** with one row per locus:

| column | meaning |
|---|---|
| exposure_CS / outcome_CS | credible sets at the reported threshold (coverage 0.95, purity 0.5) |
| outcome_CS_cov50 / outcome_CS_relaxpurity | did weaker coverage / purity reveal an outcome CS? (B vs C) |
| outcome_effects_fired | susie effects with max PIP > 0.10 before purity filtering |
| regional_min_p_outcome | smallest outcome p in the window (is there anything to fine-map?) |
| n_signal_pairs / coloc_susie_maxPP_H4 | genuine coloc.susie output (only if outcome_CS ≥ 1) |
| coloc_abf_PP_H4 | the single-variant fallback (what the table currently shows) |
| method / verdict | `susie` vs `abf`; and A / B / C |

## Send back to finalise
`TableS1_rebuilt.tsv` + every `<locus>/*_outcome_CS.tsv` + every `<locus>/*_outcome_CS_sensitivity.tsv`.
With those I rebuild Table S1 and set the p.13 / p.17 / Table S1 wording to match the verdict —
**no manuscript wording changes until then.** (v21 currently holds the conservative B-wording as a
placeholder; if the result is A, we put the SuSiE confirmation back.)

## Note on what survives regardless
The ROCK2 main null rests on `coloc.abf`, absence of regional outcome signal, cross-tissue
consistency, and the primary MR — not on SuSiE. So even a "B / not evaluable" verdict does not
break the paper; it only removes the "SuSiE independently re-confirmed it" decoration.
