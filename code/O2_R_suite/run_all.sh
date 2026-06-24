#!/usr/bin/env bash
# run_all.sh — ROCK pathway MR O2 reinforcement (v3). Run from O2_scripts/.
# Order is control-anchored: sanity (0) -> max-gain (1,1b,2) -> cheap shields (7,9) -> M2 core (4,5,6) -> tissue (3) -> lookups.
# For SLURM: prepend `sbatch --wrap="Rscript ..."` or wrap each in a job script. Interactive shown here.
set -e
cd "$(dirname "$0")"
module load R 2>/dev/null || true          # O2: load R + gcc as needed
export R="Rscript --vanilla"

echo "== Task 0: sanity ==";            $R R/task0_sanity.R
echo "== Task 1: DIAMANTE T2D ==";      $R R/task1_diamante_t2d.R
echo "== Task 1b: Chen+Ghodsian NAFLD =="; $R R/task1b_chen_nafld.R
echo "== Task 2: SORT1 dense LDL ==";   $R R/task2_sort1_ldl_dense.R
echo "== Task 7: prior sensitivity =="; $R R/task7_prior_sensitivity.R
echo "== Task 9: power/MDE ==";         $R R/task9_power_mde.R
echo "== Task 4: SuSiE coloc ==";       $R R/task4_susie_coloc.R
echo "== Task 5: regional table ==";    $R R/task5_regional_table.R
echo "== Task 6: muscle coloc ==";      $R R/task6_muscle_coloc.R
echo "== Task 3: tissue feasibility =="; $R R/task3_tissue_feasibility.R
echo "Tasks 10/11 are manual web lookups — see R/task10_11_lookups.md"
echo "DONE. Outputs in $(grep -m1 'out ' R/config.R || echo CONFIG\$out)"
