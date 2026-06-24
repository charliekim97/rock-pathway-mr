## Task 2 — SORT1 -> LDL on dense GLGC (Graham 2021, GCST90239658, EUR ~1.3M)
## Completes positive controls: KCNJ11 (disease) + SORT1 (lipid). Target PP4 0.013 -> >0.75.
## DOWNLOAD (GWAS Catalog harmonised):
##   https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90239001-GCST90240000/GCST90239658/harmonised/
source("R/config.R"); source("R/helpers.R")

## SORT1 eQTL region is staged in eqtl_region_converted.csv (GeneSymbol==SORT1)
e_sort1 <- load_eqtl_region("SORT1")
stopifnot(nrow(e_sort1) > 50)

out_ldl <- read_harmonised(CONFIG$glgc_ldl)[SNP %in% e_sort1$SNP, .(SNP, beta, se)]
cat("SORT1 eQTL SNPs:", nrow(e_sort1), " | overlap with GLGC LDL:", nrow(out_ldl), "\n")

r <- run_coloc(e_sort1, out_ldl, type2="quant")          # LDL = quantitative
res <- data.table(pair="SORT1_x_LDL_GLGC", nsnp=r$nsnp, PP0=r$PP0, PP1=r$PP1, PP2=r$PP2, PP3=r$PP3, PP4=r$PP4)
fwrite(res, file.path(CONFIG$out,"task2_sort1_ldl_coloc.csv")); print(res)
cat(ifelse(as.numeric(r$PP4)>=CONFIG$coloc_pp4_threshold,
           ">> PASS: SORT1 colocalizes on dense LDL (M3 closed).\n",
           ">> still low — check overlap N and build; consider SuSiE (task4).\n"))
## locuszoom: see task5 (reuse SORT1 region) or use locuscomparer with eQTL vs LDL p-values.
