## Task 1b — Chen/Du 2023 NAFLD (Nat Genet 2023;55:1640; GCST90271622; imaging 66,814 + dx)
## Doubles NAFLD axis (Ghodsian case-control + Chen imaging). Both null -> M2 closed at NAFLD.
## DOWNLOAD (GWAS Catalog harmonised, hg38):
##   https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90271001-GCST90272000/GCST90271622/harmonised/
## NOTE: multi-ancestry. Prefer an EUR-only file if released; else trans-ancestry (state in MS, LD caution).
source("R/config.R"); source("R/helpers.R"); source("R/driver.R")

out_chen <- read_harmonised(CONFIG$chen_nafld)          # hm_rsid etc.
out_chen <- out_chen[, chr := NA]                        # rsID matching only
mr_c  <- analyze_mr(out_chen,  "NAFLD_Chen2023")
## Chen case fraction (diagnostic component ~3,584 / 621,081); imaging is quantitative-ish but distributed as cc here.
s_chen <- 3584/(3584+621081)
col_c <- analyze_coloc(out_chen, "NAFLD_Chen2023", type2="cc", s=s_chen)

## re-run Ghodsian for side-by-side (existing file)
out_gh <- read_harmonised(CONFIG$ghodsian_nafld)
mr_g  <- analyze_mr(out_gh,  "NAFLD_Ghodsian2021")
col_g <- analyze_coloc(out_gh, "NAFLD_Ghodsian2021", type2="cc", s=8434/(8434+770180))

mr  <- rbind(mr_c, mr_g); col <- rbind(col_c, col_g)
fwrite(mr,  file.path(CONFIG$out,"task1b_nafld_dual_MR.csv"))
fwrite(col, file.path(CONFIG$out,"task1b_nafld_dual_coloc.csv"))
print(col); cat("\nBoth NAFLD definitions should be null for ROCK2 (PP4<0.1) to close M2.\n")
