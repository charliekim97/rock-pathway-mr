## Task 1 — DIAMANTE 2022 EUR T2D (Mahajan, Nat Genet 2022;54:560; 80,154 cases / 853,816 ctrl)
## Removes the "ROCK1/RHOA absent in Xue T2D" limitation. Xue result kept as Supplementary sensitivity.
## DOWNLOAD: T2D Knowledge Portal kp4cd.org (node/872) or diagram-consortium access form -> place at CONFIG$diamante_t2d
## BUILD: confirm in readme (CONFIG$diamante_build); rsID matching used regardless.
source("R/config.R"); source("R/helpers.R"); source("R/driver.R")

## --- read & map columns. EDIT col names to match the DIAMANTE file header. ---
COL <- list(SNP="rsID", chr="chromosome", EA="effect_allele", OA="other_allele",
            beta="Beta", se="SE", p="Pval")              # <-- CONFIRM against header
d <- fread(cmd=paste("gzip -dcf", shQuote(path.expand(CONFIG$diamante_t2d))))
setnames(d, unlist(COL[names(COL)]), names(COL), skip_absent=TRUE)
d <- d[get("chr") %in% c(2,3,18)]                         # ROCK2/RHOA/ROCK1 chromosomes
out <- d[SNP %like% "^rs", .(SNP, EA=toupper(EA), OA=toupper(OA), beta=as.numeric(beta), se=as.numeric(se), p=as.numeric(p))]

s_t2d <- 80154/(80154+853816)
mr  <- analyze_mr(out, "T2D_DIAMANTE2022")
col <- analyze_coloc(out, "T2D_DIAMANTE2022", type2="cc", s=s_t2d, N2=933970)

fwrite(mr,  file.path(CONFIG$out,"task1_diamante_t2d_MR.csv"))
fwrite(col, file.path(CONFIG$out,"task1_diamante_t2d_coloc.csv"))
print(mr); print(col)
cat("\nIf ROCK1/RHOA leads were absent, proxy column shows the r2>0.8 substitute used.\n")
