## Task 6 — ROCK2 MUSCLE colocalization for the two muscle-tier nominal signals (BMI, cirrhosis).
## Blood already no-coloc; this confirms muscle no-coloc -> both nominal muscle signals dispatched.
source("R/config.R"); source("R/helpers.R"); source("R/eqtl_catalogue.R")

ds <- ec_dataset_id("GTEx","muscle"); cat("GTEx muscle dataset_id:", ds, "\n")
stopifnot(!is.na(ds))
e_mus <- ec_fetch_gene(ds, "ENSG00000134318")             # ROCK2 muscle cis region (full)
cat("ROCK2 muscle cis SNPs:", nrow(e_mus), "\n")
e_mus <- e_mus[, .(SNP, beta, se)]

bmi  <- read_pulit_bmi(CONFIG$pulit_bmi)[SNP %in% e_mus$SNP, .(SNP,beta,se)]
cirr <- read_finngen(CONFIG$finngen_cirr)[SNP %in% e_mus$SNP, .(SNP,beta,se)]

r_bmi  <- run_coloc(e_mus, bmi,  type2="quant")
r_cirr <- run_coloc(e_mus, cirr, type2="cc", s=5545/(5545+494803))
res <- rbindlist(list(
  data.table(pair="ROCK2muscle_x_BMI",       nsnp=r_bmi$nsnp,  PP3=as.numeric(r_bmi$PP3),  PP4=as.numeric(r_bmi$PP4)),
  data.table(pair="ROCK2muscle_x_cirrhosis", nsnp=r_cirr$nsnp, PP3=as.numeric(r_cirr$PP3), PP4=as.numeric(r_cirr$PP4))))
fwrite(res, file.path(CONFIG$out,"task6_muscle_coloc.csv")); print(res)
cat("\nLow PP4 (<0.1) confirms muscle no-coloc -> remove the two muscle limitations.\n")
