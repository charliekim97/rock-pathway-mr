## Task 7 — coloc prior (p12) sensitivity. Cheap. Shows low PP4 holds across p12 = 1e-6..1e-4.
source("R/config.R"); source("R/helpers.R")
sweep_p12 <- function(e, od, type2, s=NULL, N2=NULL, label){
  grid <- c(1e-6,5e-6,1e-5,5e-5,1e-4)
  rbindlist(lapply(grid, function(p12){
    r <- run_coloc(e, unique(od[,.(SNP,beta,se)],by="SNP"), type2=type2, s=s, N2=N2, p12=p12)
    data.table(pair=label, p12=p12, nsnp=r$nsnp, PP4=as.numeric(r$PP4))
  }))
}
## wired pairs (outcome files in CONFIG). Add ROCK1-hepatic fat / KCNJ11-T2D once their files are placed.
res <- rbindlist(list(
  sweep_p12(load_eqtl_region("ROCK2"), read_pulit_bmi(CONFIG$pulit_bmi),       "quant",                       label="ROCK2_BMI"),
  sweep_p12(load_eqtl_region("ROCK2"), read_harmonised(CONFIG$ghodsian_nafld), "cc", s=8434/(8434+770180),    label="ROCK2_NAFLD"),
  sweep_p12(load_eqtl_region("ROCK2"), read_finngen(CONFIG$finngen_cirr),      "cc", s=5545/(5545+494803),    label="ROCK2_cirrhosis"),
  sweep_p12(load_eqtl_region("SORT1"), read_harmonised(CONFIG$glgc_ldl),       "quant",                       label="SORT1_LDL")
))
fwrite(res, file.path(CONFIG$out,"task7_prior_sensitivity.csv")); print(dcast(res, pair~p12, value.var="PP4"))
cat("\nConclusion holds if PP4 stays below threshold across the whole p12 grid.\n")
