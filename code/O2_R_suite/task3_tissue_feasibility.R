## Task 3 — liver/adipose ROCK2 eQTL FEASIBILITY (v3: gated, Supplementary only).
## Step A: is there a usable instrument (F>10) in liver/adipose? If not -> say so, stop.
## Step B (only if F>10): MR/coloc on 4 outcomes ONLY (NAFLD, cirrhosis, hepatic fat, BMI). No 7-trait expansion.
## Keeps the frozen-17 declaration intact; tissue work is sensitivity, not a new analysis axis.
source("R/config.R"); source("R/helpers.R"); source("R/eqtl_catalogue.R")

tissues <- c("liver","adipose subcutaneous","adipose visceral")
feas <- rbindlist(lapply(tissues, function(ti){
  ds <- ec_dataset_id("GTEx", ti); if (is.na(ds)) return(data.table(tissue=ti, dataset=NA, lead=NA, F=NA, usable=FALSE))
  e <- ec_fetch_gene(ds, "ENSG00000134318")
  if (!nrow(e)) return(data.table(tissue=ti, dataset=ds, lead=NA, F=NA, usable=FALSE))
  e[, Z := beta/se]; lead <- e[which.max(abs(Z))]
  data.table(tissue=ti, dataset=ds, lead=lead$SNP, F=(lead$Z)^2, usable=(lead$Z)^2 > 10)
}), fill=TRUE)
fwrite(feas, file.path(CONFIG$out,"task3_tissue_feasibility.csv")); print(feas)

## ---- GATE: only proceed for tissues with F>10 ----
usable <- feas[usable==TRUE]
if (!nrow(usable)){ cat("\nNo liver/adipose ROCK2 instrument with F>10. Report 'instrument unavailable/weak'. STOP.\n"); quit(save="no") }

## Step B — restricted to 4 outcomes, Supplementary. (wire outcome readers as available)
## Example for one usable tissue with NAFLD/cirrhosis/BMI:
for (k in seq_len(nrow(usable))){
  ds <- usable$dataset[k]; ti <- usable$tissue[k]
  e <- ec_fetch_gene(ds, "ENSG00000134318")[, .(SNP,beta,se)]
  for (oc in list(list(lab="NAFLD", dt=read_harmonised(CONFIG$ghodsian_nafld), type="cc", s=8434/(8434+770180)),
                  list(lab="cirrhosis", dt=read_finngen(CONFIG$finngen_cirr), type="cc", s=5545/(5545+494803)),
                  list(lab="BMI", dt=read_pulit_bmi(CONFIG$pulit_bmi), type="quant", s=NULL))){
    r <- run_coloc(e, oc$dt[SNP %in% e$SNP,.(SNP,beta,se)], type2=oc$type, s=oc$s)
    cat(sprintf("[Supp] ROCK2 %s x %s: PP4=%.3f (n=%d)\n", ti, oc$lab, as.numeric(r$PP4), r$nsnp))
  }
}
cat("\nAll tissue results are SUPPLEMENTARY; main conclusion stays blood+muscle.\n")
