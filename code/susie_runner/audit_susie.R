#!/usr/bin/env Rscript
## ============================================================================
## coloc.susie AUDIT — does the OUTCOME have a SuSiE credible set, or not?
## Re-runs the 7 (+1) loci on the SAME frozen inputs / same ±250 kb window /
## same 1000G-EUR LD, and saves the RAW SuSiE objects + credible sets so that
## Table S1 can be rebuilt from fact, not from a summary line.
##
## RUN ON O2 (reproduce the original environment & pinned versions):
##   module load gcc/14.2.0 R/4.5.2          # or the modules you used originally
##   Rscript audit_susie.R  [INDIR]  [OUTDIR]
##   # INDIR  default = ~/rock2_o2/susie_inputs   (the frozen *.betas.tsv/*.ld.gz/*.meta)
##   # OUTDIR default = ~/rock2_o2/susie_audit
##
## REQUIRES coloc 5.2.3 + susieR (the versions the manuscript cites). The script
## prints the loaded versions; if coloc != 5.2.3, install it first:
##   Rscript -e 'remotes::install_version("coloc","5.2.3",repos="https://cloud.r-project.org")'
## ============================================================================

suppressMessages({library(coloc); library(susieR); library(data.table)})

args  <- commandArgs(trailingOnly = TRUE)
INDIR  <- if (length(args) >= 1) args[1] else path.expand("~/rock2_o2/susie_inputs")
OUTDIR <- if (length(args) >= 2) args[2] else path.expand("~/rock2_o2/susie_audit")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

cat("======================================================================\n")
cat("coloc.susie audit\n")
cat("  coloc  :", as.character(packageVersion("coloc")), "\n")
cat("  susieR :", as.character(packageVersion("susieR")), "\n")
cat("  R      :", R.version.string, "\n")
cat("  INDIR  :", INDIR, "\n  OUTDIR :", OUTDIR, "\n")
cat("======================================================================\n\n")

## the 7 loci the audit must cover (+ ROCK2-cirrhosis, expected exposure non-conv.)
WANT <- c("KCNJ11-T2D","ROCK2-BMI","ROCK1-hepaticfat","ROCK2muscle-cirrhosis",
          "ROCK2-T2D","ROCK2-NAFLD","ROCK1-T2D","ROCK2-cirrhosis")
have <- sub("\\.betas\\.tsv$","", list.files(INDIR, pattern="\\.betas\\.tsv$"))
loci <- intersect(WANT, have)
missing <- setdiff(WANT, have)
if (length(missing)) cat("NOTE: inputs not found for:", paste(missing, collapse=", "), "\n\n")

## ---- helpers ----------------------------------------------------------------
build_data <- function(L){
  b  <- fread(file.path(INDIR, paste0(L,".betas.tsv")))
  R  <- as.matrix(fread(file.path(INDIR, paste0(L,".ld.gz")), header = FALSE))
  kv <- readLines(file.path(INDIR, paste0(L,".meta")))
  kv <- setNames(sub(".*=","",kv), sub("=.*","",kv))
  rownames(R) <- colnames(R) <- b$SNP
  D1 <- list(beta=b$b1, varbeta=b$vb1, snp=b$SNP, LD=R, type="quant", sdY=1,
             N=as.numeric(kv["Neqtl"]))
  D2 <- list(beta=b$b2, varbeta=b$vb2, snp=b$SNP, LD=R, type=unname(kv["type2"]),
             N=as.numeric(kv["N2"]))
  if (kv["type2"]=="cc") D2$s <- as.numeric(kv["s"]) else D2$sdY <- 1
  list(b=b, R=R, kv=kv, D1=D1, D2=D2)
}

## credible sets at a chosen coverage / purity, using the LD as Xcorr
cs_at <- function(S, R, coverage, min_abs_corr){
  out <- tryCatch(susie_get_cs(S, Xcorr=R, coverage=coverage, min_abs_corr=min_abs_corr),
                  error=function(e) NULL)
  if (is.null(out) || is.null(out$cs)) return(0L)
  length(out$cs)
}

## detailed CS table (size, coverage, purity) from a fitted object's own sets
cs_table <- function(S){
  cs <- S$sets$cs
  if (is.null(cs) || !length(cs))
    return(data.table(cs_id=integer(), n_snps=integer(), coverage=numeric(),
                      purity_min_abs_corr=numeric(), purity_mean_abs_corr=numeric()))
  pur <- S$sets$purity
  data.table(cs_id = seq_along(cs),
             n_snps = vapply(cs, length, integer(1)),
             coverage = if (!is.null(S$sets$coverage)) S$sets$coverage else NA_real_,
             purity_min_abs_corr  = if (!is.null(pur)) pur[,"min.abs.corr"]  else NA_real_,
             purity_mean_abs_corr = if (!is.null(pur)) pur[,"mean.abs.corr"] else NA_real_)
}

rows <- list()

for (L in loci){
  cat("==================  ", L, "  ==================\n", sep="")
  od <- file.path(OUTDIR, L); dir.create(od, showWarnings=FALSE, recursive=TRUE)
  dat <- build_data(L); R <- dat$R; D1 <- dat$D1; D2 <- dat$D2

  ## sanity / QC
  ok_snp <- identical(D1$snp, D2$snp) && nrow(R)==ncol(R) && nrow(R)==length(D1$snp)
  saveRDS(list(check1=tryCatch(check_dataset(D1, req="LD"), error=identity),
               check2=tryCatch(check_dataset(D2, req="LD"), error=identity),
               snp_ld_consistent=ok_snp),
          file.path(od, paste0(L,"_QC.rds")))

  ## regional minimum OUTCOME p (is there anything to fine-map?)
  z2 <- D2$beta/sqrt(D2$varbeta); p2 <- 2*pnorm(-abs(z2))
  min_p_out <- min(p2, na.rm=TRUE)

  ## fit SuSiE per trait — REPRODUCE the original call (defaults L=10, coverage=0.95,
  ## repeat_until_convergence=FALSE), exactly as run_susie.R used.
  S1 <- tryCatch(runsusie(D1, maxit=100, repeat_until_convergence=FALSE), error=function(e) e)
  S2 <- tryCatch(runsusie(D2, maxit=100, repeat_until_convergence=FALSE), error=function(e) e)
  saveRDS(list(S1=S1, S2=S2, coloc_version=as.character(packageVersion("coloc")),
               susieR_version=as.character(packageVersion("susieR"))),
          file.path(od, paste0(L,"_susie_objects.rds")))

  exp_ok <- !inherits(S1,"error"); out_ok <- !inherits(S2,"error")
  if (!exp_ok) writeLines(conditionMessage(S1), file.path(od,paste0(L,"_exposure_ERROR.txt")))
  if (!out_ok) writeLines(conditionMessage(S2), file.path(od,paste0(L,"_outcome_ERROR.txt")))
  if (exp_ok) try(capture.output(summary(S1), file=file.path(od,paste0(L,"_exposure_summary.txt"))), silent=TRUE)
  if (out_ok) try(capture.output(summary(S2), file=file.path(od,paste0(L,"_outcome_summary.txt"))), silent=TRUE)
  if (exp_ok) fwrite(cs_table(S1), file.path(od,paste0(L,"_exposure_CS.tsv")), sep="\t")
  if (out_ok) fwrite(cs_table(S2), file.path(od,paste0(L,"_outcome_CS.tsv")), sep="\t")

  ## --- B vs C adjudication for the OUTCOME ---
  oc_default  <- if (out_ok) cs_at(S2, R, 0.95, 0.5) else NA_integer_  # as reported
  oc_cov50    <- if (out_ok) cs_at(S2, R, 0.50, 0.5) else NA_integer_  # weaker coverage
  oc_relaxpur <- if (out_ok) cs_at(S2, R, 0.95, 0.1) else NA_integer_  # relaxed purity
  ec_default  <- if (exp_ok) cs_at(S1, R, 0.95, 0.5) else NA_integer_
  ## number of susie effects that "fired" (max per-effect PIP > 0.10) pre-purity
  eff_fired_out <- if (out_ok && !is.null(S2$alpha)) sum(apply(S2$alpha,1,max) > 0.10) else NA_integer_

  fwrite(data.table(coverage=c(0.95,0.50,0.95), min_abs_corr=c(0.5,0.5,0.1),
                    outcome_CS=c(oc_default,oc_cov50,oc_relaxpur)),
         file.path(od, paste0(L,"_outcome_CS_sensitivity.tsv")), sep="\t")

  ## --- coloc.susie (only meaningful if both have >=1 CS) ---
  n_pairs <- 0L; maxH4 <- NA_real_; maxH3 <- NA_real_; method <- "abf"
  if (exp_ok && out_ok){
    cs <- tryCatch(coloc.susie(S1,S2), error=function(e) NULL)
    if (!is.null(cs) && !is.null(cs$summary) && nrow(as.data.table(cs$summary))>0){
      su <- as.data.table(cs$summary)
      fwrite(su, file.path(od,paste0(L,"_coloc_susie_summary.tsv")), sep="\t")
      if (!is.null(cs$results)) fwrite(as.data.table(cs$results),
             file.path(od,paste0(L,"_coloc_susie_results.tsv")), sep="\t")
      n_pairs <- nrow(su); maxH4 <- max(su$PP.H4.abf,na.rm=TRUE)
      maxH3 <- if ("PP.H3.abf" %in% names(su)) max(su$PP.H3.abf,na.rm=TRUE) else NA_real_
      method <- "susie"
    }
  }
  ## coloc.abf fallback value (what the table currently reports)
  rabf <- tryCatch(coloc.abf(D1,D2), error=function(e) NULL)
  abf_H4 <- if (is.null(rabf)) NA_real_ else as.numeric(rabf$summary["PP.H4.abf"])
  abf_H3 <- if (is.null(rabf)) NA_real_ else as.numeric(rabf$summary["PP.H3.abf"])

  ## verdict
  verdict <- if (!exp_ok) "exposure SuSiE error -> abf only"
    else if (ec_default==0) "exposure has 0 CS -> abf only"
    else if (!out_ok) "outcome SuSiE error -> abf only"
    else if (oc_default>=1 && n_pairs>=1) "A: outcome CS present -> coloc.susie EVALUABLE"
    else if (oc_default==0 && oc_cov50==0 && oc_relaxpur==0) "B: no outcome CS at any threshold -> NOT evaluable (abf retained)"
    else "C: outcome CS only under relaxed coverage/purity -> borderline, interpret cautiously"

  row <- data.table(locus=L,
                    exposure_CS=ec_default, outcome_CS=oc_default,
                    outcome_CS_cov50=oc_cov50, outcome_CS_relaxpurity=oc_relaxpur,
                    outcome_effects_fired=eff_fired_out,
                    regional_min_p_outcome=signif(min_p_out,3),
                    n_signal_pairs=n_pairs,
                    coloc_susie_maxPP_H4=signif(maxH4,3), coloc_susie_maxPP_H3=signif(maxH3,3),
                    coloc_abf_PP_H4=signif(abf_H4,3), coloc_abf_PP_H3=signif(abf_H3,3),
                    method=method, verdict=verdict)
  fwrite(row, file.path(od, paste0(L,"_row.tsv")), sep="\t")
  print(row); cat("\n")
  rows[[L]] <- row
}

tab <- rbindlist(rows, fill=TRUE)
fwrite(tab, file.path(OUTDIR, "TableS1_rebuilt.tsv"), sep="\t")
cat("\n=================  REBUILT TABLE S1 (", OUTDIR, "/TableS1_rebuilt.tsv )  =================\n", sep="")
print(tab)
cat("\nPer-locus folders contain: *_susie_objects.rds, *_exposure_CS.tsv, *_outcome_CS.tsv,\n",
    "*_outcome_CS_sensitivity.tsv, *_exposure_summary.txt, *_outcome_summary.txt,\n",
    "*_coloc_susie_summary.tsv/_results.tsv (if evaluable), *_QC.rds, *_row.tsv\n", sep="")
cat("\nSEND BACK to finalise Table S1:  TableS1_rebuilt.tsv  +  every */*_outcome_CS.tsv\n",
    "                                 +  every */*_outcome_CS_sensitivity.tsv\n", sep="")
