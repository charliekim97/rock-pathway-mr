## SuSiE conditional colocalization runner (O2). Reads pre-built z/beta + 1000G-EUR LD
## (allele-aligned in sandbox) and runs runsusie + coloc.susie per locus.
## Run on a COMPUTE node:  srun -t 0:30:00 --mem=6G -c 2 Rscript susie_inputs/run_susie.R
.libPaths("~/R/libs")
suppressMessages({library(susieR); library(coloc); library(data.table)})
DIR  <- path.expand("~/rock2_o2/susie_inputs")
ODIR <- path.expand("~/rock2_o2/out"); dir.create(ODIR, showWarnings=FALSE, recursive=TRUE)
loci <- sub("\\.betas\\.tsv$","", list.files(DIR, pattern="\\.betas\\.tsv$"))
cat("Loci:", paste(loci, collapse=", "), "\n\n")

top_pip <- function(S){ p <- tryCatch(susie_get_pip(S), error=function(e) NULL)
  if (is.null(p) || is.null(names(p)) || !length(p)) return(NA_character_); names(p)[which.max(p)] }
csizes  <- function(S){ if(!is.null(S$sets$cs) && length(S$sets$cs)) paste(sapply(S$sets$cs,length),collapse="/") else "0" }

run_one <- function(L){
  b  <- fread(file.path(DIR, paste0(L,".betas.tsv")))
  R  <- as.matrix(fread(file.path(DIR, paste0(L,".ld.gz")), header=FALSE))
  kv <- readLines(file.path(DIR, paste0(L,".meta"))); kv <- setNames(sub(".*=","",kv), sub("=.*","",kv))
  rownames(R) <- colnames(R) <- b$SNP
  D1 <- list(beta=b$b1, varbeta=b$vb1, snp=b$SNP, LD=R, type="quant", sdY=1, N=as.numeric(kv["Neqtl"]))
  D2 <- list(beta=b$b2, varbeta=b$vb2, snp=b$SNP, LD=R, type=unname(kv["type2"]), N=as.numeric(kv["N2"]))
  if (kv["type2"]=="cc") D2$s <- as.numeric(kv["s"]) else D2$sdY <- 1
  S1 <- tryCatch(runsusie(D1, maxit=100, repeat_until_convergence=FALSE), error=function(e) NULL)
  S2 <- tryCatch(runsusie(D2, maxit=100, repeat_until_convergence=FALSE), error=function(e) NULL)
  cs1 <- if(is.null(S1)) "NA" else csizes(S1); cs2 <- if(is.null(S2)) "NA" else csizes(S2)
  t1  <- if(is.null(S1)) NA_character_ else top_pip(S1); t2 <- if(is.null(S2)) NA_character_ else top_pip(S2)
  if (!is.null(S1) && !is.null(S2)) {
    cs <- tryCatch(coloc.susie(S1, S2), error=function(e) NULL)
    if (!is.null(cs) && !is.null(cs$summary) && nrow(as.data.table(cs$summary))>0) {
      su <- as.data.table(cs$summary)[order(-PP.H4.abf)]
      return(data.table(locus=L, idx1=su$idx1, idx2=su$idx2, PP.H4=su$PP.H4.abf, nsnps=su$nsnps,
                        eQTL_CS=cs1, out_CS=cs2, eQTL_topPIP=t1, out_topPIP=t2, method="susie"))
    }
    note <- if (cs2=="0") "susie_no_outcome_CS" else if (cs1=="0") "susie_no_eqtl_CS" else "susie_no_pair"
  } else note <- "susie_failed"
  ## fallback / no-pair: report coloc.abf PP.H4 for completeness
  r <- tryCatch(coloc.abf(D1, D2), error=function(e) NULL)
  data.table(locus=L, idx1=NA, idx2=NA,
             PP.H4=if(is.null(r)) NA_real_ else as.numeric(r$summary["PP.H4.abf"]),
             nsnps=nrow(b), eQTL_CS=cs1, out_CS=cs2, eQTL_topPIP=t1, out_topPIP=t2,
             method=paste0(note,"+abf"))
}
res <- rbindlist(lapply(loci, function(L){ cat("==", L, "==\n"); x <- run_one(L); print(x); flush.console(); x }), fill=TRUE)
fwrite(res, file.path(ODIR, "susie_coloc_results.csv"))
cat("\n=== SUMMARY (out/susie_coloc_results.csv) ===\n"); print(res)
