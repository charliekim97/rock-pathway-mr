## ============================================================
## Task 4 — SuSiE conditional colocalization (relaxes the single-causal-variant
## assumption of coloc.abf). Run on O2/HPC where susieR + coloc + plink are available.
##
## Loci (5 + 1 control; RHOA excluded — single weak instrument):
##   1 ROCK2-BMI          (eQTL ROCK2 blood  × Pulit BMI)
##   2 ROCK2-NAFLD        (eQTL ROCK2 blood  × Ghodsian NAFLD, cc)
##   3 ROCK2-cirrhosis    (eQTL ROCK2 blood  × FinnGen cirrhosis, cc)
##   4 ROCK1-hepatic fat  (eQTL ROCK1 blood  × Haas hepatic fat)
##   5 SORT1-LDL          (eQTL SORT1 *liver* × GLGC LDL dense)   <- positive control (multi-signal)
##   6 KCNJ11-T2D         (eQTL KCNJ11 blood × DIAMANTE-EUR T2D, cc) <- positive control (sensitivity)
##
## EXPECTATION (defends the expression-null + justifies restoring "robust"):
##   loci 1-4 -> low PP.H4 at every signal pair (no colocalization under multi-signal model),
##              matching coloc.abf; loci 5-6 -> high PP.H4 (method recovers true shared signals).
##
## Output order (packet B3): (1) PP.H4 per signal pair  (2) credible-set size
##   (3) top-PIP SNP  (4) overlap. Falls back to coloc.abf if SuSiE does not converge.
##
## ---- FILES TO STAGE on O2 (edit paths in config.R) ----
##   data/eqtl_region_converted.csv   (ROCK2/KCNJ11/SORT1 eQTL: SNP,SNPChr,SNPPos,AssessedAllele,OtherAllele,eaf,NrSamples,Zscore,beta,se,varbeta,GeneSymbol)
##   data/rho_eqtl_region_conv.csv    (ROCK1/RHOA eQTL beta/se)
##   data/rho_region_full_eqtl.txt    (ROCK1/RHOA eQTL WITH alleles; CONFIG$rho_full)
##   DIAMANTE-EUR, Ghodsian NAFLD, FinnGen cirrhosis, Haas hepatic fat, GLGC LDL, Pulit BMI (CONFIG paths)
##   For SORT1-LDL liver: place a LIVER SORT1 eQTL file as eqtl_region_converted.csv GeneSymbol=="SORT1_liver"
##     (GTEx liver via eQTL Catalogue; columns as above). If only blood SORT1 is staged the control will be
##     mis-tissued (PP4~0) — use liver, per the main-text finding.
##   ref/1000G_EUR.{bed,bim,fam}      (CONFIG$ld_plink_bfile)  + plink 1.9 on PATH
## ============================================================
source("R/config.R"); source("R/helpers.R")
if (!requireNamespace("susieR", quietly=TRUE)) stop("install.packages('susieR')")
if (!requireNamespace("coloc",  quietly=TRUE)) stop("install.packages('coloc')")
suppressMessages({library(susieR); library(coloc); library(data.table)})

## ---- LD: signed correlation matrix + A1 alleles, from 1000G EUR plink bfile ----
## Returns list(R=matrix[snps x snps], a1=named vector of plink A1 allele).
get_ld <- function(snps){
  bf <- path.expand(CONFIG$ld_plink_bfile); tmp <- tempfile()
  writeLines(unique(snps), paste0(tmp,".snps"))
  system(sprintf("%s --bfile %s --extract %s --r square --write-snplist --out %s",
                 CONFIG$plink_bin, bf, paste0(tmp,".snps"), tmp),
         ignore.stdout=TRUE, ignore.stderr=TRUE)
  if (!file.exists(paste0(tmp,".ld"))) return(NULL)
  kept <- readLines(paste0(tmp,".snplist"))
  M <- as.matrix(fread(paste0(tmp,".ld"))); dimnames(M) <- list(kept, kept)
  bim <- fread(paste0(bf,".bim"), header=FALSE)          # V2=SNP, V5=A1, V6=A2
  a1  <- setNames(toupper(bim$V5), bim$V2)[kept]
  list(R=M, a1=a1)
}

## ---- one locus: allele-harmonise eQTL & outcome to LD A1, run SuSiE on both, coloc.susie ----
susie_pair <- function(eq, od, type2, s=NULL, N2, Neqtl=31684, label){
  ## eq and od each carry columns: SNP, EA, OA, beta, se (readers in helpers.R)
  m <- merge(unique(eq, by="SNP"), unique(od[,.(SNP,EA,OA,beta,se)], by="SNP"),
             by="SNP", suffixes=c(".e",".o"))
  m <- m[se.e>0 & se.o>0 & is.finite(beta.e) & is.finite(beta.o)]
  ## harmonise outcome to eQTL effect allele (flip sign if swapped; drop mismatches)
  same <- toupper(m$EA.o)==toupper(m$EA.e) & toupper(m$OA.o)==toupper(m$OA.e)
  swap <- toupper(m$EA.o)==toupper(m$OA.e) & toupper(m$OA.o)==toupper(m$EA.e)
  m <- m[same | swap]; m[swap, beta.o := -beta.o]
  if (nrow(m) < 20) return(list(summary=data.table(locus=label, note="too few SNP after harmonise", PP.H4=NA_real_)))
  ld <- get_ld(m$SNP); if (is.null(ld)) return(list(summary=data.table(locus=label, note="no LD", PP.H4=NA_real_)))
  keep <- intersect(m$SNP, rownames(ld$R)); m <- m[SNP %in% keep]
  R <- ld$R[m$SNP, m$SNP]; a1 <- ld$a1[m$SNP]
  ## align eQTL effect allele to LD A1 (flip BOTH traits' beta where eQTL EA != A1; LD then consistent)
  flip <- toupper(m$EA.e) != a1
  m[flip, `:=`(beta.e=-beta.e, beta.o=-beta.o)]
  z1 <- m$beta.e/m$se.e; z2 <- m$beta.o/m$se.o
  D1 <- list(beta=m$beta.e, varbeta=m$se.e^2, snp=m$SNP, LD=R, type="quant", sdY=1, N=Neqtl)
  D2 <- list(beta=m$beta.o, varbeta=m$se.o^2, snp=m$SNP, LD=R, type=type2, N=N2)
  if (type2=="cc") D2$s <- s else D2$sdY <- 1
  out <- tryCatch({
    S1 <- runsusie(D1); S2 <- runsusie(D2)
    cs <- coloc.susie(S1, S2)
    cs1 <- if (!is.null(S1$sets$cs)) sapply(S1$sets$cs, length) else integer(0)
    cs2 <- if (!is.null(S2$sets$cs)) sapply(S2$sets$cs, length) else integer(0)
    top1 <- names(which.max(susie_get_pip(S1))); top2 <- names(which.max(susie_get_pip(S2)))
    su <- if (!is.null(cs$summary)) as.data.table(cs$summary) else
          data.table(idx1=NA,idx2=NA,PP.H4.abf=NA_real_,nsnps=nrow(m))
    su[, `:=`(locus=label, eqtl_CSsize=paste(cs1,collapse="/"), out_CSsize=paste(cs2,collapse="/"),
              eqtl_topPIP=top1, out_topPIP=top2, method="susie")]
    list(summary=su[order(-PP.H4.abf)])
  }, error=function(e){
    r <- run_coloc(dt1=m[,.(SNP,beta=beta.e,se=se.e)], dt2=m[,.(SNP,beta=beta.o,se=se.o)],
                   type2=type2, s=s, N2=N2)
    list(summary=data.table(locus=label, idx1=NA, idx2=NA, PP.H4.abf=as.numeric(r$PP4),
                            nsnps=r$nsnp, method=paste0("abf_fallback(",conditionMessage(e),")")))
  })
  out
}

## ---- wire loci ----
sort1_liver <- function(){ d <- fread(file.path(path.expand(CONFIG$data),"eqtl_region_converted.csv"))
  g <- if ("SORT1_liver" %in% d$GeneSymbol) "SORT1_liver" else "SORT1"
  d[GeneSymbol==g, .(SNP, EA=toupper(AssessedAllele), OA=toupper(OtherAllele), beta, se)] }

loci <- list(
  list("ROCK2-BMI",       load_eqtl_region2("ROCK2"),  read_pulit_bmi(CONFIG$pulit_bmi),     "quant", NULL, 806834),
  list("ROCK2-NAFLD",     load_eqtl_region2("ROCK2"),  read_harmonised(CONFIG$ghodsian_nafld),"cc", 8434/778614, 778614),
  list("ROCK2-cirrhosis", load_eqtl_region2("ROCK2"),  read_finngen(CONFIG$finngen_cirr),    "cc", 5545/500348, 500348),
  list("ROCK1-hepaticfat",load_eqtl_region2("ROCK1"),  read_haas(CONFIG$haas_hepfat),        "quant", NULL, 32974),
  list("SORT1-LDL(liver)",sort1_liver(),               read_harmonised(CONFIG$glgc_ldl),     "quant", NULL, 1320016),
  list("KCNJ11-T2D",      load_eqtl_region2("KCNJ11"), read_diamante(CONFIG$diamante_t2d),   "cc", 80154/933970, 933970)
)
res <- rbindlist(lapply(loci, function(L)
         susie_pair(L[[2]], L[[3]], L[[4]], L[[5]], L[[6]], label=L[[1]])$summary), fill=TRUE)
fwrite(res, file.path(path.expand(CONFIG$out),"task4_susie_coloc.csv")); print(res)
cat("\nInterpretation: report PP.H4 per signal pair FIRST; credible-set size, top-PIP SNP, then overlap.\n",
    "Expect ROCK loci low PP.H4 (no coloc under multi-signal); KCNJ11-T2D & SORT1-LDL(liver) high PP.H4 (controls).\n")
