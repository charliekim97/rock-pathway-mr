## ============================================================
## helpers.R — shared functions. source("config.R") first.
## deps: data.table, coloc, susieR (task4), ggplot2 (figs), R.utils (gz)
## ============================================================
need <- c("data.table","coloc")
for (p in need) if (!requireNamespace(p, quietly=TRUE)) stop("install.packages('",p,"')")
suppressMessages({library(data.table); library(coloc)})

z_from_p <- function(p){ p[p<=0] <- 1e-300; qnorm(1 - p/2) }
p_from_z <- function(z) 2*pnorm(-abs(z))

## ---- eQTLGen Z -> beta/se (Zhu 2016). n default 31684 ----
zhu_convert <- function(Z, eaf, n=31684){
  den <- sqrt(2*eaf*(1-eaf)*(n + Z^2))
  list(beta = Z/den, se = 1/den)
}

## ---- frozen instruments ----
load_instruments <- function() fread(file.path(path.expand(CONFIG$data),"instruments_frozen.csv"))
load_eqtl_region <- function(gene){
  f <- if (gene %in% c("ROCK1","RHOA")) "rho_eqtl_region_conv.csv" else "eqtl_region_converted.csv"
  d <- fread(file.path(path.expand(CONFIG$data), f))
  d[GeneSymbol==gene, .(SNP, beta, se)]
}

## ---- harmonise one outcome row to exposure allele (returns aligned by) ----
## exp: list(EA,OA,beta,se); out: list(EA,OA,beta,se). match by allele; flip if swapped.
harmonise <- function(exp, out){
  if (is.na(out$beta) || is.na(out$se) || out$se<=0) return(NULL)
  bo <- out$beta
  if (toupper(out$EA)==toupper(exp$EA) && toupper(out$OA)==toupper(exp$OA)) {}
  else if (toupper(out$EA)==toupper(exp$OA) && toupper(out$OA)==toupper(exp$EA)) bo <- -bo
  else return(NULL)        # allele mismatch (do not strand-flip blindly; log)
  list(bx=exp$beta, sx=exp$se, by=bo, sy=out$se)
}
mr_wald <- function(h){ b<-h$by/h$bx; se<-abs(h$sy/h$bx); c(b=b, se=se, p=p_from_z(b/se)) }
mr_ivw  <- function(hl){ # list of harmonised pairs
  hl <- Filter(Negate(is.null), hl); if (length(hl)<2) return(NULL)
  bx<-sapply(hl,`[[`,"bx"); by<-sapply(hl,`[[`,"by"); sy<-sapply(hl,`[[`,"sy")
  w<-bx^2/sy^2; b<-sum(bx*by/sy^2)/sum(w); se<-sqrt(1/sum(w))
  Q<-sum((by-b*bx)^2/sy^2); df<-length(hl)-1
  c(b=b, se=se, p=p_from_z(b/se), Q=Q, Qp=pchisq(Q,df,lower.tail=FALSE), nsnp=length(hl))
}

## ---- coloc.abf wrapper. d1=eQTL(quant), d2=outcome. type2 'cc' or 'quant'. ----
## dt1/dt2: data.table(SNP, beta, se [, MAF, N]). returns summary + PP4.
run_coloc <- function(dt1, dt2, type2="cc", s=NULL, sdY1=1, N2=NULL, p12=CONFIG$p12){
  m <- merge(unique(dt1,by="SNP"), unique(dt2,by="SNP"), by="SNP", suffixes=c("1","2"))
  m <- m[se1>0 & se2>0 & is.finite(beta1) & is.finite(beta2)]
  if (nrow(m) < 5) return(list(nsnp=nrow(m), PP3=NA, PP4=NA, note="too few SNP"))
  D1 <- list(beta=m$beta1, varbeta=m$se1^2, snp=m$SNP, type="quant", sdY=sdY1, N=31684)
  D2 <- list(beta=m$beta2, varbeta=m$se2^2, snp=m$SNP, type=type2)
  if (type2=="cc") D2$s <- s else D2$sdY <- 1
  if (!is.null(N2)) D2$N <- N2
  r <- suppressWarnings(coloc.abf(D1, D2, p1=CONFIG$p1, p2=CONFIG$p2, p12=p12))
  s <- r$summary
  list(nsnp=as.integer(s["nsnps"]), PP0=s["PP.H0.abf"], PP1=s["PP.H1.abf"],
       PP2=s["PP.H2.abf"], PP3=s["PP.H3.abf"], PP4=s["PP.H4.abf"], res=r)
}

## ---- outcome region readers (return data.table SNP,EA,OA,beta,se[,p,eaf]) ----
read_finngen <- function(path){   # hg38 cols: #chrom pos ref alt rsids ... pval mlogp beta sebeta af_alt
  d <- fread(cmd=paste("gzip -dc", shQuote(path.expand(path))))
  setnames(d, old=c("rsids","alt","ref","beta","sebeta","pval"), new=c("SNP","EA","OA","beta","se","p"), skip_absent=TRUE)
  d[SNP %like% "^rs", .(SNP, EA, OA, beta, se, p)]
}
read_harmonised <- function(path){ # GWAS Catalog harmonised: hm_rsid hm_effect_allele hm_other_allele hm_beta standard_error p_value
  d <- fread(cmd=paste("gzip -dc", shQuote(path.expand(path))))
  setnames(d, old=c("hm_rsid","hm_effect_allele","hm_other_allele","hm_beta","standard_error","p_value"),
           new=c("SNP","EA","OA","beta","se","p"), skip_absent=TRUE)
  d[!is.na(beta) & !is.na(se), .(SNP, EA, OA, beta, se, p)]
}
read_pulit_bmi <- function(path){  # SNP=rsid:A1:A2 Allele1 Allele2 Effect StdErr P-value
  d <- fread(cmd=paste("gzip -dc", shQuote(path.expand(path))))
  d[, SNP:=sub(":.*","",SNP)]
  d[, .(SNP, EA=toupper(Allele1), OA=toupper(Allele2), beta=Effect, se=StdErr, p=`P-value`)]
}
## generic: subset an outcome dt to instrument rsIDs and return named list per SNP
outcome_lookup <- function(out_dt, snp){ r<-out_dt[SNP==snp]; if(!nrow(r)) return(NULL)
  list(EA=r$EA[1], OA=r$OA[1], beta=r$beta[1], se=r$se[1], p=if("p"%in%names(r)) r$p[1] else NA) }

read_diamante <- function(path){ # DIAMANTE-EUR (Mahajan 2022, hg19), space-delim:
  # chromosome(b37) position(b37) chrposID rsID effect_allele other_allele EAF Fixed-effects_beta SE p
  d <- fread(cmd=paste("gzip -dc", shQuote(path.expand(path))))
  setnames(d, old=c("rsID","effect_allele","other_allele","Fixed-effects_beta","Fixed-effects_SE","Fixed-effects_p-value"),
           new=c("SNP","EA","OA","beta","se","p"), skip_absent=TRUE)
  d[SNP %like% "^rs" & !is.na(beta) & se>0, .(SNP, EA=toupper(EA), OA=toupper(OA), beta, se, p)]
}
read_haas <- function(path){ # Haas 2021 hepatic fat — GWAS Catalog harmonised (reuse harmonised reader)
  read_harmonised(path)
}
## allele-aware eQTL region loader: returns SNP,EA,OA,beta,se.
## ROCK2/KCNJ11/SORT1/GIPR from eqtl_region_converted.csv (has alleles+beta/se);
## ROCK1/RHOA from rho_eqtl_region_conv.csv (beta/se) merged with rho_region_full_eqtl.txt (alleles).
load_eqtl_region2 <- function(gene){
  if (gene %in% c("ROCK1","RHOA")) {
    bs <- fread(file.path(path.expand(CONFIG$data),"rho_eqtl_region_conv.csv"))[GeneSymbol==gene,.(SNP,beta,se)]
    al <- fread(path.expand(CONFIG$rho_full))[GeneSymbol==gene,.(SNP,EA=AssessedAllele,OA=OtherAllele)]
    merge(bs, al, by="SNP")[, .(SNP,EA=toupper(EA),OA=toupper(OA),beta,se)]
  } else {
    d <- fread(file.path(path.expand(CONFIG$data),"eqtl_region_converted.csv"))[GeneSymbol==gene]
    d[, .(SNP, EA=toupper(AssessedAllele), OA=toupper(OtherAllele), beta, se)]
  }
}

message("helpers loaded.")
