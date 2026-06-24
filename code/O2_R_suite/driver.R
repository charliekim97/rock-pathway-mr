## ============================================================
## driver.R — generic MR + coloc drivers. source config.R, helpers.R first.
## ============================================================

## plink-based LD proxy: returns proxy rsID (r2>thresh) present in out_dt, or NA
find_proxy <- function(lead, out_dt, r2=CONFIG$proxy_r2){
  bf <- path.expand(CONFIG$ld_plink_bfile)
  tmp <- tempfile()
  cmd <- sprintf("%s --bfile %s --ld-snp %s --r2 --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 %s --out %s",
                 CONFIG$plink_bin, bf, lead, r2, tmp)
  ok <- tryCatch(system(cmd, ignore.stdout=TRUE, ignore.stderr=TRUE)==0, error=function(e) FALSE)
  if (!ok || !file.exists(paste0(tmp,".ld"))) return(NA_character_)
  ld <- fread(paste0(tmp,".ld"))
  cand <- ld[SNP_B %in% out_dt$SNP][order(-R2)]
  if (nrow(cand)) cand$SNP_B[1] else NA_character_
}

## MR across frozen tiers for one outcome. out_dt: data.table(SNP,EA,OA,beta,se,p)
analyze_mr <- function(out_dt, label, allow_proxy=TRUE){
  inst <- load_instruments(); res <- list()
  pull <- function(snp){
    o <- outcome_lookup(out_dt, snp)
    if (is.null(o) && allow_proxy){ px <- find_proxy(snp, out_dt)
      if (!is.na(px)){ o <- outcome_lookup(out_dt, px); attr(o,"proxy") <- px } }
    o }
  exp_row <- function(snp){ r<-inst[SNP==snp]; list(EA=r$EA,OA=r$OA,beta=r$beta,se=r$se) }
  wald1 <- function(snp,tag){ o<-pull(snp); if(is.null(o)) return(NULL)
    h<-harmonise(exp_row(snp),o); if(is.null(h)) return(NULL)
    w<-mr_wald(h); data.frame(outcome=label,test=tag,snp=snp,proxy=ifelse(is.null(attr(o,"proxy")),NA,attr(o,"proxy")),
                              b=w["b"],se=w["se"],p=w["p"],nsnp=1) }
  ivwset <- function(snps,tag){ hl<-lapply(snps,function(s){o<-pull(s);if(is.null(o))return(NULL);harmonise(exp_row(s),o)})
    iv<-mr_ivw(hl); if(is.null(iv))return(NULL)
    data.frame(outcome=label,test=tag,snp=paste(snps,collapse=";"),proxy=NA,b=iv["b"],se=iv["se"],p=iv["p"],nsnp=iv["nsnp"]) }
  res[["r2lead"]]   <- wald1("rs12468344","ROCK2_blood_lead_Wald")
  res[["r2muscle"]] <- ivwset(c("rs4335920","rs12996712"),"ROCK2_muscle_IVW")
  res[["r2set"]]    <- ivwset(c("rs12468344","rs13034888","rs7574150"),"ROCK2_blood_r2.01_IVW")
  res[["rock1"]]    <- ivwset(c("rs142716063","rs72881399"),"ROCK1_blood_IVW")
  res[["rock1l"]]   <- wald1("rs142716063","ROCK1_blood_lead_Wald")
  res[["rhoa"]]     <- wald1("rs115431681","RHOA_blood_Wald")
  do.call(rbind, res)
}

## coloc for one outcome across genes. out_dt needs beta,se (allele-agnostic for coloc).
analyze_coloc <- function(out_dt, label, type2="cc", s=NULL, N2=NULL,
                          genes=c("ROCK2","ROCK1","RHOA"), p12=CONFIG$p12){
  od <- unique(out_dt[SNP %like% "^rs", .(SNP, beta, se)], by="SNP")
  rbindlist(lapply(genes, function(g){
    e <- load_eqtl_region(g); if(!nrow(e)) return(NULL)
    r <- run_coloc(e, od, type2=type2, s=s, N2=N2, p12=p12)
    data.table(outcome=label, gene=g, nsnp=r$nsnp, PP3=as.numeric(r$PP3), PP4=as.numeric(r$PP4))
  }), fill=TRUE)
}
message("driver loaded.")
