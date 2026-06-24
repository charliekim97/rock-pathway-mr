## Task 5 (v3 reduced) — regional descriptive table + locuszoom. NO permutation.
## Distinguishes "no outcome signal in region" (coloc uninformative) from "signal present but different variant".
## Per locus: regional SNP count, regional min outcome p, LD r2(eQTL lead, outcome lead), coloc result.
source("R/config.R"); source("R/helpers.R")
suppressMessages({library(ggplot2)})

ld_r2 <- function(a,b){ bf<-path.expand(CONFIG$ld_plink_bfile); tmp<-tempfile()
  system(sprintf("%s --bfile %s --ld %s %s --out %s", CONFIG$plink_bin, bf, a, b, tmp),
         ignore.stdout=TRUE, ignore.stderr=TRUE)
  l <- tryCatch(readLines(paste0(tmp,".log")), error=function(e) "")
  m <- regmatches(l, regexpr("R-sq = [0-9.]+", l)); if(length(m)) as.numeric(sub("R-sq = ","",m[1])) else NA }

regional_row <- function(e, od, eqtl_lead, label){
  m <- merge(unique(e,by="SNP"), unique(od,by="SNP"), by="SNP")
  if (!nrow(m) || !"p" %in% names(od)) return(data.table(locus=label, nSNP=nrow(m), min_p=NA, lead_LD=NA))
  out_lead <- od[SNP %in% m$SNP][which.min(p)]$SNP
  data.table(locus=label, nSNP=nrow(m), min_p=min(od[SNP %in% m$SNP]$p, na.rm=TRUE),
             eqtl_lead=eqtl_lead, out_lead=out_lead, lead_LD_r2=ld_r2(eqtl_lead, out_lead))
}
rows <- rbindlist(list(
  regional_row(load_eqtl_region("ROCK2"), read_pulit_bmi(CONFIG$pulit_bmi),        "rs12468344","ROCK2_BMI"),
  regional_row(load_eqtl_region("ROCK2"), read_harmonised(CONFIG$ghodsian_nafld),  "rs12468344","ROCK2_NAFLD"),
  regional_row(load_eqtl_region("ROCK2"), read_finngen(CONFIG$finngen_cirr),       "rs12468344","ROCK2_cirrhosis"),
  regional_row(load_eqtl_region("SORT1"), read_harmonised(CONFIG$glgc_ldl),        "rs657420","SORT1_LDL")
), fill=TRUE)
fwrite(rows, file.path(CONFIG$out,"task5_regional_table.csv")); print(rows)
cat("\nInterpretation: min_p > 1e-5 => no outcome signal in region (coloc uninformative).",
    "\n min_p genome-wide-sig AND lead_LD_r2 low => signal present but distinct from eQTL (true no-coloc).\n")

## locuscompare-style panel (eQTL -log10p vs outcome -log10p), one example (ROCK2 x cirrhosis)
e <- load_eqtl_region("ROCK2"); e[, p_eqtl := p_from_z(beta/se)]
o <- read_finngen(CONFIG$finngen_cirr)[, .(SNP,p)]
mm <- merge(e[,.(SNP,p_eqtl)], o, by="SNP")
g <- ggplot(mm, aes(-log10(p_eqtl), -log10(p)))+geom_point(alpha=.5,size=.8)+theme_bw(base_size=11)+
  labs(x="eQTL  -log10 P (ROCK2 blood)", y="Cirrhosis  -log10 P", title="ROCK2 region: eQTL vs cirrhosis")
ggsave(file.path(CONFIG$out,"task5_locuscompare_ROCK2_cirrhosis.png"), g, width=5, height=5, dpi=300)
