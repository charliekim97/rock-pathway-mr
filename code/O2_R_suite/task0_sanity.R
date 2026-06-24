## Task 0 — pre-submission sanity checks (run FIRST). v3: 0a split into internal + direction.
## Not a re-run of the main result; integrity checks on existing conversion + harmonisation logic.
source("R/config.R"); source("R/helpers.R"); source("R/eqtl_catalogue.R")
options(width=120)

cat("==== 0a-1: Z->beta/se internal round-trip (eQTLGen) ====\n")
e <- fread(file.path(path.expand(CONFIG$data),"rock2_eqtlgen_converted.txt"))
## property of Zhu conversion: beta/se == Z exactly
e[, Z_recon := beta/se]
maxerr <- e[, max(abs(Z_recon - Zscore))]
cat("max |beta/se - Zscore| =", signif(maxerr,3), if (maxerr<1e-6) "  PASS\n" else "  FAIL (conversion bug)\n")
## recompute beta from Z, eaf, n and compare
e[, beta_re := Zscore/sqrt(2*eaf*(1-eaf)*(n+Zscore^2))]
cat("max |beta - beta_recompute| =", signif(e[,max(abs(beta-beta_re))],3),"  (should be ~0)\n")
## allele/EAF preserved (no swaps / NA)
cat("EAF in (0,1):", e[, all(eaf>0 & eaf<1)], " | alleles non-missing:", e[, all(AssessedAllele!="" & OtherAllele!="")],"\n")

cat("\n==== 0a-2: cross-study DIRECTION concordance (eQTLGen blood vs GTEx muscle) — direction only ====\n")
ds <- ec_dataset_id("GTEx","muscle")
if (!is.na(ds)){
  g <- ec_fetch_gene(ds, "ENSG00000134318")
  m <- merge(e[,.(SNP, EA=AssessedAllele, OA=OtherAllele, b_blood=beta)],
             g[,.(SNP, EA_g=EA, OA_g=OA, slope=beta)], by="SNP")
  ## align allele
  m[, slope_al := ifelse(EA_g==EA, slope, ifelse(EA_g==OA, -slope, NA))]
  m <- m[!is.na(slope_al)]
  conc <- m[, mean(sign(b_blood)==sign(slope_al))]
  cat("common SNP:", nrow(m), " | sign concordance:", round(conc,3),
      "\n  (magnitude differs by tissue/normalization — DO NOT use size as pass/fail)\n")
  pdf(file.path(CONFIG$out,"task0_blood_vs_muscle_direction.pdf"),5,5)
  plot(m$b_blood, m$slope_al, pch=20, col="#00000060", xlab="eQTLGen blood beta", ylab="GTEx muscle slope (aligned)",
       main="ROCK2 cross-study direction (sanity)"); abline(h=0,v=0,lty=3); dev.off()
} else cat("GTEx muscle dataset not resolved; skip direction check.\n")

cat("\n==== 0b: harmonisation / build log check (manual review prompts) ====\n")
cat(" - Confirm <100% allele-mismatch in each outcome harmonise (100% => build error).\n",
    "- Palindromic rs4335920: confirm EAF ~0.34 on BOTH sides (logged), clean resolve.\n",
    "- All 17 frozen instruments matched by rsID (no position fallback).\n")
cat("\n==== 0c: coloc input SNP counts -> see task5 regional table (nSNP column). ====\n")
cat("\nIf 0a PASS and 0b logs clean -> proceed to tasks 1..11.\n")
