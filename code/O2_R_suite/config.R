## ============================================================
## config.R  — ROCK pathway MR O2 reinforcement suite (v3)
## Edit the paths in CONFIG, then source this from every task script.
## Principle: rsID matching across builds (eQTLGen hg19, GTEx/FinnGen/GWAS-Catalog-harmonised hg38).
## ============================================================
CONFIG <- list(
  base   = "~/rock2_o2",                 # working root on O2 (edit)
  data   = "~/rock2_o2/data",            # staged frozen data (instruments_frozen.csv, eqtl_region_converted.csv, rho_eqtl_region_conv.csv)
  out    = "~/rock2_o2/out",             # outputs
  # ---- LD reference: 1000G phase3 EUR plink bfiles (prefix, no extension) ----
  ld_plink_bfile = "~/ref/1000G_EUR",    # EDIT: plink .bed/.bim/.fam prefix
  plink_bin      = "plink",              # plink 1.9 on PATH (for clumping/proxies/SuSiE LD)
  # ---- outcome summary-stat files (download or place; see README) ----
  diamante_t2d   = "~/rock2_o2/data/DIAMANTE_EUR_T2D.txt.gz",      # Mahajan 2022 EUR (access form; place file)
  diamante_build = "hg19",                                         # CONFIRM from readme
  chen_nafld     = "~/rock2_o2/data/GCST90271622.h.tsv.gz",        # Chen 2023 NAFLD (GWAS Catalog harmonised, hg38)
  ghodsian_nafld = "~/rock2_o2/data/GCST90091033.h.tsv.gz",        # Ghodsian 2021 (existing)
  glgc_ldl       = "~/rock2_o2/data/GCST90239658.h.tsv.gz",        # GLGC Graham 2021 LDL EUR dense
  finngen_cirr   = "~/rock2_o2/data/finngen_R12_CIRRHOSIS_BROAD.gz",
  pulit_bmi      = "~/rock2_o2/data/bmi_pulit.tbl.gz",
  haas_hepfat    = "~/rock2_o2/data/GCST90029073.h.tsv.gz",        # Haas 2021 hepatic fat (GWAS Catalog harmonised)
  rho_full       = "~/rock2_o2/data/rho_region_full_eqtl.txt",     # ROCK1/RHOA eQTL region WITH alleles (AssessedAllele/OtherAllele)
  # ---- coloc priors ----
  p1=1e-4, p2=1e-4, p12=1e-5, coloc_pp4_threshold=0.75,
  win_kb = 500,                          # +-500kb cis window for coloc/regional
  proxy_r2 = 0.8                         # LD proxy threshold when lead absent
)
## gene coordinates (for region extraction). ROCK2/ROCK1/RHOA TSS approx.
GENES <- data.frame(
  gene = c("ROCK2","ROCK1","RHOA","SORT1","KCNJ11"),
  chr  = c(2,18,3,1,11),
  hg19 = c(11183000,18529000,49396000,109822000,17409000),
  hg38 = c(11128000,21272000,49359000,109274000,17387000),
  ensg = c("ENSG00000134318","ENSG00000067900","ENSG00000067560","ENSG00000134243","ENSG00000187486"),
  stringsAsFactors=FALSE)
dir.create(path.expand(CONFIG$out), showWarnings=FALSE, recursive=TRUE)
message("config loaded. out = ", CONFIG$out)
