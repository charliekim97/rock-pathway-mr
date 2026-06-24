# Data

`processed/` holds small, reproducibility-grade intermediate files committed to the repo:
- `instruments_frozen.csv` — frozen instruments (gene, SNP, tier, tissue, alleles, β, se, eaf)
- `eqtl_region_converted.csv` — ROCK2/KCNJ11/SORT1/GIPR cis-region eQTL (alleles, β/se)
- `rho_eqtl_region_conv.csv`, `rho_region_full_eqtl.txt` — ROCK1/RHOA cis-region eQTL
- `rock2_eqtlgen_converted.txt` — ROCK2 eQTLGen Z→β/se conversion
- `region_full_eqtl.txt` — full ROCK2 cis-region eQTL (with alleles)
- `finngen_*_region.tsv`, `nafld_meta_region.tsv`, `fibrosis_region_snps.txt` — outcome cis-region extracts
- `susie_inputs/` — per-locus coloc.susie inputs (z/β + 1000G-EUR LD + meta) for the 8 loci

## Raw summary statistics (NOT in repo — download from source)
Too large to host; fetch from the original providers:

| Dataset | Source |
|---|---|
| eQTLGen cis-eQTL (blood) | https://www.eqtlgen.org |
| GTEx v8 / eQTL Catalogue (muscle, liver, adipose) | https://www.ebi.ac.uk/eqtl |
| MAGIC glycemic traits | https://magicinvestigators.org |
| Pulit 2019 BMI (GIANT+UKB) | Zenodo 1251813 |
| GLGC LDL (Graham 2021) | GWAS Catalog GCST90239658 |
| FinnGen R12 (cirrhosis, K74, NASH) | https://www.finngen.fi |
| DIAMANTE T2D (Mahajan 2022; EUR/TA) | https://diagram-consortium.org / Broad portal |
| Ghodsian / Chen NAFLD | GWAS Catalog GCST90091033 / GCST90271622 |
| Genebass (UKB exome) | https://app.genebass.org |
| 1000 Genomes phase 3 (LD reference) | https://www.internationalgenome.org |

All builds/accessions and processing steps are described in the manuscript Methods.
