## eqtl_catalogue.R — fetch GTEx-uniform cis-eQTL from eQTL Catalogue REST API v2.
## Avoids the 462GB GTEx download. AC/AN swap bug is fixed in current files.
## deps: jsonlite
if (!requireNamespace("jsonlite", quietly=TRUE)) stop("install.packages('jsonlite')")
suppressMessages(library(jsonlite))

## find dataset_id for a study/tissue (quant_method 'ge' = gene expression)
ec_dataset_id <- function(study="GTEx", tissue, quant="ge"){
  u <- sprintf("https://www.ebi.ac.uk/eqtl/api/v2/datasets/?study_label=%s&quant_method=%s&size=1000", study, quant)
  d <- fromJSON(u)
  hit <- d[grepl(tissue, d$sample_group, ignore.case=TRUE) | grepl(tissue, d$tissue_label, ignore.case=TRUE),]
  if (!nrow(hit)) return(NA_character_)
  hit$dataset_id[1]
}

## fetch all cis associations for a gene (molecular_trait_id = ENSG) in a dataset
ec_fetch_gene <- function(dataset_id, ensg, max_pages=40){
  base <- sprintf("https://www.ebi.ac.uk/eqtl/api/v2/datasets/%s/associations?molecular_trait_id=%s&size=1000",
                  dataset_id, ensg)
  out <- list(); start <- 0
  for (i in seq_len(max_pages)){
    u <- paste0(base, "&start=", start)
    pg <- tryCatch(fromJSON(u), error=function(e) NULL)
    if (is.null(pg) || length(pg)==0 || (is.data.frame(pg) && nrow(pg)==0)) break
    out[[length(out)+1]] <- pg
    if (is.data.frame(pg) && nrow(pg) < 1000) break
    start <- start + 1000
  }
  if (!length(out)) return(data.table())
  d <- rbindlist(lapply(out, as.data.table), fill=TRUE)
  # columns: rsid, beta, se, position, ref, alt, maf, pvalue, an, ac ...
  d <- d[!is.na(beta) & !is.na(se) & se>0]
  d[, .(SNP=rsid, EA=alt, OA=ref, beta=as.numeric(beta), se=as.numeric(se),
        eaf=as.numeric(maf), p=as.numeric(pvalue))][SNP %like% "^rs"]
}
message("eqtl_catalogue loaded.")
