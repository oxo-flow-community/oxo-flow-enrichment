# load libraries
library("GenomicRanges")
library("rGREAT")
library("data.table")
library("rtracklayer")

# functions
annot_terms_with_features <- function(res, df) {
  regions <- vector("character", nrow(df))
  genes <- vector("character", nrow(df))

  map_associated_regions <- as.integer(commandArgs(trailingOnly = TRUE)[12])
  adjp_col <- commandArgs(trailingOnly = TRUE)[13]
  adjp_th <- as.numeric(commandArgs(trailingOnly = TRUE)[14])
  needs_annotation <- rep(FALSE, nrow(df))
  if (map_associated_regions != 0 && adjp_col %in% names(df)) {
    significant <- which(!is.na(df[[adjp_col]]) & df[[adjp_col]] <= adjp_th)
    if (map_associated_regions == -1) {
      needs_annotation[significant] <- TRUE
    } else if (map_associated_regions > 0 && length(significant) > 0) {
      top_terms <- significant[order(df[[adjp_col]][significant])][seq_len(min(map_associated_regions, length(significant)))]
      needs_annotation[top_terms] <- TRUE
    }
  }


  term_ids <- unique(df$id[needs_annotation])
  annotations <- setNames(vector("list", length(term_ids)), term_ids)

  for (term in term_ids) {
    gr <- getRegionGeneAssociations(res, term_id = term)
    annotations[[term]] <- list(
      regions = paste(
        paste0(seqnames(gr), ":", pmax(0L, start(gr) - 1L), "-", end(gr)), # export in BED-style coordinates: 0-based, start-inclusive, end-exclusive
        collapse = ","
      ),
      annotated_genes = paste(unique(unlist(gr$annotated_genes)), collapse = ",")
    )
  }

  for (i in which(needs_annotation)) {
    annotation <- annotations[[df$id[i]]]
    if (!is.null(annotation)) {
      regions[i] <- annotation$regions
      genes[i] <- annotation$annotated_genes
    }
  }

  df$regions <- regions
  df$annotated_genes <- genes

  return(df)
}

# Port of upstream workflow/scripts/region_enrichment_analysis_GREAT.R.
# snakemake@input[["regions"|"background"|"database"]],
# snakemake@output[["result"]], snakemake@config[["genome"]], snakemake@threads
# and the rule params become positional CLI arguments:
# <regions> <background> <database> <result> <genome> <cores>
# <min_gene_set_size> <mode> <basal_upstream> <basal_downstream> <extension>
# <map_associated_regions> <adjp_col> <adjp_th>

#input
regions_file <- commandArgs(trailingOnly = TRUE)[1]
background_file <- commandArgs(trailingOnly = TRUE)[2]
database_path <- commandArgs(trailingOnly = TRUE)[3]

# output
result_path <- commandArgs(trailingOnly = TRUE)[4]

# parameters
genome <- commandArgs(trailingOnly = TRUE)[5]
cores_n <- as.integer(commandArgs(trailingOnly = TRUE)[6])

# set genome
if (genome=="hg19" | genome=="hg38"){
    orgdb <- "org.Hs.eg.db"
}else if(genome=="mm9" | genome=="mm10"){
    orgdb <- "org.Mm.eg.db"
}

# stop early for empty query or background input
if (file.size(regions_file) == 0L || file.size(background_file) == 0L){
    file.create(result_path)
    quit(save = "no", status = 0)
}

# load query and background/universe region sets (e.g., consensus region set)
regionSet_query <- import(regions_file, format = "BED")
regionSet_background <- import(background_file, format = "BED")

# load database
database = read_gmt(file.path(database_path), from = "SYMBOL", to = "ENTREZ", orgdb = orgdb)

###### GREAT

# run GREAT
res <- great(gr = regionSet_query,
      gene_sets = database,
      tss_source = genome,
      biomart_dataset = NULL,
      min_gene_set_size = as.integer(commandArgs(trailingOnly = TRUE)[7]), #default: 5
      mode = commandArgs(trailingOnly = TRUE)[8],
      basal_upstream = as.numeric(commandArgs(trailingOnly = TRUE)[9]),
      basal_downstream = as.numeric(commandArgs(trailingOnly = TRUE)[10]),
      extension = as.numeric(commandArgs(trailingOnly = TRUE)[11]),
      extended_tss = NULL,
      background = regionSet_background, #default: NULL
      exclude = "gap",
      cores = cores_n, #default: 1
      verbose = TRUE #default: great_opt$verbose
     )

# get & save result table
tb <- getEnrichmentTable(res, min_region_hits = 0)
tb$description <- paste(tb$description, tb$id)
# annotate (near) significant enrichments with features
tb <- annot_terms_with_features(res, tb)
fwrite(as.data.frame(tb), file=file.path(result_path), row.names=FALSE)
