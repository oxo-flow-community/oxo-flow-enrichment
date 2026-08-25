# Port of upstream workflow/scripts/gene_enrichment_analysis_RcisTarget.R.
# snakemake@input[["genes"|"background_genes"|"database"|"motif2tf"]],
# snakemake@output[["result"]], snakemake@wildcards[["gene_set"]],
# snakemake@config[["rcistarget_parameters"]] and snakemake@threads become
# positional CLI arguments:
# <genes> <background_genes> <database> <motif2tf> <result> <gene_set>
# <motifAnnot_highConfCat(comma-joined)> <motifAnnot_lowConfCat(comma-joined)>
# <nesThreshold> <aucMaxRank_factor> <geneErnMethod> <geneErnMaxRank> <cores>
# RcisTarget's own empty-result handling is kept verbatim: it creates the
# output file and exits 0 (upstream behavior).

# load libraries
library("RcisTarget")
library("data.table")

process_genes <- function(input) {
  # Split the input string into individual gene names
  genes <- unlist(strsplit(input, "; "))

  # Remove bracketed text ending with a dot
  genes <- gsub("\\s*\\([^\\)]+\\)\\.", "", genes)

  # Remove duplicate gene names
  unique_genes <- unique(genes)

  # Join the unique gene names back into a single string
  output <- paste(unique_genes, collapse = "; ")

  return(output)
}

# configs

# input
genes_file <- commandArgs(trailingOnly = TRUE)[1]
background_file <- commandArgs(trailingOnly = TRUE)[2]
database_path <- commandArgs(trailingOnly = TRUE)[3]
motif2tf_path <- commandArgs(trailingOnly = TRUE)[4]

# output
result_path <- commandArgs(trailingOnly = TRUE)[5]

# parameters
gene_set_name <- commandArgs(trailingOnly = TRUE)[6]
motifAnnot_highConfCat <- strsplit(commandArgs(trailingOnly = TRUE)[7], ",")[[1]]
motifAnnot_lowConfCat <- strsplit(commandArgs(trailingOnly = TRUE)[8], ",")[[1]]
nesThreshold <- as.numeric(commandArgs(trailingOnly = TRUE)[9])
aucMaxRank_factor <- as.numeric(commandArgs(trailingOnly = TRUE)[10])
geneErnMethod <- commandArgs(trailingOnly = TRUE)[11]
geneErnMaxRank <- as.numeric(commandArgs(trailingOnly = TRUE)[12])
cores_n <- as.numeric(commandArgs(trailingOnly = TRUE)[13])

# load query and background gene sets
geneSets <- list()
geneSets[[gene_set_name]] <- readLines(genes_file)

background <- readLines(background_file)

# load database, filter for background and re-rank
rankingsDb <- importRankings(database_path, columns = background)
motifRankings <- reRank(rankingsDb)
ranking_df <- getRanking(motifRankings)

# subset gene list for supported genes
geneSets[[gene_set_name]] <- intersect(colnames(ranking_df), geneSets[[gene_set_name]])

# quit early if there is no overlap
if (length(geneSets[[gene_set_name]]) == 0) {
  print("No overlap between ranking database and query genes.")
  file.create(result_path)
  quit(save = "no", status = 0)
}

# load the motif to TF annotation
motifAnnot <- importAnnotations(motif2tf_path, motifsInRanking = ranking_df$features)

###### RcisTarget

# run RcisTarget with try/catch exception handling
tryCatch({
  motifEnrichmentTable_wGenes <- cisTarget(
    geneSets = geneSets,
    motifRankings = motifRankings,
    motifAnnot = motifAnnot,
    motifAnnot_highConfCat = c(motifAnnot_highConfCat),
    motifAnnot_lowConfCat = c(motifAnnot_lowConfCat),
    highlightTFs = NULL,
    nesThreshold = nesThreshold,
    aucMaxRank = aucMaxRank_factor * ncol(motifRankings),
    geneErnMethod = geneErnMethod,
    geneErnMaxRank = geneErnMaxRank,
    nCores = cores_n,
    verbose = TRUE
  )

  # format result table
  motifEnrichmentTable_wGenes$description <- sapply(motifEnrichmentTable_wGenes$TF_highConf, process_genes)
  motifEnrichmentTable_wGenes$description <- paste0(motifEnrichmentTable_wGenes$motif, " (", motifEnrichmentTable_wGenes$description, ")")
  motifEnrichmentTable_wGenes$name <- gene_set_name

  # save result table
  fwrite(as.data.frame(motifEnrichmentTable_wGenes), file = file.path(result_path), row.names = FALSE) # quote=FALSE
}, error = function(e) {
  print("An error occurred during the cisTarget analysis.")
  overlap_percentage <- round(length(intersect(geneSets[[gene_set_name]], background)) / length(background) * 100, 2)
  print(paste("Overlap between query and background gene set might be too high with ", overlap_percentage, "%."))
  print(e)
  file.create(result_path)
  quit(save = "no", status = 0)
})
