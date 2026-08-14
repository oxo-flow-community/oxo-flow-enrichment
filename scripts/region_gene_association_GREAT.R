# load libraries
# library("LOLA")
library("GenomicRanges")
library("rGREAT")
library("data.table")
library("rtracklayer")

# Port of upstream workflow/scripts/region_gene_association_GREAT.R.
# snakemake@input[["regions"|"database"]], snakemake@output[["genes"|
# "associations_table"|"associations_plot"]], snakemake@config[["genome"]],
# snakemake@config[["great_parameters"]] and snakemake@threads become
# positional CLI arguments:
# <regions> <database> <genes> <associations_table> <associations_plot>
# <genome> <cores> <min_gene_set_size> <mode> <basal_upstream>
# <basal_downstream> <extension>

#input
regions_file <- commandArgs(trailingOnly = TRUE)[1]
database_path <- commandArgs(trailingOnly = TRUE)[2]

# output
gene_path <- commandArgs(trailingOnly = TRUE)[3]
associations_table_path <- commandArgs(trailingOnly = TRUE)[4]
associations_plot_path <- commandArgs(trailingOnly = TRUE)[5]

# parameters
genome <- commandArgs(trailingOnly = TRUE)[6]
cores_n <- as.integer(commandArgs(trailingOnly = TRUE)[7])

# set genome
if (genome=="hg19" | genome=="hg38"){
    orgdb <- "org.Hs.eg.db"
}else if(genome=="mm9" | genome=="mm10"){
    orgdb <- "org.Mm.eg.db"
}

# stop early for empty region input
if (file.size(regions_file) == 0L){
    file.create(gene_path)
    file.create(associations_table_path)
    file.create(associations_plot_path)
    quit(save = "no", status = 0)
}

# load query region set
regionSet_query <- import(regions_file, format = "BED")

# load database
database = read_gmt(file.path(database_path), from = "SYMBOL", to = "ENTREZ", orgdb = orgdb)

###### GREAT

# run GREAT
res <- great(gr = regionSet_query,
      gene_sets = database,
      tss_source = genome,
      biomart_dataset = NULL,
      min_gene_set_size = as.integer(commandArgs(trailingOnly = TRUE)[8]), #default: 5
      mode = commandArgs(trailingOnly = TRUE)[9],
      basal_upstream = as.numeric(commandArgs(trailingOnly = TRUE)[10]),
      basal_downstream = as.numeric(commandArgs(trailingOnly = TRUE)[11]),
      extension = as.numeric(commandArgs(trailingOnly = TRUE)[12]),
      extended_tss = NULL,
      background = NULL,
      exclude = "gap",
      cores = cores_n, #default: 1
      verbose = TRUE #default: great_opt$verbose
     )

# plot gene-region association
pdf(file=file.path(associations_plot_path), width=12, height=4)
plotRegionGeneAssociations(res)
dev.off()

# get and save gene-region association
associations <- getRegionGeneAssociations(res)
associations_df <- as.data.frame(associations)
# BED uses 0-based starts; export start in BED-style coordinates for direct comparability to input BED.
associations_df$start <- associations_df$start - 1
fwrite(associations_df, file=file.path(associations_table_path), row.names=TRUE)

# save unique associated genes by using mcols(), which returns a DataFrame object containing the metadata columns.
genes <- unique(unlist(mcols(associations)$annotated_genes))
write(genes, file.path(gene_path))
