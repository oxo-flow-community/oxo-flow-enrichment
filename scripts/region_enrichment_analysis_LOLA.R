# load libraries
library("LOLA")
library("GenomicRanges")
library("data.table")

# Port of upstream workflow/scripts/region_enrichment_analysis_LOLA.R.
# snakemake@input[["regions"|"background"|"database"]],
# snakemake@output[["result"]], snakemake@wildcards[["database"|"region_set"]]
# and snakemake@config[["genome"]] become positional CLI arguments
# <regions> <background> <database> <result> <database_name> <genome> <region_set>.

# input
query_regions <- commandArgs(trailingOnly = TRUE)[1]
background_regions <- commandArgs(trailingOnly = TRUE)[2]
database_path <- commandArgs(trailingOnly = TRUE)[3]

# output
result_path <- commandArgs(trailingOnly = TRUE)[4]

# parameters
database_name <- commandArgs(trailingOnly = TRUE)[5]
genome <- commandArgs(trailingOnly = TRUE)[6] #"hg38" "mm10"
region_set <- commandArgs(trailingOnly = TRUE)[7]

# stop early for empty query or background input
if (file.size(query_regions) == 0L || file.size(background_regions) == 0L) {
    file.create(result_path)
    quit(save = "no", status = 0)
}

### load data

# load query region sets
regionSet_query <- readBed(query_regions)

# load background/universe region sets (e.g., consensus region set)
regionSet_background <- readBed(background_regions)

# requires resources downloaded from: https://databio.org/regiondb
# requires simpleCache package installed
database <- loadRegionDB(file.path(database_path))

###### LOLA

# run LOLA
res <- runLOLA(regionSet_query, regionSet_background, database, cores=1)

# make description more descriptive
if (database_name=='LOLACore'){
    res$description <- paste(res$description, res$cellType, res$antibody, sep='.')
}else{
    res$description <- paste(res$description, res$filename, sep='.')
}

# format description that values are unique
res$description <- make.names(res$description, unique=TRUE)

# determine raw p-value
res$pValue <- ('^'(10,-1*res[['pValueLog']]))

# enrichment columns consumed downstream (upstream enrichment table contract)
res$oddsRatio <- res$b / res$c
res$qValue <- p.adjust(res$pValue, method = "BH")

# save results
fwrite(as.data.frame(res), file=file.path(result_path), row.names=FALSE)
