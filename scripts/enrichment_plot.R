# load libraries
library("ggplot2")
library("svglite")
library("data.table")

# Port of upstream workflow/scripts/enrichment_plot.R.
# snakemake@input[["enrichment_result"]], snakemake@output[["enrichment_plot"]],
# snakemake@wildcards[["tool"|"db"|"feature_set"]], snakemake@params[["utils_path"]]
# and the config["column_names"] entries of the tool become positional CLI
# arguments:
# <enrichment_result> <enrichment_plot> <tool> <db> <feature_set> <utils_path>
# <top_n> <p_value> <adj_pvalue> <effect_size> <overlap> <term>

# source utility functions
# source("workflow/scripts/utils.R")
# snakemake@source("./utils.R") # does not work when loaded as module (https://github.com/snakemake/snakemake/issues/2205)
source(commandArgs(trailingOnly = TRUE)[6])

# draw enrichment plot
do_enrichment_plot <- function(plot_data, title, x, y, size, colorBy, font.size, path, filename, top_n){
    plot_height <- max(80, 10 * top_n)

    enr_p <- ggplot(plot_data, aes_string(x=x, y=y, size=size, color=colorBy))  +
        geom_point() +
        scale_color_continuous(low="red", high="blue", name = colorBy, guide=guide_colorbar(reverse=TRUE)) +
        ggtitle(title) +
#         theme_dose(font.size) +
        scale_size(range=c(3, 8)) +
        scale_y_discrete(label=addline_format, limits=rev) +
        theme(axis.text.y=element_text(vjust=0.6))

    ggsave_new(filename = filename,
           results_path=path,
           plot=enr_p,
           width=200,
           height=plot_height,
              units = "mm")
}

# configs

# input
enrichment_result_path <- commandArgs(trailingOnly = TRUE)[1]

# output
enrichment_plot_path <- commandArgs(trailingOnly = TRUE)[2]

# parameters
tool <- commandArgs(trailingOnly = TRUE)[3]
database <- commandArgs(trailingOnly = TRUE)[4]
feature_set <- commandArgs(trailingOnly = TRUE)[5]
plot_context <- paste(tool, database, feature_set, sep = " | ")

top_n <- as.numeric(commandArgs(trailingOnly = TRUE)[7])
pval_col <- commandArgs(trailingOnly = TRUE)[8]
adjp_col <- commandArgs(trailingOnly = TRUE)[9]
effect_col <- commandArgs(trailingOnly = TRUE)[10]
overlap_col <- commandArgs(trailingOnly = TRUE)[11]
term_col <- commandArgs(trailingOnly = TRUE)[12]

# load enrichment result
if (file.size(enrichment_result_path) != 0L){
    enrichment_result <- data.frame(fread(file.path(enrichment_result_path), header=TRUE))
}else{
    make_message_plot(enrichment_plot_path, "No results found\nenrichment result file is empty", plot_context)
    quit(save = "no", status = 0)
}

# stop early for header-only or empty result tables
if (nrow(enrichment_result) == 0L){
    make_message_plot(enrichment_plot_path, "No results found\nenrichment result has no rows", plot_context)
    quit(save = "no", status = 0)
}

top_n <- min(top_n, nrow(enrichment_result))

# evaluate overlap numerically if necessary
if(class(enrichment_result[[overlap_col]])=="character"){
    enrichment_result[[overlap_col]] <- as.numeric(lapply(enrichment_result[[overlap_col]], evaltext))
}

# calculate comparable effect size either NES or odds-ratio/fold based
if (tool!="preranked_GSEApy" & tool!="pycisTarget" & tool!="RcisTarget"){
    # calculate log2(effect-size) and put in new column
    effect_col_new <- paste0("log2_",effect_col)
    enrichment_result[[effect_col_new]] <- log2(enrichment_result[[effect_col]])
    effect_col <- effect_col_new
}

# determine ranks
enrichment_result$PValue_Rnk <- if (tool!="pycisTarget" & tool!="RcisTarget") rank(enrichment_result[[pval_col]]) else rank(-enrichment_result[[pval_col]])
enrichment_result$Fold_Rnk <- rank(-abs(enrichment_result[[effect_col]]))
enrichment_result$Coverage_Rnk <- rank(-enrichment_result[[overlap_col]])
# calculate and sort by mean rank
enrichment_result$meanRnk <- rowMeans(enrichment_result[,c('PValue_Rnk', 'Fold_Rnk','Coverage_Rnk')])
enrichment_result <- enrichment_result[order(enrichment_result$meanRnk, decreasing=FALSE),]

# format term column that order is kept and values are unique
enrichment_result[[term_col]] <- make.unique(as.character(enrichment_result[[term_col]]), sep = "_") #make.names(enrichment_result[[term_col]], unique=TRUE)
enrichment_result[[term_col]] <- factor(enrichment_result[[term_col]], levels = enrichment_result[[term_col]])

# plot top_n terms by mean_rnk
do_enrichment_plot(plot_data=enrichment_result[1:top_n,],
               title=paste0(tool, ' results of ',feature_set,'\nin ',database),
               x=effect_col,
               y=term_col,
               size=overlap_col,
               colorBy=adjp_col,
               font.size=10,
               path=file.path(dirname(enrichment_plot_path)),
               filename=paste0(feature_set,"_",database),
                   top_n = top_n
              )
