#!/usr/bin/env Rscript
# Regenerate the RcisTarget fixture kit:
#   test/fixtures/test/data/CorcesTxt/{Bcell_TF_genes,Ery_TF_genes,background_genes}.txt
#   test/fixtures/test/resources/enrichment_analysis/rcistarget_test_db.feather
#
# The upstream annotation rows ending in ".txt" are raw gene lists (one gene
# per line). The port fans the RcisTarget rules over config.txt_gene_sets;
# these fixtures are the smallest set that exercises the full path:
#   - Bcell_TF_genes / Ery_TF_genes: the two query .txt gene lists;
#   - background_genes.txt: the shared background (query genes + padding).
#
# The rankings feather is a synthetic gene-motif ranking DB in the RcisTarget
# format (rows = motifs from the real 89_motifs_test.tbl, first column
# "features", remaining columns = genes). Two motifs carry designed signal:
# bergman__Stat92E ranks the B-cell genes on top, bergman__Su_H_ the
# erythroid genes, so the default nesThreshold=3 run returns real enrichment
# rows for each query set. The remaining ~5954 padding genes are needed
# because RcisTarget's aprox gene-enrichment method requires
# maxRank + nMean <= ncol(rankings) with the port's defaults
# (geneErnMaxRank=5000, nMean=50), i.e. a gene universe > 5050 entries.
suppressMessages(library(arrow))

bcell <- c("CD19","CD79A","MS4A1","CD79B","BLNK","PAX5","IGHM","PTPRC",
           "TCF3","EBF1","IRF4","BCL6","SPIB","MYC","BTK","IL7R","LEF1","TCF7")
ery <- c("HBB","GATA1","KLF1","ALAS2","EPB42","GATA2","TAL1","NFE2",
         "KLF4","SPI1","FLI1","ERG","FOS","JUN","ATF3")
shared <- c("EGFR","MAPK1","CD14","LYZ","S100A8","RUNX1","ETS1",
            "FOXO1","ZEB2","CCND3","GABPA","ELF1","GABPB1")
genes <- c(bcell, ery, shared)
genes <- c(genes, paste0("PAD", sprintf("%04d", 1:5954)))

dir.create("test/fixtures/test/data/CorcesTxt", recursive = TRUE, showWarnings = FALSE)
for (nm in c("Bcell_TF_genes", "Ery_TF_genes")) {
  gl <- if (nm == "Bcell_TF_genes") bcell else ery
  writeLines(gl, sprintf("test/fixtures/test/data/CorcesTxt/%s.txt", nm))
}
writeLines(genes, "test/fixtures/test/data/CorcesTxt/background_genes.txt")

tbl <- read.delim("test/fixtures/test/resources/enrichment_analysis/89_motifs_test.tbl",
                  sep = "\t", header = TRUE, comment.char = "", stringsAsFactors = FALSE)
motifs <- unique(tbl[[1]])
m <- matrix(0L, nrow = length(motifs), ncol = length(genes))
for (i in seq_along(motifs)) {
  for (j in seq_along(genes)) {
    if (motifs[i] == "bergman__Stat92E" && is.element(genes[j], bcell))
      m[i, j] <- match(genes[j], bcell) - 1L
    else if (motifs[i] == "bergman__Su_H_" && is.element(genes[j], ery))
      m[i, j] <- match(genes[j], ery) - 1L
    else m[i, j] <- (7L * (i - 1L) + 13L * (j - 1L)) %% length(genes) + length(genes) %/% 2L
  }
}
df <- data.frame(features = motifs, stringsAsFactors = FALSE)
df[genes] <- as.data.frame(m)
arrow::write_feather(df, "test/fixtures/test/resources/enrichment_analysis/rcistarget_test_db.feather")
cat("wrote", nrow(df), "x", ncol(df), "feather and", length(genes), "background genes\n")
