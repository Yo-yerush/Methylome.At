# configurations
var1 <- "wt"
var2 <- "suvh8"
var1_path <- c(
  "/home/yoyerush/yo/methylome_pipeline/other_mutants/stroud_et_al_2013/bismark_results/wt_2_bismark_se.CX_report.txt.gz",
  "/home/yoyerush/yo/methylome_pipeline/other_mutants/stroud_et_al_2013/bismark_results/wt_3_bismark_se.CX_report.txt.gz"
)
var2_path <- "/home/yoyerush/yo/methylome_pipeline/other_mutants/stroud_et_al_2013/bismark_results/suvh8_bismark_se.CX_report.txt.gz"
Methylome.At_path <- "/home/yoyerush/yo/methylome_pipeline/Methylome.At_020126/Methylome.At/"
annotation_file <- "/home/yoyerush/yo/methylome_pipeline/Methylome.At_020126/Methylome.At/annotation_files/Methylome.At_annotations.csv.gz"
description_file <- "/home/yoyerush/yo/methylome_pipeline/Methylome.At_020126/Methylome.At/annotation_files/Methylome.At_description_file.csv.gz"
TEs_file <- "/home/yoyerush/yo/methylome_pipeline/Methylome.At_020126/Methylome.At/annotation_files/TAIR10_Transposable_Elements.txt"
minProportionDiff <- c(0.4, 0.2, 0.1)
binSize <- 100
minCytosinesCount <- 4
minReadsPerCytosine <- 4
pValueThreshold <- 0.05
methyl_files_type <- "CX_report"
img_type <- "svg"
n.cores <- 30
GO_analysis <- FALSE
KEGG_pathways <- FALSE
analyze_dH <- TRUE
TE_metaPlots <- TRUE
GeneBody_metaPlots <- TRUE
GeneFeatures_metaPlots <- TRUE
gene_features_binSize <- 10
metaPlot.random.genes <- 5000

### load libreries
lib_packages <- c(
  "dplyr", "tidyr", "ggplot2", "DMRcaller", "rtracklayer", "lattice",
  "PeakSegDisk", "topGO", "KEGGREST", "Rgraphviz", "org.At.tair.db",
  "GenomicFeatures", "geomtextpath", "plyranges", "parallel",
  "RColorBrewer", "circlize", "cowplot", "knitr", "data.table"
)
for (n.pkg in seq(lib_packages)) {
  tryCatch(
    {
      suppressWarnings(suppressMessages(library(lib_packages[n.pkg], character.only = TRUE)))
      perc_val <- (n.pkg / length(lib_packages)) * 100
      cat(paste0("\rloading libraries [", round(perc_val, 1), "%] "))
    },
    error = function(e) {
      cat(paste0("\nError loading ", lib_packages[n.pkg], ": ", e$message, "\n"))
      message(paste0("\n* Error loading ", lib_packages[n.pkg], "package\n"))
    }
  )
}