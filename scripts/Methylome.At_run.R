suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(DMRcaller))
suppressMessages(library(rtracklayer))
suppressMessages(library(lattice))
suppressMessages(library(PeakSegDisk))
suppressMessages(library(topGO))
suppressMessages(library(KEGGREST))
suppressMessages(library(Rgraphviz))
suppressMessages(library(org.At.tair.db))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(geomtextpath))
suppressMessages(library(plyranges))
suppressMessages(library(parallel))
suppressMessages(library("RColorBrewer"))
suppressMessages(library(circlize))

# configuration from the command line arguments (in Bash)
configs <- commandArgs(trailingOnly = TRUE)

# upload samples file
var_sep_format = ifelse(grepl("\\.csv$", configs[1]), ",", "\t")
var_table = suppressWarnings(read.csv(configs[1], header = F, sep = var_sep_format))
vars_vector = unique(var_table[,1])
var1_path = var_table[grep(vars_vector[1], var_table[,1]), 2]
var2_path = var_table[grep(vars_vector[2], var_table[,1]), 2]

# remove "/" from the end of a path strings
rmv_d <- function(x) {
  if (substr(x, nchar(x), nchar(x)) == "/") {
    x = substr(x, 1, nchar(x)-1)
  }
  return(x)
}

# run the main function
source(paste0(rmv_d(configs[2]),"/scripts/Methylome.At_main.R"))
try({
  Methylome.At_main(var1 = vars_vector[1],
                    var2 = vars_vector[2],
                    var1_path = var1_path,
                    var2_path = var2_path,
                    Methylome.At_path = rmv_d(configs[2]),
                    annotation_file = configs[3],
                    description_file = configs[4],
                    TEs_file = configs[5],
                    minProportionDiff = as.numeric(configs[6:8]), # CG, CHG, CHH
                    binSize = as.numeric(configs[9]),
                    minCytosinesCount = as.numeric(configs[10]),
                    minReadsPerCytosine = as.numeric(configs[11]),
                    pValueThreshold = as.numeric(configs[12]),
                    n.cores = as.numeric(configs[13]),
                    GO_analysis = configs[14],
                    KEGG_pathways = configs[15]
  )
  #message("\n\nwarnings:\n", warnings())
})