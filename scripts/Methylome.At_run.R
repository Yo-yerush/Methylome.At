options(width = 100, warn = -1)

# upload libraries
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
cat("\n")

# configuration from the command line arguments (in Bash)
configs <- commandArgs(trailingOnly = TRUE)

# upload samples file
suppressWarnings({
  var_sep_format <- ifelse(
    grepl("\\.csv$", configs[1]) | grepl(",", readLines(configs[1])[1]),
    ",", "\t"
  )
  var_table <- read.csv(configs[1], header = F, sep = var_sep_format)
})
vars_vector <- unique(var_table[, 1])
var1_path <- var_table[grep(vars_vector[1], var_table[, 1]), 2]
var2_path <- var_table[grep(vars_vector[2], var_table[, 1]), 2]

# remove "/" from the end of a path strings
rmv_d <- function(x) {
  if (substr(x, nchar(x), nchar(x)) == "/") {
    x <- substr(x, 1, nchar(x) - 1)
  }
  return(x)
}

# run the main function
source(paste0(rmv_d(configs[2]), "/scripts/Methylome.At_main.R"))
try({
  Methylome.At_main(
    var1 = vars_vector[1],
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
    methyl_files_type = configs[13],
    img_type = configs[14],
    n.cores = as.numeric(configs[15]),
    GO_analysis = configs[16],
    KEGG_pathways = configs[17],
    analyze_dH = configs[18],
    TE_metaPlots = configs[19],
    GeneBody_metaPlots = configs[20],
    GeneFeatures_metaPlots = configs[21],
    gene_features_binSize = configs[22],
    metaPlot.random.genes = configs[23]
  )
  # message("\n\nwarnings:\n", warnings())
})
