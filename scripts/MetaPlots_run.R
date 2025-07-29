start_time <- Sys.time()

# upload libraries
lib_packages <- c(
  "dplyr", "ggplot2", "DMRcaller", "org.At.tair.db",
  "GenomicFeatures", "plyranges", "parallel"
  )
for (n.pkg in seq(lib_packages)) {
    suppressWarnings(suppressMessages(library(lib_packages[n.pkg], character.only = TRUE)))
    perc_val <- (n.pkg / length(lib_packages)) * 100
    cat(paste0("\rloading libraries [", round(perc_val, 1), "%] "))
}
cat("\n")

# remove "/" from the end of a path strings
rmv_d <- function(x) {
  if (substr(x, nchar(x), nchar(x)) == "/") {
    x = substr(x, 1, nchar(x)-1)
  }
  return(x)
}

########################################################################### 

# configuration from the command line arguments (in Bash)
configs <- commandArgs(trailingOnly = TRUE)

# upload samples file
var_sep_format = ifelse(grepl("\\.csv$", configs[1]), ",", "\t")
var_table = suppressWarnings(read.csv(configs[1], header = F, sep = var_sep_format))
vars_vector = unique(var_table[,1])
var1_path = var_table[grep(vars_vector[1], var_table[,1]), 2]
var2_path = var_table[grep(vars_vector[2], var_table[,1]), 2]

# config variables
var1 = vars_vector[1] # control
var2 = vars_vector[2] # treatment
Methylome.At_path = configs[2]
annotation_file = configs[3]
TEs_file = configs[4]
minReadsPerCytosine = as.numeric(configs[5])
metaPlot.random.genes = configs[6]
n.cores = as.numeric(configs[7])
binSize = as.numeric(configs[8])
methyl_files_type = configs[9]
analyze_Gene_n_TEs = as.logical(configs[10])
analyze_GeneFeatures = as.logical(configs[11])

########################################################################### 

# scripts directory path
scripts_dir = paste0(rmv_d(Methylome.At_path),"/scripts/")

# metaPlot results directory path
dir.create(paste0(rmv_d(Methylome.At_path),"/results/",var2,"_vs_",var1), showWarnings = F)

metaPlot_path = paste0(rmv_d(Methylome.At_path),"/results/",var2,"_vs_",var1,"/metaPlots/")
dir.create(metaPlot_path, showWarnings = F)

source(paste0(scripts_dir,"trimm_and_rename_seq.R"))

########################################################################### 

##### Read annotation files #####
cat("\rload annotations files [0/2]")
# annotation file
tryCatch({
  # if its 'csv' file
  if (grepl("\\.csv$|\\.csv\\.gz$", annotation_file)) {
    annotation.gr = read.csv(annotation_file) %>% 
      makeGRangesFromDataFrame(., keep.extra.columns = T) %>%
      trimm_and_rename()
    # if its 'gtf'/'gff'/'gff3' file
  } else if (grepl("\\.gtf$|\\.gff$|\\.gff3$|\\.gtf\\.gz$|\\.gff\\.gz$|\\.gff3\\.gz$", tolower(annotation_file))) { 
    annotation.gr = import.gff3(annotation_file) %>%
      trimm_and_rename()
  }
  message("load annotation file")
}, error = function(cond) {
  cat("\n*\n",as.character(cond),"\n*\n")
  message("load 'annotation' file: fail")
})
cat("\rload annotations files [1/2]")

 # TAIR10 Transposable Elements file
tryCatch({
  source(paste0(scripts_dir,"edit_TE_file.R"))
  TE.gr = read.csv(TEs_file, sep = "\t")
  TE.gr = edit_TE_file(TE.gr)
  message("load Transposable Elements file")
}, error = function(cond) {
  cat("\n*\n",as.character(cond),"\n*\n")
  message("load Transposable Elements file: fail")
})
cat("\rload annotations files [2/2]")

cat("\n\n")

########################################################################### 

var_args = list(
  list(path = var1_path, name = var1),
  list(path = var2_path, name = var2)
)

##### load the data for replicates ##### 
message("load CX methylation data...")
source(paste0(scripts_dir,"load_replicates.R"))
tryCatch({
  # load 'CX_reports'
  n.cores.load = ifelse(n.cores > 1, 2, 1)
  load_vars = mclapply(var_args, function(x) load_replicates(x$path, n.cores, x$name, T, methyl_files_type), mc.cores = n.cores.load)
  
  # trimm seqs objects (rename if not 'TAIR10' Chr seqnames)
  meth_var1 = trimm_and_rename(load_vars[[1]])
  meth_var2 = trimm_and_rename(load_vars[[2]])
  message("load and pool replicates data: successfully\n") 
},
error = function(cond) {stop("load and pool replicates data: fail\n")}
)

cat("\n")

###########################################################################

# run metPlot function for coding-Genes and TEs
if (analyze_Gene_n_TEs) {
  ##### calculate metaPlot for Genes and TEs
  source(paste0(scripts_dir,"Genes_metaPlot_fun.R"))
  
  tryCatch({
    message("generate metaPlot to protein-coding Genes...")
    setwd(metaPlot_path)
    Genes_metaPlot(meth_var1,meth_var2,var1,var2,annotation.gr,metaPlot.random.genes,minReadsPerCytosine,n.cores)
  }, error = function(cond) {message(paste0("process average metaPlot to ",metaPlot.random.genes," Protein Coding Genes: fail"))})
  
  tryCatch({
    message("\ngenerate metaPlot to Transposable Elements...")
    setwd(metaPlot_path)
    Genes_metaPlot(meth_var1,meth_var2,var1,var2,TE.gr,metaPlot.random.genes,minReadsPerCytosine,n.cores,is_TE=T)
  }, error = function(cond) {message(paste0("process average metaPlot to ",metaPlot.random.genes," Transposable Elements: fail\n"))})
}

###########################################################################

# run metPlot function for coding-Gene features
if (analyze_GeneFeatures) {
  source(paste0(scripts_dir,"Gene_features_metaPlot_fun.R"))
  
  tryCatch({
    message("\ngenerate metaPlot to protein-coding Gene Features...")
    setwd(metaPlot_path)
    Genes_features_metaPlot(meth_var1,meth_var2,var1,var2,annotation.gr,metaPlot.random.genes,minReadsPerCytosine,binSize,n.cores)
  }, error = function(cond) {message(paste0("process average metaPlot to ",metaPlot.random.genes," Protein Coding Genes: fail"))})
}

message(paste0("**\t",var2," vs ",var1,": done\n"))

###########################################################################

# date and time of the end
message(paste0("\n**\t",
               paste(format(Sys.time(),"%d"),format(Sys.time(),"%m"),format(Sys.time(),"%Y"), sep = "-"),
               " ",format(Sys.time(),"%H:%M")))
end_time <- Sys.time()
message(paste0("**\ttime: ", round(difftime(end_time,start_time, units = "hours"),2)," hours\n"))
