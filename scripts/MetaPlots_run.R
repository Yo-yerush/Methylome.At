suppressMessages(library(DMRcaller))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(org.At.tair.db))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(plyranges))
suppressMessages(library(parallel))

start_time <- Sys.time()

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
analyze_Gene_n_TEs = as.logical(configs[9])
analyze_GeneFeatures = as.logical(configs[10])

########################################################################### 

# scripts directory path
scripts_dir = paste0(rmv_d(Methylome.At_path),"/scripts/")

# metaPlot results directory path
dir.create(paste0(rmv_d(Methylome.At_path),"/results/"), showWarnings = F)
dir.create(paste0(rmv_d(Methylome.At_path),"/results/",var2,"_vs_",var1), showWarnings = F)

metaPlot_path = paste0(rmv_d(Methylome.At_path),"/results/",var2,"_vs_",var1,"/metaPlots/")
dir.create(metaPlot_path, showWarnings = F)

########################################################################### 

##### load the data for replicates ##### 
message("load replicates data...")
source(paste0(scripts_dir,"load_replicates.R"))
tryCatch({
  # load 'CX_reports'
  n.cores.load = ifelse(n.cores > 1, 2, 1)
  load_vars = mclapply(list(var1_path,var2_path), function(x) load_replicates(x,n.cores,T), mc.cores = n.cores.load)
  
  # trimm seqs objects (rename if not 'TAIR10' Chr seqnames)
  source(paste0(scripts_dir,"trimm_and_rename_seq.R"))
  meth_var1 = trimm_and_rename(load_vars[[1]])
  meth_var2 = trimm_and_rename(load_vars[[2]])
  message("load and pool replicates data: successfully\n") 
},
error = function(cond) {stop("load and pool replicates data: fail\n")}
)

###########################################################################

##### import annotation files ##### 
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
  message("import annotation file")
}, error = function(cond) {message("import 'annotation' file: fail")})

# TAIR10 Transposable Elements file
tryCatch({
  source(paste0(scripts_dir,"edit_TE_file.R"))
  TE.gr = read.csv(TEs_file, sep = "\t")
  TE.gr = edit_TE_file(TE.gr)
  message("import Transposable Elements file\n")
}, error = function(cond) {message("import Transposable Elements file: fail\n")})

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
    message("generate metaPlot to Transposable Elements...")
    setwd(metaPlot_path)
    Genes_metaPlot(meth_var1,meth_var2,var1,var2,TE.gr,metaPlot.random.genes,minReadsPerCytosine,n.cores,is_TE=T)
  }, error = function(cond) {message(paste0("process average metaPlot to ",metaPlot.random.genes," Transposable Elements: fail\n"))})
}

###########################################################################

# run metPlot function for coding-Gene features
if (analyze_GeneFeatures) {
  source(paste0(scripts_dir,"Gene_features_metaPlot_fun.R"))
  
  tryCatch({
    message("generate metaPlot to protein-coding Gene Features...")
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