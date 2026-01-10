library(DMRcaller)
library(rtracklayer)
library(dplyr)
library(parallel)
library(RColorBrewer)
library(data.table)


source("https://raw.githubusercontent.com/Yo-yerush/Methylome.At/main/scripts/trimm_and_rename_seq.R")
###########################################################################


##### Read annotation and description files #####

annotation.gr <- read.csv("C:/Users/YonatanY/Migal/Rachel Amir Team - General/Arabidopsis_db/RA_costume_annotations_files/Methylome.At_annotations.csv.gz") %>%
  makeGRangesFromDataFrame(., keep.extra.columns = T) %>%
  trimm_and_rename()

# TAIR10 Transposable Elements file

source("https://raw.githubusercontent.com/Yo-yerush/Methylome.At/main/scripts/edit_TE_file.R")
TE_file.df <- read.csv("https://raw.githubusercontent.com/Yo-yerush/Methylome.At/main/annotation_files/TAIR10_Transposable_Elements.txt", sep = "\t")
TE_file <- edit_TE_file(TE_file.df)


# upload description file

description_df <- read.csv("/home/yoyerush/yo/methylome_pipeline/leaf_senescence_vatov_22/Methylome.At/annotation_files/Methylome.At_description_file.csv.gz")
names(description_df)[1] <- "gene_id"


###########################################################################
var1 <- "wt"
var2 <- "suvh8"
var1_path <- c(
  "/home/yoyerush/yo/methylome_pipeline/other_mutants/stroud_et_al_2013/bismark_results/wt_2_bismark_se.CX_report.txt.gz",
  "/home/yoyerush/yo/methylome_pipeline/other_mutants/stroud_et_al_2013/bismark_results/wt_3_bismark_se.CX_report.txt.gz"
)
var2_path <- "/home/yoyerush/yo/methylome_pipeline/other_mutants/stroud_et_al_2013/bismark_results/suvh8_bismark_se.CX_report.txt.gz"


is_single <- (length(var1_path) == 1 & length(var2_path) == 1) # both genotypes includes 1 sample
is_Replicates <- (length(var1_path) > 1 & length(var2_path) > 1) # both genotypes includes >1 samples

var_args <- list(
  list(path = var1_path, name = var1),
  list(path = var2_path, name = var2)
)

##### load methylation data ('CX_report' file) #####
message("load CX methylation data...")
source("https://raw.githubusercontent.com/Yo-yerush/Methylome.At/main/scripts/load_replicates.R")
n.cores <- 10
methyl_files_type <- "CX_report"

# load 'CX_reports'
n.cores.load <- ifelse(n.cores > 1, 2, 1)
load_vars <- mclapply(var_args, function(x) load_replicates(x$path, n.cores, x$name, F, methyl_files_type), mc.cores = n.cores.load)

# trimm seqs objects (rename if not 'TAIR10' Chr seqnames)
meth_var1 <- trimm_and_rename(load_vars[[1]]$methylationData_pool)
meth_var2 <- trimm_and_rename(load_vars[[2]]$methylationData_pool)
meth_var1_replicates <- trimm_and_rename(load_vars[[1]]$methylationDataReplicates)
meth_var2_replicates <- trimm_and_rename(load_vars[[2]]$methylationDataReplicates)

cat("\nPooled ", var1, " object for example:\n\n", sep = "")
capture.output(print(meth_var1), file = NULL)[-c(1, 16)] %>% cat(sep = "\n")
cat(paste("seq-levels:", paste(seqlevels(meth_var1), collapse = " ")), "  -------", sep = "\n")

print(meth_var1)
