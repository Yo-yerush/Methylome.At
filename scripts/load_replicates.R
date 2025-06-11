load_replicates <- function(var_path, n.cores, var_name, metaPlot_script=F) {
  
  n.cores.f = ifelse(n.cores > 2, length(var_path), 1)
  ### load the data for replicates
  methylationDataList <- mclapply(var_path,
                                  function(p) {
                                    msg_0 <- capture.output(rb_out <- readBismark(p), type = "message")
                                    msg <- strsplit(msg_0, " ")[[1]][2]
                                    message(paste0("\tread ", msg, " C's [",var_name,"]"))
                                    rb_out
                                    },
                                  mc.cores = n.cores.f)
  
  # Pooling Methylation Datasets
  methylationData_pool <- poolMethylationDatasets(GRangesList(methylationDataList))
  
  # Joining Replicates
  if (!metaPlot_script) {

    if (length(var_path) == 1) {
      methylationDataReplicates = methylationDataList[[1]] # 1 sample
      names(mcols(methylationDataReplicates))[names(mcols(methylationDataReplicates)) %in% c("readsM", "readsN")] = c("readsM1", "readsN1")
    } else {
      message(paste0("\t----"))
      message(paste0("\tjoin replicates   [", var_name, "]"))
      methylationDataReplicates = joinReplicates(methylationDataList[[1]], methylationDataList[[2]]) # 2 samples
      if (length(var_path) == 3) {
        methylationDataReplicates = joinReplicates(methylationDataReplicates, methylationDataList[[3]]) # 3 samples
      }
    }
    
    return(list(methylationData_pool = methylationData_pool,
                methylationDataReplicates = methylationDataReplicates))
  } else {
    
    return(methylationData_pool)
  }
}