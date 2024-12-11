load_replicates <- function(var_path, n.cores, metPlot_script=F) {
  
  n.cores.f = ifelse(n.cores > 2, length(var_path), 1)
  ### load the data for replicates
  methylationDataList <- mclapply(var_path, readBismark, mc.cores = n.cores.f)
  
  # Pooling Methylation Datasets
  methylationData_pool <- poolMethylationDatasets(GRangesList(methylationDataList))
  
  # Joining Replicates
  if (!metPlot_script) {
    methylationDataReplicates = joinReplicates(methylationDataList[[1]],methylationDataList[[2]])
    if (length(var_path) == 3) {
      methylationDataReplicates = joinReplicates(methylationDataReplicates,methylationDataList[[3]])
    }
    
    return(list(methylationData_pool = methylationData_pool,
                methylationDataReplicates = methylationDataReplicates))
  } else {
    
    return(methylationData_pool)
  }
}