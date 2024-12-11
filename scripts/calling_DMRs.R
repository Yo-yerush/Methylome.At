calling_DMRs <- function(methylationDataReplicates_joints, var1, var2,
                         var1_path, var2_path, comparison_name, context,
                         minProportionDiff, binSize, pValueThreshold,
                         minCytosinesCount, minReadsPerCytosine, n.cores) {
  
  condition = c(rep(var1, length(var1_path)),
                rep(var2, length(var2_path)))
  
  if (context == "CG") {
    minProportionDifference_var = minProportionDiff[1]
  } else if (context == "CHG") {
    minProportionDifference_var = minProportionDiff[2]
  } else if (context == "CHH") {
    minProportionDifference_var = minProportionDiff[3]
  }
  message(paste0("\tmin difference in methylation proportion: ",minProportionDifference_var))
  
  
  # get the range for each chromosome (needed to run in parallel)
  methData_split <- split(methylationDataReplicates_joints, seqnames(methylationDataReplicates_joints))
  start_min <- min(start(methData_split))
  end_max <- max(end(methData_split))
  chromosome_ranges = GRanges(seqnames = names(methData_split), IRanges(start = start_min, end = end_max))
  
  DMRsReplicates.0 = GRanges()
  for (i_chr in 1:length(chromosome_ranges)) {
    tryCatch({
    DMRsReplicates.loop <- computeDMRsReplicates(methylationDataReplicates_joints,
                                              condition = condition,
                                              regions = chromosome_ranges[i_chr],
                                              context = context,
                                              method = "bins",
                                              binSize = binSize,
                                              test = "betareg",
                                              pseudocountM = 1,
                                              pseudocountN = 2,
                                              pValueThreshold = 0.99,
                                              minCytosinesCount = minCytosinesCount,
                                              minProportionDifference = minProportionDifference_var,
                                              minGap = 0,
                                              minSize = 1,
                                              minReadsPerCytosine = minReadsPerCytosine,
                                              cores = n.cores)
    DMRsReplicates.0 = c(DMRsReplicates.0, DMRsReplicates.loop)
           }, error = function(cond) {
      DMRsReplicates.0 = c(DMRsReplicates.0, NULL)
      message("\t* fail to calculate DMRs in chromosome ", i_chr)
    })
  }

  DMRsReplicates = DMRsReplicates.0[DMRsReplicates.0$pValue <= pValueThreshold,]
  mcols(DMRsReplicates)[, paste0("proportionsR", 1:5)] <- NULL # remove 'NA' cols
  DMRsReplicates$log2FC = log2(DMRsReplicates$proportion2/DMRsReplicates$proportion1)
  
  DMRsReplicates = as.data.frame(DMRsReplicates) %>%
    dplyr::relocate(pValue, log2FC, context, .after = strand) %>%
    makeGRangesFromDataFrame(keep.extra.columns = T) %>%
    sort()
  
  write.csv(DMRsReplicates,
            paste0("DMRs_",context,"_",comparison_name,".csv"),
            row.names = F)
  
  return(DMRsReplicates)
}