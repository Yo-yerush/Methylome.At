calling_DMRs <- function(methylationDataReplicates_joints, meth_var1, meth_var2,
                         var1, var2, var1_path, var2_path, comparison_name,
                         context, minProportionDiff, binSize, pValueThreshold,
                         minCytosinesCount, minReadsPerCytosine, n.cores, is_Replicates) {
  condition <- c(
    rep(var1, length(var1_path)),
    rep(var2, length(var2_path))
  )

  if (context == "CG") {
    minProportionDifference_var <- minProportionDiff[1]
  } else if (context == "CHG") {
    minProportionDifference_var <- minProportionDiff[2]
  } else if (context == "CHH") {
    minProportionDifference_var <- minProportionDiff[3]
  }
  message(paste0("\tmin difference in methylation proportion: ", minProportionDifference_var))


  # get the range for each chromosome (needed to run in parallel)
  methData_split <- split(methylationDataReplicates_joints, seqnames(methylationDataReplicates_joints))
  methData_split <- methData_split[order(names(methData_split))]
  start_min <- min(start(methData_split))
  end_max <- max(end(methData_split))
  chromosome_ranges <- GRanges(seqnames = names(methData_split), IRanges(start = start_min, end = end_max))

  DMRs_gr <- GRanges()
  for (i_chr in 1:length(chromosome_ranges)) {
    tryCatch(
      {
        # both treatment are in replicates
        if (is_Replicates) {
          # Define a function for replicates
          runReplicates <- function(cores) {
            invisible(
              capture.output(
                x <- computeDMRsReplicates(
                  methylationDataReplicates_joints,
                  condition = condition,
                  regions = chromosome_ranges[i_chr],
                  context = context,
                  method = "bins",
                  binSize = binSize,
                  test = "betareg", # for replicates
                  pseudocountM = 1,
                  pseudocountN = 2,
                  pValueThreshold = 0.05,
                  minCytosinesCount = minCytosinesCount,
                  minProportionDifference = minProportionDifference_var,
                  minGap = 0,
                  minSize = 1,
                  minReadsPerCytosine = minReadsPerCytosine,
                  cores = cores
                )
              )
            )
            return(x)
          }

          DMRs_gr.loop <- tryCatch(
            {
              runReplicates(n.cores)
            },
            error = function(e) {
              if (n.cores > 10) {
                cat("Error encountered. Retrying with cores = 10\n")
                runReplicates(10)
              } else if (n.cores > 1) {
                cat("Error encountered. Retrying with cores = 1\n")
                runReplicates(1)
              } else {
                stop(e)
              }
            }
          )
        } else {
          # one or both of treatment are single sample
          # use pooled data
          runPooled <- function(cores) {
            invisible(
              capture.output(
                x <- computeDMRs(
                  meth_var1,
                  meth_var2,
                  regions = chromosome_ranges[i_chr],
                  context = context,
                  method = "bins",
                  binSize = binSize,
                  test = "fisher", # for single samples
                  pValueThreshold = 0.05,
                  minCytosinesCount = minCytosinesCount,
                  minProportionDifference = minProportionDifference_var,
                  minGap = 0,
                  minSize = 1,
                  minReadsPerCytosine = minReadsPerCytosine,
                  cores = cores
                )
              )
            )
            return(x)
          }

          DMRs_gr.loop <- tryCatch(
            {
              runPooled(n.cores)
            },
            error = function(e) {
              if (n.cores > 10) {
                cat("Error encountered. Retrying with cores = 10\n")
                runPooled(10)
              } else if (n.cores > 1) {
                cat("Error encountered. Retrying with cores = 1\n")
                runPooled(1)
              } else {
                stop(e)
              }
            }
          )
        }

        DMRs_gr <- c(DMRs_gr, DMRs_gr.loop)
      },
      error = function(cond) {
        DMRs_gr <- c(DMRs_gr, NULL)
        message("\t* fail to calculate DMRs in chromosome ", i_chr)
      }
    )
  }

  if (length(DMRs_gr) != 0) {
    mcols(DMRs_gr)[, paste0("proportionsR", 1:length(condition))] <- NULL # remove 'NA' cols
    
    # normelized log2FC (to not get INF values)
    DMRs_gr$log2FC <- log2((DMRs_gr$proportion2 + 1e-5) / (DMRs_gr$proportion1 + 1e-5))

    DMRs_gr <- as.data.frame(DMRs_gr) %>%
      dplyr::relocate(pValue, log2FC, context, .after = strand) %>%
      makeGRangesFromDataFrame(keep.extra.columns = T) %>%
      sort()
  }

  write.csv(DMRs_gr,
    paste0("DMRs_", context, "_", comparison_name, ".csv"),
    row.names = F
  )

  return(DMRs_gr)
}
