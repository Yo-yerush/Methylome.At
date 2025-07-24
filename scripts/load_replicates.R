load_replicates <- function(var_path, n.cores, var_name, metaPlot_script=F, file_type="CX_report") {
  n.cores.f <- ifelse(n.cores > 2, length(var_path), 1)
  ### load the data for replicates
  methylationDataList <- mclapply(var_path,
    function(cx_4_read_0) {
      if (file_type != "CX_report") {
        dir.create("CX_report_files", showWarnings = F)
        prefix_4_cx <- sub(paste0("\\.", file_type, "$"), "", basename(cx_4_read_0))

        # convert to 'CX_report' if the files type is 'bedMeyhl' or 'CGmap'
        if (file_type == "bed") {
          cat("convert 'bedMethyl' to 'CX_report'...\n")
          system(paste0("scripts/bedmethyl_2_cx.sh -i ", cx_4_read_0, " -o CX_report_files/", prefix_4_cx))
        } else if (file_type == "CGmap") {
          cat("convert 'CGmap' to 'CX_report'...\n")
          system(paste0("scripts/cgmap_2_cx.sh -i ", cx_4_read_0, " -o CX_report_files/", prefix_4_cx))
        }
        cx_4_read <- paste0("CX_report_files/", prefix_4_cx, ".CX_report.txt")
        cat("done:", cx_4_read, "\n")
      } else {
        cx_4_read <- cx_4_read_0
      }

      msg_0 <- capture.output(rb_out <- readBismark(cx_4_read), type = "message")
      if (length(msg_0) > 0 && length(strsplit(msg_0, " ")[[1]]) >= 2) {
        msg <- strsplit(msg_0, " ")[[1]][2]
      } else {
        msg <- "unknown"
      }
      message(paste0("\tread ", msg, " C's [", var_name, "]"))
      
      rb_out
    },
    mc.cores = n.cores.f
  )

  # Pooling Methylation Datasets
  methylationData_pool <- poolMethylationDatasets(GRangesList(methylationDataList))

  # Joining Replicates
  if (!metaPlot_script) {
    if (length(var_path) == 1) {
      methylationDataReplicates <- methylationDataList[[1]] # 1 sample
      names(mcols(methylationDataReplicates))[names(mcols(methylationDataReplicates)) %in% c("readsM", "readsN")] <- c("readsM1", "readsN1")
    } else {
      message(paste0("\t----"))
      message(paste0("\tjoin replicates   [", var_name, "]"))
      methylationDataReplicates <- joinReplicates(methylationDataList[[1]], methylationDataList[[2]]) # 2 samples
      if (length(var_path) == 3) {
        methylationDataReplicates <- joinReplicates(methylationDataReplicates, methylationDataList[[3]]) # 3 samples
      }
    }

    return(list(
      methylationData_pool = methylationData_pool,
      methylationDataReplicates = methylationDataReplicates
    ))
  } else {
    return(methylationData_pool)
  }
}
