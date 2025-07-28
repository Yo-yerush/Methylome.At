Methylome.At_main <- function(var1, # control
                              var2, # treatment
                              var1_path,
                              var2_path,
                              Methylome.At_path,
                              annotation_file,
                              description_file,
                              TEs_file,
                              minProportionDiff, # CG, CHG, CHH
                              binSize,
                              minCytosinesCount,
                              minReadsPerCytosine,
                              pValueThreshold,
                              methyl_files_type,
                              n.cores,
                              GO_analysis,
                              KEGG_pathways) {
  ###########################################################################

  start_time <- Sys.time()
  scripts_dir <- paste0(Methylome.At_path, "/scripts/")
  source(paste0(scripts_dir, "trimm_and_rename_seq.R"))

  ###########################################################################

  setwd(Methylome.At_path)

  ##### Read annotation and description files #####
  cat("\rload annotations and description files [0/3]")

  # annotation file
  tryCatch(
    {
      # if its 'csv' file
      if (grepl("\\.csv$|\\.csv\\.gz$", annotation_file)) {
        annotation.gr <- read.csv(annotation_file) %>%
          makeGRangesFromDataFrame(., keep.extra.columns = T) %>%
          trimm_and_rename()
        # if its 'gtf'/'gff'/'gff3' file
      } else if (grepl("\\.gtf$|\\.gff$|\\.gff3$|\\.gtf\\.gz$|\\.gff\\.gz$|\\.gff3\\.gz$", tolower(annotation_file))) {
        annotation.gr <- import.gff3(annotation_file) %>%
          trimm_and_rename()
      }
      message("load annotation file")
    },
    error = function(cond) {
      cat("\n*\n", as.character(cond), "\n*\n")
      message("load 'annotation' file: fail")
    }
  )
  cat("\rload annotations and description files [1/3]")

  # TAIR10 Transposable Elements file
  tryCatch(
    {
      source(paste0(scripts_dir, "edit_TE_file.R"))
      TE_file.df <- read.csv(TEs_file, sep = "\t")
      TE_file <- edit_TE_file(TE_file.df)
      message("load Transposable Elements file")
    },
    error = function(cond) {
      cat("\n*\n", as.character(cond), "\n*\n")
      message("load Transposable Elements file: fail")
    }
  )
  cat("\rload annotations and description files [2/3]")

  # upload description file
  tryCatch(
    {
      des_file_sep <- ifelse(grepl("\\.csv$|\\.csv\\.gz$", description_file), ",", "\t")
      description_df <- read.csv(description_file, sep = des_file_sep)
      names(description_df)[1] <- "gene_id"
      message("load description file\n")
    },
    error = function(cond) {
      cat("\n*\n", as.character(cond), "\n*\n")
      message("load description file: fail\n")
    }
  )
  cat("\rload annotations and description files [3/3]")

  cat("\n\n")

  ###########################################################################
  setwd(Methylome.At_path)

  is_single <- (length(var1_path) == 1 & length(var2_path) == 1) # both genotypes includes 1 sample
  is_Replicates <- (length(var1_path) > 1 & length(var2_path) > 1) # both genotypes includes >1 samples

  var_args <- list(
    list(path = var1_path, name = var1),
    list(path = var2_path, name = var2)
  )

  ##### load methylation data ('CX_report' file) #####
  message("load CX methylation data...")
  source(paste0(scripts_dir, "load_replicates.R"))
  tryCatch(
    {
      # load 'CX_reports'
      n.cores.load <- ifelse(n.cores > 1, 2, 1)
      load_vars <- mclapply(var_args, function(x) load_replicates(x$path, n.cores, x$name, F, methyl_files_type), mc.cores = n.cores.load)

      # trimm seqs objects (rename if not 'TAIR10' Chr seqnames)
      meth_var1 <- trimm_and_rename(load_vars[[1]]$methylationData_pool)
      meth_var2 <- trimm_and_rename(load_vars[[2]]$methylationData_pool)
      meth_var1_replicates <- trimm_and_rename(load_vars[[1]]$methylationDataReplicates)
      meth_var2_replicates <- trimm_and_rename(load_vars[[2]]$methylationDataReplicates)

      if (!is_single) {
        # join replicates
        methylationDataReplicates_joints <- joinReplicates(
          meth_var1_replicates,
          meth_var2_replicates
        )
        message("load and join replicates data: successfully")
      } else {
        # join singles
        methylationDataReplicates_joints <- joinReplicates(
          meth_var1,
          meth_var2
        )
        message("load and join single-samples data: successfully")
      }
    },
    error = function(cond) {
      cat("\n*\n", as.character(cond), "\n*\n")
      stop("load and join CX methylation data: fail")
    }
  )

  ###########################################################################

  # new folders path names
  comparison_name <- paste0(var2, "_vs_", var1)
  exp_path <- paste0(Methylome.At_path, "/results/", comparison_name)
  ChrPlot_CX_path <- paste0(exp_path, "/ChrPlot_CX")
  ChrPlot_subCX_path <- paste0(exp_path, "/ChrPlot_subCX")
  PCA_plots_path <- paste0(exp_path, "/PCA_plots")
  meth_levels_path <- paste0(exp_path, "/methylation_levels")
  metaPlot_path <- paste0(exp_path, "/metaPlots")
  gainORloss_path <- paste0(exp_path, "/gain_OR_loss")
  genome_ann_path <- paste0(exp_path, "/genome_annotation")
  DMRs_bigWig_path <- paste0(exp_path, "/DMRs_bigWig")
  ChrPlots_DMRs_path <- paste0(exp_path, "/ChrPlot_DMRs")

  dir.create(exp_path, showWarnings = F)
  setwd(exp_path)

  ###########################################################################

  ##### calculate the conversion rate by the chloroplast chromosome (ChrC)
  message("\nconversion rate (C->T) along the Chloroplast genome:", appendLF = F)
  tryCatch(
    {
      message("")
      source(paste0(scripts_dir, "ChrC_conversionRate.R"))
      conR_var1 <- conversionRate(load_vars[[1]]$methylationDataReplicates, var1)
      conR_var2 <- conversionRate(load_vars[[2]]$methylationDataReplicates, var2)
      conR_b <- rbind(conR_var1, conR_var2)
      write.csv(conR_b, paste0(exp_path, "/conversion_rate.csv"), row.names = F)
      print(kable(conR_b))
    },
    error = function(cond) {
      cat("\n\n* conversion rate:\n", as.character(cond), "\n*\n")
      message("fail")
      cat(" fail\n")
    }
  )

  cat("\n")
  rm(load_vars)

  ##### PCA plot to total methlyation in all contexts
  if (!is_single) {
    dir.create(PCA_plots_path, showWarnings = F)
    setwd(PCA_plots_path)

    message(paste0("\ngenerating PCA plots of total methylation levels: "), appendLF = F)
    cat("\nPCA plots...")
    tryCatch(
      {
        source(paste0(scripts_dir, "pca_plot.R"))
        pca_plot(methylationDataReplicates_joints, var1, var2, var1_path, var2_path, "CG")
        pca_plot(methylationDataReplicates_joints, var1, var2, var1_path, var2_path, "CHG")
        pca_plot(methylationDataReplicates_joints, var1, var2, var1_path, var2_path, "CHH")
        pca_plot(methylationDataReplicates_joints, var1, var2, var1_path, var2_path, "all_contexts")
        message("done")
      },
      error = function(cond) {
        cat("\n\n* PCA plot:\n", as.character(cond), "\n*\n")
        message("fail")
      }
    )
    setwd(exp_path)
    cat(" done\n")
  } else {
    message("\n* skipping PCA plots for single-samples data")
    cat(" fail\n")
  }
  ###########################################################################

  ##### calculate and plot total methylation levels (%)
  dir.create(meth_levels_path, showWarnings = F)
  setwd(meth_levels_path)

  message("plotting total methylation levels (5-mC%): ", appendLF = F)
  cat("total methylation levels...")
  tryCatch(
    {
      source(paste0(scripts_dir, "total_meth_levels.R"))
      total_meth_levels <- total_meth_levels(meth_var1_replicates, meth_var2_replicates, var1, var2)
      message("done")
      cat(" done\n")
    },
    error = function(cond) {
      cat("\n\n* total methylation levels:\n", as.character(cond), "\n*\n")
      message("fail")
      cat(" fail\n")
    }
  )

  ##### ChrPlots for CX methylation #####
  dir.create(ChrPlot_CX_path, showWarnings = F)
  setwd(ChrPlot_CX_path)
  source(paste0(scripts_dir, "ChrPlots_CX.R"))

  message("generating chromosome methylation plots (ChrPlots): ", appendLF = F)
  tryCatch(
    {
      # suppressWarnings(ChrPlots_CX(comparison_name, meth_var1, meth_var2, var1, var2, scripts_dir, n.cores))
      # run_ChrPlots_CX(comparison_name, meth_var1, meth_var2, var1, var2, TE_file)
      run_ChrPlots_CX(var1, var2, meth_var1, meth_var2, TE_file, n.cores)
      message("done")
    },
    error = function(cond) {
      cat("\n\n* ChrPlots:\n", as.character(cond), "\n*\n")
      message("fail")
    }
  )

  setwd(exp_path)

  ##### ChrPlots for sub-CX methylation #####
  dir.create(ChrPlot_subCX_path, showWarnings = F)
  setwd(ChrPlot_subCX_path)

  message("generating chromosome methylation plots to sub-contexts (ChrPlots): ", appendLF = F)
  tryCatch(
    {
      source(paste0(scripts_dir, "ChrPlots_sub_CX.R"))
      chr_names <- unique(as.character(seqnames(annotation.gr)))
      run_ChrPlots_sub_CX(var1, var2, meth_var1, meth_var2, n.cores, chr_names)
      message("done\n")
    },
    error = function(cond) {
      cat("\n\n* sub-context ChrPlots:\n", as.character(cond), "\n*\n")
      message("fail\n")
    }
  )

  setwd(exp_path)

  ###########################################################################

  ##### call DMRs for replicates/single data
  message(paste0("call DMRs for replicates data: ", is_Replicates))

  if (is_Replicates) {
    message(paste0("compute ", binSize, "bp DMRs using: beta regression\n"))
  } else {
    message(paste0("compute ", binSize, "bp DMRs using: fisher’s exact test\n"))
  }

  ###########################################################################

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### main loop for 'DMRs' and its downstream results ### ###

  for (context in c("CG", "CHG", "CHH")) {
    cat("\n*\n", context, "\n*\n")
    message(paste0("DMRs in ", context, " context..."))

    ##############################
    ##### Calling DMRs in Replicates #####
    source(paste0(scripts_dir, "calling_DMRs.R"))
    DMRs_bins <- calling_DMRs(
      methylationDataReplicates_joints, meth_var1, meth_var2,
      var1, var2, var1_path, var2_path, comparison_name,
      context, minProportionDiff, binSize, pValueThreshold,
      minCytosinesCount, minReadsPerCytosine, n.cores, is_Replicates
    )
    message(paste0("\tstatistically significant DMRs: ", length(DMRs_bins)))
    message(paste0("\tDMRs caller in ", context, " context: done"))

    ##############################
    #####  Gain or Loss - DMRs #####
    dir.create(gainORloss_path, showWarnings = FALSE)
    setwd(gainORloss_path)
    source(paste0(scripts_dir, "gainORloss.R"))

    ##### Pie chart
    tryCatch(
      {
        gainORloss(DMRs_bins, context, add_count = T)
        message("\tpie chart (gain or loss): done")
      },
      error = function(cond) {
        cat("\n*\n", as.character(cond), "\n*\n")
        message("\tpie chart (gain or loss): fail\n")
      }
    )

    ##### Ratio distribution
    tryCatch(
      {
        ratio.distribution(DMRs_bins, var1, var2, context, comparison_name)
        message("\tratio distribution (gain or loss): done")
      },
      error = function(cond) {
        cat("\n*\n", as.character(cond), "\n*\n")
        message("\tratio distribution (gain or loss): fail\n")
      }
    )
    ##############################

    ##### ChrPlots for DMRs #####
    tryCatch(
      {
        dir.create(ChrPlots_DMRs_path, showWarnings = FALSE)
        setwd(ChrPlots_DMRs_path)
        source(paste0(scripts_dir, "ChrPlots_DMRs.R"))
        ChrPlots_DMRs(comparison_name, DMRs_bins, var1, var2, context, scripts_dir)
        message("\tgenerated ChrPlots for all DMRs: done")
      },
      error = function(cond) {
        cat("\n*\n", as.character(cond), "\n*\n")
        message("\tgenerated ChrPlots for all DMRs: fail\n")
      }
    )
    setwd(exp_path)

    ##### Annotate DMRs and total-methylations #####
    message("\tgenome annotations for DMRs...")
    source(paste0(scripts_dir, "genome_ann.R"))
    source(paste0(scripts_dir, "DMRs_ann.R"))
    source(paste0(scripts_dir, "CX_ann.R"))
    source(paste0(scripts_dir, "DMRs_ann_plots.R"))
    source(paste0(scripts_dir, "TE_ann_plots.R"))
    source(paste0(scripts_dir, "TE_Super_Family_Frequency.R"))
    dir.create(genome_ann_path, showWarnings = FALSE)
    setwd(genome_ann_path)

    # genome annotations
    tryCatch(
      {
        ann_list <- genome_ann(annotation.gr, TE_file) # create annotations from annotation file as a list
        DMRs_ann(ann_list, DMRs_bins, context, description_df) # save tables of annotate DMRs
        CX_ann(ann_list, var1, var2, meth_var1, meth_var2, context) # save tables of annotate CX
        DMRs_ann_plots(var1, var2, context)
        message("\tgenome annotations for DMRs: done")
      },
      error = function(cond) {
        cat("\n*\n", as.character(cond), "\n*\n")
        message("\tgenome annotations for DMRs: fail")
      }
    )

    # additional TE annotations results
    tryCatch(
      {
        TE_ann_plots(context, TE_file)
        TE_Super_Family_Frequency(context, TE_file)
      },
      error = function(cond) {
        cat("\n*\n", as.character(cond), "\n*\n")
        message("\tTE families plots: fail")
      }
    )
    setwd(exp_path)

    ##### save DMRs as bigWig file #####
    dir.create(DMRs_bigWig_path, showWarnings = FALSE)
    source(paste0(scripts_dir, "DMRs_2_bigWig.R"))
    ann_res_files <- list.files(paste0(genome_ann_path, "/", context))
    for (ann.loop.bigWig in c("all", ann_res_files)) {
      suppressWarnings(try(
        {
          DMRs_2_bigWig(var1, var2, context, ann.loop.bigWig)
        },
        silent = T
      ))
    }
    message("\tsaved all DMRs also as bigWig files\n")
  }

  ### ### # finish main loop  ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  ###########################################################################

  ##### DMRs density - curcular plot #####
  tryCatch(
    {
      setwd(exp_path)
      # setwd(ChrPlots_DMRs_path)
      source(paste0(scripts_dir, "circular_plot_DMRs.R"))
      DMRs_circular_plot(annotation.gr, TE_file, comparison_name)
      DMRs_circular_plot_legends()
      message("generated DMRs density plot for all contexts: done\n")
    },
    error = function(cond) {
      cat("\n*\n", as.character(cond), "\n*\n")
      message("generated DMRs density plot for all contexts: fail\n")
    }
  )

  ###########################################################################

  source(paste0(scripts_dir, "multiplot_ggplot2.R"))

  ##### GO analysis for annotated DMRs #####
  if (GO_analysis) {
    tryCatch(
      {
        source(paste0(scripts_dir, "go_script.R"))
        GO_path <- paste0(exp_path, "/GO_analysis")
        dir.create(GO_path, showWarnings = FALSE)

        message("GO analysis for annotated DMRs...")
        run_GO(comparison_name, genome_ann_path, GO_path, n.cores)
        message("GO analysis for annotated DMRs: done\n")
      },
      error = function(cond) {
        cat("\n*\n", as.character(cond), "\n*\n")
        message("GO analysis for annotated DMRs: fail\n")
      }
    )
  }

  ##### KEGG pathways for annotated DMRs #####
  if (KEGG_pathways) {
    tryCatch(
      {
        source(paste0(scripts_dir, "kegg_script.R"))
        KEGG_path <- paste0(exp_path, "/KEGG_pathway")
        dir.create(KEGG_path, showWarnings = FALSE)

        message("KEGG pathways for annotated DMRs...")
        run_KEGG(comparison_name, genome_ann_path, KEGG_path, n.cores)
        message("KEGG pathways for annotated DMRs: done\n")
      },
      error = function(cond) {
        cat("\n*\n", as.character(cond), "\n*\n")
        message("KEGG pathways for annotated DMRs: fail\n")
      }
    )
  }

  ###########################################################################

  setwd(Methylome.At_path)
  message(paste0("**\t", var2, " vs ", var1, ": done\n"))

  ###########################################################################

  # date and time of the end
  message(paste0(
    "\n**\t",
    paste(format(Sys.time(), "%d"), format(Sys.time(), "%m"), format(Sys.time(), "%Y"), sep = "-"),
    " ", format(Sys.time(), "%H:%M")
  ))
  end_time <- Sys.time()
  time_diff <- difftime(end_time, start_time, units = "mins") %>% as.numeric()
  end_msg <- paste(
    "**\ttime:",
    paste0(floor(time_diff / 60), "hr"),
    paste0(floor(time_diff %% 60), "min\n")
  )
  cat(paste0("\n\n**\tDone!\n", end_msg))
  message(end_msg)
}
