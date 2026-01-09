Methylome.At_main <- function(var1, # control
                              var2, # treatment
                              var1_path,
                              var2_path,
                              Methylome.At_path = ".",
                              annotation_file = "./annotation_files/Methylome.At_annotations.csv.gz",
                              description_file = "./annotation_files/Methylome.At_annotations.csv.gz",
                              TEs_file = "./annotation_files/Methylome.At_annotations.csv.gz",
                              minProportionDiff = c(0.4, 0.2, 0.1), # CG, CHG, CHH
                              binSize = 100,
                              minCytosinesCount = 4,
                              minReadsPerCytosine = 4,
                              pValueThreshold = 0.05,
                              methyl_files_type = "CX_report",
                              img_type = "pdf",
                              n.cores = 8,
                              GO_analysis = FALSE,
                              KEGG_pathways = FALSE,
                              analyze_dH = FALSE,
                              TE_metaPlots = FALSE,
                              GeneBody_metaPlots = FALSE,
                              GeneFeatures_metaPlots = FALSE,
                              gene_features_binSize = 10,
                              metaPlot.random.genes = 10000) {
  ###########################################################################

  start_time <- Sys.time()
  sep_cat <- paste0("\n", paste(rep("-", 50), collapse = ""), "\n")
  scripts_dir <- paste0(Methylome.At_path, "/scripts")

  # source all R scripts
  script_files <- list.files(scripts_dir, pattern = "\\.R$", full.names = TRUE)
  script_files <- script_files[!grepl("install_R_packages\\.R$", script_files)]
  script_files <- script_files[!grepl("Methylome\\.At_run\\.R$", script_files)]
  script_files <- script_files[!grepl("Methylome\\.At_main\\.R$", script_files)]
  script_files <- script_files[!grepl("MetaPlots_run\\.R$", script_files)]
  script_files <- script_files[!grepl("mean_deltaH_CX\\.R", script_files)] # have match functions with 'ChrPlots' functions
  script_files <- script_files[!grepl("ChrPlots_", script_files)]
  invisible(lapply(script_files, source))

  ###########################################################################

  # image device function
  img_device <<- function(filename, w, h) {
    dev_call <- ifelse(img_type == "pdf", "cairo_pdf", img_type)
    img_env <- get(dev_call, envir = asNamespace("grDevices"))
    full_file_name <- paste0(filename, ".", img_type)

    if (img_type == "svg" | img_type == "pdf") {
      img_env(full_file_name, width = w, height = h, family = "serif")
    } else if (img_type == "tiff") {
      img_env(full_file_name, width = w, height = h, units = "in", res = 900, family = "serif", compression = "lzw")
    } else {
      img_env(full_file_name, width = w, height = h, units = "in", res = 900, family = "serif")
    }
  }

  ###########################################################################

  setwd(Methylome.At_path)

  ##### read annotation and description files #####
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
      cat("\n*\n load annotation file:\n", as.character(cond), "*\n")
      message("load 'annotation' file: fail")
    }
  )
  cat("\rload annotations and description files [1/3]")

  # TAIR10 Transposable Elements file
  tryCatch(
    {
      TE_file.df <- read.csv(TEs_file, sep = "\t")
      TE_file <- edit_TE_file(TE_file.df)
      message("load Transposable Elements file")
    },
    error = function(cond) {
      cat("\n*\n load TE file:\n", as.character(cond), "*\n")
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
      cat("\n*\n load description file:\n", as.character(cond), "*\n")
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
  tryCatch(
    {
      # load 'CX_reports'
      n.cores.load <- ifelse(n.cores > 1, 2, 1)
      load_vars <- mclapply(var_args, function(x) {
        Sys.sleep(ifelse(x$name == var1, 0, 6))
        load_replicates(x$path, n.cores, x$name, F, methyl_files_type)
      },
      mc.cores = n.cores.load
      )

      # trimm seqs objects (rename if not 'TAIR10' Chr seqnames)
      meth_var1 <- trimm_and_rename(load_vars[[1]]$methylationData_pool)
      meth_var2 <- trimm_and_rename(load_vars[[2]]$methylationData_pool)
      meth_var1_replicates <- trimm_and_rename(load_vars[[1]]$methylationDataReplicates)
      meth_var2_replicates <- trimm_and_rename(load_vars[[2]]$methylationDataReplicates)
      cat("\nPooled ", var1, " object for example:\n\n", sep = "")
      capture.output(print(meth_var1), file = NULL)[-c(1, 16)] %>% cat(sep = "\n")
      cat(paste("seq-levels:", paste(seqlevels(meth_var1), collapse = " ")), "\n\n")

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
      cat("\n*\n load and join CX methylation data:\n", as.character(cond), "*\n")
      stop("load and join CX methylation data: fail")
    }
  )

  ###########################################################################

  # new folders path names
  comparison_name <- paste0(var2, "_vs_", var1)
  exp_path <- paste0(Methylome.At_path, "/results/", comparison_name)
  ChrPlot_CX_path <- paste0(exp_path, "/ChrPlot_CX")
  ChrPlot_subCX_path <- paste0(exp_path, "/ChrPlot_CX/subCX")
  PCA_plots_path <- paste0(exp_path, "/PCA_plots")
  meth_levels_path <- paste0(exp_path, "/methylation_levels")
  gainORloss_path <- paste0(exp_path, "/gain_OR_loss")
  genome_ann_path <- paste0(exp_path, "/genome_annotation")
  DMRs_bigWig_path <- paste0(exp_path, "/DMRs_bigWig")
  ChrPlots_DMRs_path <- paste0(exp_path, "/ChrPlot_DMRs")
  dH_CX_path <- paste0(exp_path, "/deltaH")
  dH_CX_ann_path <- paste0(exp_path, "/deltaH/genome_annotation")
  metaPlot_path <- paste0(exp_path, "/metaPlots")

  dir.create(exp_path, showWarnings = F)
  setwd(exp_path)

  ###########################################################################

  cat(sep_cat)

  ##### calculate the conversion rate by the chloroplast chromosome (ChrC)
  message("\nconversion rate (C->T) along the Chloroplast genome:", appendLF = F)
  cat("\nconversion rate (C->T) along the Chloroplast genome:") # "\n"
  tryCatch(
    {
      message("")
      conR_var1 <- conversionRate(load_vars[[1]]$methylationDataReplicates, var1)
      conR_var2 <- conversionRate(load_vars[[2]]$methylationDataReplicates, var2)
      conR_b <- rbind(conR_var1, conR_var2)
      write.csv(conR_b, paste0(exp_path, "/conversion_rate.csv"), row.names = F)
      print(kable(conR_b))
    },
    error = function(cond) {
      cat("\n*\n conversion rate:\n", as.character(cond), "*\n")
      message("fail")
      cat(" fail\n")
    }
  )

  cat(sep_cat)
  rm(load_vars)

  ##### PCA plot to total methlyation in all contexts
  if (!is_single) {
    dir.create(PCA_plots_path, showWarnings = F)
    setwd(PCA_plots_path)

    message(paste0("\ngenerating PCA plots of total methylation levels: "), appendLF = F)
    cat("\nPCA plots...")
    tryCatch(
      {
        pca_plot(methylationDataReplicates_joints, var1, var2, var1_path, var2_path, "CG")
        pca_plot(methylationDataReplicates_joints, var1, var2, var1_path, var2_path, "CHG")
        pca_plot(methylationDataReplicates_joints, var1, var2, var1_path, var2_path, "CHH")
        pca_plot(methylationDataReplicates_joints, var1, var2, var1_path, var2_path, "all_contexts")
        message("done")
      },
      error = function(cond) {
        cat("\n*\n PCA plot:\n", as.character(cond), "*\n")
        message("fail")
      }
    )
    setwd(exp_path)
    cat(" done\n")
  } else {
    message("\n* skipping PCA plots for single-samples data")
    cat("skipping PCA plots for single-samples data\n")
  }

  ##### calculate and plot total methylation levels (%)
  dir.create(meth_levels_path, showWarnings = F)
  setwd(meth_levels_path)

  message("plotting total methylation levels (5-mC%): ", appendLF = F)
  cat("total methylation levels...")
  tryCatch(
    {
      total_meth_levels <- total_meth_levels(meth_var1_replicates, meth_var2_replicates, var1, var2)
      message("done")
      cat(" done\n")
    },
    error = function(cond) {
      cat("\n*\n total methylation levels:\n", as.character(cond), "*\n")
      message("fail")
      cat(" fail\n")
    }
  )

  ###########################################################################

  ##### ChrPlots for CX methylation #####
  dir.create(ChrPlot_CX_path, showWarnings = F)
  dir.create(ChrPlot_subCX_path, showWarnings = F)
  setwd(ChrPlot_CX_path)

  cat("\n* chromosome methylation plots (ChrPlots):")
  message("generating chromosome methylation plots (ChrPlots): ", appendLF = F)
  tryCatch(
    {
      source(paste0(scripts_dir, "/ChrPlots_CX.R"))
      suppressWarnings(run_ChrPlots_CX(var1, var2, meth_var1, meth_var2, TE_file, n.cores))
      message("done")
    },
    error = function(cond) {
      cat("\n*\n ChrPlots:\n", as.character(cond), "*\n")
      message("fail")
    }
  )

  setwd(exp_path)

  ###########################################################################

  cat(sep_cat)

  ##### call DMRs for replicates/single data
  cat("\ncall DMRs for replicates data:", is_Replicates, "\n")
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
    cat("\n-------------------------\n")
    cat("\ncalculate DMRs in", context, "context...\n")
    message(paste0("DMRs in ", context, " context..."))

    ##############################
    ##### Calling DMRs in Replicates #####
    tryCatch(
      {
        DMRs_bins <- calling_DMRs(
          methylationDataReplicates_joints, meth_var1, meth_var2,
          var1, var2, var1_path, var2_path, comparison_name,
          context, minProportionDiff, binSize, pValueThreshold,
          minCytosinesCount, minReadsPerCytosine, n.cores, is_Replicates
        )
        cat(paste0("statistically significant DMRs: ", length(DMRs_bins), "\nDMRs plots...\n"))
        message(paste0("\tstatistically significant DMRs: ", length(DMRs_bins)))
        message(paste0("\tDMRs caller in ", context, " context: done"))
      },
      error = function(cond) {
        cat(paste0("\n*\n Calling DMRs in", context, ":\n"), as.character(cond), "*\ncontinue without calling DMRs!\n\n")
        message("\tCalling DMRs: fail\n")
        GO_analysis <- FALSE
        KEGG_pathways <- FALSE
        break
      }
    )

    ##############################
    #####  Gain or Loss - DMRs #####
    dir.create(gainORloss_path, showWarnings = FALSE)
    setwd(gainORloss_path)

    ##### Pie chart
    tryCatch(
      {
        gainORloss(DMRs_bins, context, add_count = T)
        message("\tpie chart (gain or loss): done")
      },
      error = function(cond) {
        cat("\n*\n pie chart (gain or loss):\n", as.character(cond), "*\n")
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
        cat("\n*\n ratio distribution (gain or loss):\n", as.character(cond), "*\n")
        message("\tratio distribution (gain or loss): fail\n")
      }
    )
    ##############################

    ##### ChrPlots for DMRs #####
    tryCatch(
      {
        source(paste0(scripts_dir, "/ChrPlots_DMRs.R"))
        dir.create(ChrPlots_DMRs_path, showWarnings = FALSE)
        setwd(ChrPlots_DMRs_path)
        ChrPlots_DMRs(comparison_name, DMRs_bins, var1, var2, context, scripts_dir)
        message("\tgenerated ChrPlots for all DMRs: done")
      },
      error = function(cond) {
        cat("\n*\n tgenerated ChrPlots for all DMRs:\n", as.character(cond), "*\n")
        message("\tgenerated ChrPlots for all DMRs: fail\n")
      }
    )
    setwd(exp_path)

    ##### Annotate DMRs and total-methylations #####
    cat("genome annotations for DMRs...\n")
    message("\tgenome annotations for DMRs...")
    dir.create(genome_ann_path, showWarnings = FALSE)
    setwd(genome_ann_path)

    # genome annotations
    tryCatch(
      {
        ann_list <- genome_ann(annotation.gr, TE_file) # create annotations from annotation file as a list
        DMRs_ann(ann_list, DMRs_bins, context, description_df) # save tables of annotate DMRs. have to run after 'genome_ann'
        CX_ann(ann_list, var1, var2, meth_var1, meth_var2, context) # save tables of annotate CX
        DMRs_ann_plots(var1, var2, context)
        message("\tgenome annotations for DMRs: done")
      },
      error = function(cond) {
        cat("\n*\n tgenome annotations for DMRs:\n", as.character(cond), "*\n")
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
        cat("\n*\n TE families plots:\n", as.character(cond), "*\n")
        message("\tTE families plots: fail")
      }
    )
    setwd(exp_path)

    ##### save DMRs as bigWig file #####
    dir.create(DMRs_bigWig_path, showWarnings = FALSE)
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
    cat("done\n")
  }

  setwd(exp_path)

  ### ### # finish main loop  ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  ###########################################################################

  cat(sep_cat)

  ##### DMRs density - curcular plot #####
  tryCatch(
    {
      setwd(exp_path)
      cat("\ngenerated DMRs density plot for all contexts: ")
      # setwd(ChrPlots_DMRs_path)
      DMRs_circular_plot(annotation.gr, TE_file, comparison_name)
      DMRs_circular_plot_legends()
      cat("done\n")
      message("generated DMRs density plot for all contexts: done\n")
    },
    error = function(cond) {
      cat("\n*\n DMRs density plot:\n", as.character(cond), "*\n")
      message("generated DMRs density plot for all contexts: fail\n")
    }
  )

  ###########################################################################

  ##### TEs superfamily - curcular plot #####
  tryCatch(
    {
      setwd(genome_ann_path)
      cat("\ngenerated DMRs density plots over TEs superfamilies: ")
      TEs_superfamily_circular_plot(annotation.gr)
      cat("done\n")
      message("generated DMRs over TEs superfamilies: done\n")
    },
    error = function(cond) {
      cat("\n*\n DMRs over TEs superfamilies:\n", as.character(cond), "*\n")
      message("generated DMRs over TEs superfamilies: fail\n")
    }
  )

  setwd(exp_path)

  ###########################################################################

  ##### DMRs within Genes-groups #####
  func_groups_path <- paste0(genome_ann_path, "/functional_groups")
  setwd(exp_path)
  dir.create(func_groups_path, showWarnings = FALSE)
  cat("Annotate DMRs into functional groups... ")
  for (ann.l in c("Genes", "Promoters")) {
    groups_results <- c()
    for (cntx.l in c("CG", "CHG", "CHH", "all")) {
      tryCatch(
        {
          cat(".")
          groups_results <- rbind(
            groups_results,
            DMRs_into_groups(treatment = comparison_name, ann = ann.l, context = cntx.l)
          )
        },
        error = function(cond) {
          cat("\n*\n Annotate *", cntx.l, "* - *", ann.l, "* DMRs into functional groups:\n", as.character(cond), "*\n")
          message(paste("Annotate", cntx.l, "-", ann.l, "DMRs into functional groups: fail\n"))
        }
      )
    }
    tryCatch(
      {
        img_device(paste0(func_groups_path, "/", ann.l, "_groups_barPlots_", comparison_name), w = 14, h = 4.5)
        print(groups_barPlots(groups_results))
        dev.off()
        cat(" done\n")
      },
      error = function(cond) {
        cat("\n*\n Bar-plot of *", cntx.l, "* - *", ann.l, "* annotated DMRs into functional groups:\n", as.character(cond), "*\n")
      }
    )
  }
  message("Annotate DMRs into functional groups: done\n")

  ###########################################################################

  ##### TEs - size and distance plots
  tryCatch(
    {
      TE_context_list <- TE_delta_meth(list(meth_var1, meth_var2), TE_file)

      dir.create(paste0(genome_ann_path, "/TEs_addiotionnal_results/TE_size_n_distance/"), showWarnings = FALSE)

      ## TE methylation levels (delta) and size
      cat("\nTE delta-methylation vs. TE-size\n")
      for (te_sz_cntx in c("CG", "CHG", "CHH")) {
        ggsave(
          filename = paste0(genome_ann_path, "/TEs_addiotionnal_results/TE_size_n_distance/", te_sz_cntx, "_TE_size_delta_scatter.png"),
          plot = te_size_plot(TE_context_list, te_sz_cntx),
          width = 3750,
          height = 2500,
          units = "px",
          dpi = 1200
        )
      }
      message("TE delta-methylation vs. TE-size: done")

      ## TE methylation levels (delta) and distance from centromer
      cat("TE delta-methylation vs. distance from centromere\n")
      TE_distance <- distance_from_centromer(TE_context_list, TE_file, window_size = 1e6)

      ggsave(
        filename = paste0(genome_ann_path, "/TEs_addiotionnal_results/TE_size_n_distance/TE_centromere_distance_delta.png"),
        plot = TE_distance$plot,
        width = 3750,
        height = 2500,
        units = "px",
        dpi = 1200
      )
      message("TE delta-methylation vs. distance from centromer: done\n")
    },
    error = function(cond) {
      cat("\n*\n TEs size and distance plots:\n", as.character(cond), "*\n")
    }
  )

  ###########################################################################


  ##### GO analysis for annotated DMRs
  if (GO_analysis) {
    cat(sep_cat)
    tryCatch(
      {
        GO_path <- paste0(exp_path, "/GO_analysis")
        dir.create(GO_path, showWarnings = FALSE)

        message("GO analysis for annotated DMRs...")
        run_GO(comparison_name, genome_ann_path, GO_path, n.cores)
        message("GO analysis for annotated DMRs: done\n")
      },
      error = function(cond) {
        cat("\n*\n GO analysis:\n", as.character(cond), "*\n")
        message("GO analysis for annotated DMRs: fail\n")
      }
    )
  }

  ##### KEGG pathways for annotated DMRs
  if (KEGG_pathways) {
    cat(sep_cat)
    tryCatch(
      {
        KEGG_path <- paste0(exp_path, "/KEGG_pathway")
        dir.create(KEGG_path, showWarnings = FALSE)

        message("KEGG pathways for annotated DMRs...")
        run_KEGG(comparison_name, genome_ann_path, KEGG_path, n.cores)
        message("KEGG pathways for annotated DMRs: done\n")
      },
      error = function(cond) {
        cat("\n*\n KEGG pathways:\n", as.character(cond), "*\n")
        message("KEGG pathways for annotated DMRs: fail\n")
      }
    )
  }

  ###########################################################################

  ##### dH analysis over CX methylation #####
  if (analyze_dH) {
    cat(sep_cat)
    dir.create(dH_CX_path, showWarnings = F)
    dir.create(dH_CX_ann_path, showWarnings = F)
    setwd(dH_CX_path)

    cat("\n* chromosome dH plots (ChrPlots):")
    message("generating ChrPlot (mean of dH) and scatter-plot (dH vs dmC): ", appendLF = F)
    tryCatch(
      {
        source(paste0(scripts_dir, "/mean_deltaH_CX.R"))
        suppressWarnings(run_mean_deltaH_CX(var1, var2, meth_var1, meth_var2, TE_file, n.cores))
        message("done")
      },
      error = function(cond) {
        cat("\n*\n mean dH:\n", as.character(cond), "*\n")
        message("fail")
      }
    )

    message("generating sum dH analysis:\n", appendLF = F) # , rep("-", 29)
    tryCatch(
      {
        suppressWarnings(run_sum_deltaH_CX(var1, var2, meth_var1, meth_var2, annotation.gr, TE_file, description_df, n.cores, fdr = 0.95))
      },
      error = function(cond) {
        cat("\n*\n sum dH:\n", as.character(cond), "*\n")
        message("fail\n")
      }
    )
  }

  setwd(exp_path)

  ###########################################################################

  ##### run metPlot function for coding-Genes and TEs
  if (TE_metaPlots | GeneBody_metaPlots | GeneFeatures_metaPlots) {
    cat(sep_cat, "\nmetaPlots:\n----------")
    dir.create(metaPlot_path, showWarnings = F)

    # calculate metaPlot for genes bodies
    if (TE_metaPlots) {
      tryCatch(
        {
          message(paste("generate metaPlot to", metaPlot.random.genes, "protein-coding Genes..."))
          setwd(metaPlot_path)
          Genes_metaPlot(meth_var1, meth_var2, var1, var2, annotation.gr, metaPlot.random.genes, minReadsPerCytosine, n.cores, is_TE = F)
          setwd(metaPlot_path)
          delta_metaplot("Genes", var1, var2)
        },
        error = function(cond) {
          cat("\n*\n TEs metaPlots:\n", as.character(cond), "*\n")
          message(paste0("process average metaPlot to ", metaPlot.random.genes, " Protein Coding Genes: fail"))
        }
      )
    }

    # calculate metaPlot for TEs
    if (GeneBody_metaPlots) {
      tryCatch(
        {
          message(paste("\ngenerate metaPlot to", metaPlot.random.genes, " Transposable Elements..."))
          setwd(metaPlot_path)
          Genes_metaPlot(meth_var1, meth_var2, var1, var2, TE_file, metaPlot.random.genes, minReadsPerCytosine, n.cores, is_TE = T)
          setwd(metaPlot_path)
          delta_metaplot("TEs", var1, var2)
        },
        error = function(cond) {
          cat("\n*\n Transposable Elements metaPlots:\n", as.character(cond), "*\n")
          message(paste0("process average metaPlot to ", metaPlot.random.genes, " Transposable Elements: fail\n"))
        }
      )
    }

    # calculate metaPlot for coding-Gene features (CDS, introns, etc.)
    if (GeneFeatures_metaPlots) {
      tryCatch(
        {
          message(paste("\ngenerate metaPlot to", metaPlot.random.genes, "protein-coding Gene Features..."))
          setwd(metaPlot_path)
          Genes_features_metaPlot(meth_var1, meth_var2, var1, var2, annotation.gr, metaPlot.random.genes, minReadsPerCytosine, gene_features_binSize, n.cores)
          # delta_metaplot("Gene_features", var1, var2, is_geneFeature = TRUE)
        },
        error = function(cond) {
          cat("\n*\n Gene features metaPlots:\n", as.character(cond), "*\n")
          message(paste0("process average metaPlot to ", metaPlot.random.genes, " Protein Coding Gene Features: fail"))
        }
      )
    }
  }

  setwd(exp_path)

  ###########################################################################

  setwd(Methylome.At_path)
  message(paste0("**\t", var2, " vs ", var1, ": done\n"))
  cat(sep_cat)

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
