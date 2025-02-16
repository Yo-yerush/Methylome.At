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
                              metaPlot.random.genes,
                              pValueThreshold,
                              n.cores,
                              GO_analysis,
                              KEGG_pathways)
{
  ########################################################################### 
  start_time <- Sys.time()
  
  # scripts directory path
  scripts_dir = paste0(Methylome.At_path,"/scripts/")
  
  ########################################################################### 
  
  ##### load the data for replicates ##### 
  message("load replicates data...")
  source(paste0(scripts_dir,"load_replicates.R"))
  tryCatch({
    # load 'CX_reports'
    if (n.cores > 1) {n.cores.load = 2} else {n.cores.load = 1}
    load_vars = mclapply(list(var1_path,var2_path), function(x) load_replicates(x, n.cores), mc.cores = n.cores.load)
    
    # trimm seqs objects (rename if not 'TAIR10' Chr seqnames)
    source(paste0(scripts_dir,"trimm_and_rename_seq.R"))
    meth_var1 = trimm_and_rename(load_vars[[1]]$methylationData_pool)
    meth_var2 = trimm_and_rename(load_vars[[2]]$methylationData_pool)
    meth_var1_replicates = trimm_and_rename(load_vars[[1]]$methylationDataReplicates)
    meth_var2_replicates = trimm_and_rename(load_vars[[2]]$methylationDataReplicates)
    
    # join replicates
    message("join replicates data...")
    methylationDataReplicates_joints = joinReplicates(meth_var1_replicates,
                                                      meth_var2_replicates)
    message("load and join replicates data: successfully") 
  },
  error = function(cond) {
    cat("\n*\n",as.character(cond),"\n*\n")
    stop("load and join replicates data: fail")
  })
  
  ###########################################################################
  
  # results directory path
  Methylome.At_res_dir = paste0(Methylome.At_path,"/results/")
  
  # new folders path names
  comparison_name = paste0(var2,"_vs_",var1)
  exp_path = paste0(Methylome.At_res_dir,comparison_name)
  ChrPlot_CX_path = paste0(Methylome.At_res_dir,comparison_name,"/ChrPlot_CX")
  meth_levels_path = paste0(Methylome.At_res_dir,comparison_name,"/methylation_levels")
  metaPlot_path = paste0(Methylome.At_res_dir,comparison_name,"/metaPlots")
  gainORloss_path = paste0(Methylome.At_res_dir,comparison_name,"/gain_OR_loss")
  genome_ann_path = paste0(Methylome.At_res_dir,comparison_name,"/genome_annotation")
  DMRs_bedGragh_path = paste0(Methylome.At_res_dir,comparison_name,"/DMRs_bedGragh")
  ChrPlots_DMRs_path = paste0(Methylome.At_res_dir,comparison_name,"/ChrPlot_DMRs")
  
  dir.create(exp_path)
  setwd(exp_path)
  
  ###########################################################################
  
  # calculate the conversion rate by the chloroplast chromosome (ChrC)
  message("\nconversion rate (C->T) along the Chloroplast genome:")
  source(paste0(scripts_dir,"ChrC_conversionRate.R"))
  conR_var1 = conversionRate(load_vars[[1]]$methylationDataReplicates, var1)
  conR_var2 = conversionRate(load_vars[[2]]$methylationDataReplicates, var2)
  write.csv(rbind(conR_var1,conR_var2), paste0(exp_path,"/conversion_rate.csv"), row.names = F)
  
  ###########################################################################
  
  ##### calculate and plot total methylation levels (%)
  dir.create(meth_levels_path)
  setwd(meth_levels_path)
  
  tryCatch({
    source(paste0(scripts_dir,"total_meth_levels.R"))
    total_meth_levels = total_meth_levels(meth_var1_replicates, meth_var2_replicates, var1, var2)
    message(paste0("\nplot total methylation levels (5-mC%)"))
  },
  error = function(cond) {
    cat("\n*\n",as.character(cond),"\n*\n")
    message("\nplot total methylation levels (5-mC%): fail")
  })
  
  ##### ChrPlots for CX methylation #####
  dir.create(ChrPlot_CX_path)
  setwd(ChrPlot_CX_path)
  
  tryCatch({
    source(paste0(scripts_dir,"ChrPlots_CX.R"))
    suppressWarnings(ChrPlots_CX(comparison_name, meth_var1, meth_var2, var1, var2, scripts_dir, n.cores))
    message("generated ChrPlots to methylation levels\n")
  },
  error = function(cond) {
    cat("\n*\n",as.character(cond),"\n*\n")
    message("generated ChrPlots to methylation levels: fail\n")
  })
  setwd(exp_path)
  
  ###########################################################################
  
  ##### import annotation and description files ##### 
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
  }, error = function(cond) {
    cat("\n*\n",as.character(cond),"\n*\n")
    message("import 'annotation' file: fail")
  })
  
  # TAIR10 Transposable Elements file
  tryCatch({
    source(paste0(scripts_dir,"edit_TE_file.R"))
    TE_file.df = read.csv(TEs_file, sep = "\t")
    TE_file = edit_TE_file(TE_file.df)
    message("import Transposable Elements file")
  }, error = function(cond) {
    cat("\n*\n",as.character(cond),"\n*\n")
    message("import Transposable Elements file: fail")
  })
  
  # upload description file
  tryCatch({
    des_file_sep = ifelse(grepl("\\.csv$|\\.csv\\.gz$",description_file), ",", "\t")
    description_df = read.csv(description_file, sep = des_file_sep)
    names(description_df)[1] = "gene_id"
    message("import description file\n")
  }, error = function(cond) {
    cat("\n*\n",as.character(cond),"\n*\n")
    message("import description file: fail\n")
  })
  
  setwd(exp_path)
  
  ###########################################################################
  
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### main loop for 'DMRs' and its downstream results ### ### 
  
  for (context in c("CG","CHG","CHH")) {
    cat("\n*\n",context,"\n*\n")
    message(paste0("DMRs in ",context ," context..."))
    
    ##############################
    ##### Calling DMRs in Replicates ##### 
    source(paste0(scripts_dir,"calling_DMRs.R"))
    DMRsReplicates = calling_DMRs(methylationDataReplicates_joints, var1, var2,
                                  var1_path, var2_path, comparison_name, context,
                                  minProportionDiff, binSize, pValueThreshold, 
                                  minCytosinesCount, minReadsPerCytosine, n.cores)
    message(paste0("\tstatistically significant DMRs: ", length(DMRsReplicates)))
    message(paste0("\tDMRs caller in ",context ," context: done"))
    
    ##############################
    #####  Gain or Loss - DMRs ##### 
    dir.create(gainORloss_path, showWarnings = FALSE)
    setwd(gainORloss_path)
    source(paste0(scripts_dir,"gainORloss.R"))
    
    ##### Pie chart
    tryCatch({
      gainORloss(DMRsReplicates, context)
      message("\tpie chart (gain or loss): done")
    }, error = function(cond) {
      cat("\n*\n",as.character(cond),"\n*\n")
      message("\tpie chart (gain or loss): fail\n")
    })
    
    ##### Ratio distribution 
    tryCatch({
      ratio.distribution(DMRsReplicates, var1, var2, context, comparison_name)
      message("\tratio distribution (gain or loss): done")
    }, error = function(cond) {
      cat("\n*\n",as.character(cond),"\n*\n")
      message("\tratio distribution (gain or loss): fail\n")
    })
    ##############################
    
    ##### ChrPlots for DMRs #####
    tryCatch({
      dir.create(ChrPlots_DMRs_path, showWarnings = FALSE)
      setwd(ChrPlots_DMRs_path)
      source(paste0(scripts_dir,"ChrPlots_DMRs.R"))
      ChrPlots_DMRs(comparison_name, DMRsReplicates, var1, var2, context, scripts_dir)
      message("\tgenerated ChrPlots for all DMRs: done")
    }, error = function(cond) {
      cat("\n*\n",as.character(cond),"\n*\n")
      message("\tgenerated ChrPlots for all DMRs: fail\n")
    })
    setwd(exp_path)
    
    ##### Annotate DMRs and total-methylations ##### 
    message("\tgenome annotations for DMRs...")
    source(paste0(scripts_dir,"genome_ann.R"))
    source(paste0(scripts_dir,"DMRs_ann.R"))
    source(paste0(scripts_dir,"CX_ann.R"))
    source(paste0(scripts_dir,"DMRs_ann_plots.R"))
    source(paste0(scripts_dir,"TE_ann_plots.R"))
    source(paste0(scripts_dir,"TE_Super_Family_Frequency.R"))
    dir.create(genome_ann_path, showWarnings = FALSE)
    setwd(genome_ann_path)
    
    # genome annotations
    tryCatch({
      ann_list = genome_ann(annotation.gr, TE_file) # create annotations from annotation file as a list
      DMRs_ann(ann_list, DMRsReplicates, context, description_df) # save tables of annotate DMRs
      CX_ann(ann_list, var1, var2, meth_var1, meth_var2, context) # save tables of annotate CX
      DMRs_ann_plots(var1, var2, context)
      message("\tgenome annotations for DMRs: done")
    }, 
    error = function(cond) {
      cat("\n*\n",as.character(cond),"\n*\n")
      message("\tgenome annotations for DMRs: fail")
    })
    
    # additional TE annotations results
    tryCatch({
      TE_ann_plots(context, TE_file)
      TE_Super_Family_Frequency(context, TE_file)
    }, 
    error = function(cond) {
      cat("\n*\n",as.character(cond),"\n*\n")
      message("\tTE families plots: fail")
    })
    setwd(exp_path)
    
    ##### save DMRs as bedGragh file ##### 
    dir.create(DMRs_bedGragh_path, showWarnings = FALSE)
    source(paste0(scripts_dir,"DMRs_2_bedGragh.R"))
    ann_res_files = list.files(paste0(genome_ann_path,"/",context))
    for (ann.loop.bedGragh in c("all", ann_res_files)) {
      suppressWarnings(try({DMRs_2_bedGragh(var1,var2,context,ann.loop.bedGragh)}, silent = T))
    }
    message("\tsaved all DMRs also as bedGragh files\n")
  }
  
  ### ### # finish main loop  ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  
  ###########################################################################
  
  ##### DMRs density - curcular plot #####
  tryCatch({
    setwd(exp_path)
    #setwd(ChrPlots_DMRs_path)
    source(paste0(scripts_dir,"circular_plot_DMRs.R"))
    DMRs_circular_plot(annotation.gr, TE_file, comparison_name)
    DMRs_circular_plot_legends()
    message("generated DMRs density plot for all contexts: done\n")
  }, error = function(cond) {
    cat("\n*\n",as.character(cond),"\n*\n")
    message("generated DMRs density plot for all contexts: fail\n")
  })
  
  ###########################################################################
  
  source(paste0(scripts_dir,"multiplot_ggplot2.R"))
  
  ##### GO analysis for annotated DMRs #####
  if (GO_analysis) {
    tryCatch({
      source(paste0(scripts_dir,"go_script.R"))
      GO_path = paste0(Methylome.At_res_dir,comparison_name,"/GO_analysis")
      dir.create(GO_path, showWarnings = FALSE)
      
      message("GO analysis for annotated DMRs...")
      run_GO(comparison_name, genome_ann_path, GO_path)
      message("GO analysis for annotated DMRs: done\n")
      
    }, error = function(cond) {
      cat("\n*\n",as.character(cond),"\n*\n")
      message("GO analysis for annotated DMRs: fail\n")
    })
  }
  
  ##### KEGG pathways for annotated DMRs ##### 
  if (KEGG_pathways) {
    tryCatch({
      source(paste0(scripts_dir,"kegg_script.R"))
      KEGG_path = paste0(Methylome.At_res_dir,comparison_name,"/KEGG_pathway")
      dir.create(KEGG_path, showWarnings = FALSE)
      
      message("KEGG pathways for annotated DMRs...")
      run_KEGG(comparison_name, genome_ann_path, KEGG_path)
      message("KEGG pathways for annotated DMRs: done\n")    
      
    }, error = function(cond) {
      cat("\n*\n",as.character(cond),"\n*\n")
      message("KEGG pathways for annotated DMRs: fail\n")
    })
  }
  
  ###########################################################################
  
  setwd(Methylome.At_res_dir)
  message(paste0("**\t",var2," vs ",var1,": done\n"))
  
  ###########################################################################
  
  # date and time of the end
  message(paste0("\n**\t",
                 paste(format(Sys.time(),"%d"),format(Sys.time(),"%m"),format(Sys.time(),"%Y"), sep = "-"),
                 " ",format(Sys.time(),"%H:%M")))
  end_time <- Sys.time()
  message(paste0("**\ttime: ", round(difftime(end_time,start_time, units = "hours"),2)," hours\n"))
}
