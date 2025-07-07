Genes_metaPlot <- function(methylationPool_var1,methylationPool_var2,var1,var2,annotations_file,n.random,minReadsC,n.cores,is_TE=F) {
  
  if (is_TE) {
    new_path.f = "TEs"
    coding.Genes = annotations_file
  } else {
    new_path.f = "Genes"
    coding.Genes.0 <- annotations_file[which(annotations_file$type == "gene")]
    coding.Genes <- coding.Genes.0[which(coding.Genes.0$gene_model_type == "protein_coding")]
  }
  
  if (n.random[1] != "all") {
    if (length(n.random) == 1) {
      rndm_genes = sample(1:length(coding.Genes), as.numeric(n.random)) # random genes
      coding.Genes = coding.Genes[rndm_genes]
    } else {
      coding.Genes = coding.Genes[which(coding.Genes$gene_id %in% n.random)] # list of TAIRs
    }
  }
  
  dir.create(new_path.f, showWarnings = F)
  setwd(new_path.f)
  
  
  # make windowSize ranges with the average value
  genes_metaPlot_fun <- function(methylationData, ann.obj, n.cores.f = n.cores) {
    
    # filter CX position below 6 reads in total
    methylationData = methylationData[which(methylationData$readsN >= minReadsC)]
    
    methylationData$Proportion = methylationData$readsM/methylationData$readsN
    
    gene_2_bins_run <- function(gene.num) {
      tryCatch({
        region = ann.obj[gene.num][,0]
        seqname = seqnames(region)
        minPos = start(region)
        maxPos = end(region)
        
        # filter before the loop for faster running time
        localMethylationData = methylationData %>% 
          subset(seqnames == as.character(seqnames(region)) & start >= (minPos-2000) & end <= (maxPos+2000))
        
        gene_2_bins <- function(cntx.fun, stream_pos, lMd = localMethylationData) {

          # Pre-filter by context first (moved up)
          lMd = lMd[which(lMd$context == cntx.fun)]

          # Early return if no data
          if (length(lMd) == 0) {
            ranges = GRanges(seqname, IRanges(1:20, 1:20))
            ranges$Proportion = rep(NA, 20)
            return(ranges)
          }

          if (stream_pos == "up.stream") {
            lMd = lMd %>% filter(start >= (minPos-2000), end <= minPos)
            windowSize = 100
            seqs = seq((minPos-2000), minPos-windowSize, windowSize)

          } else if (stream_pos == "gene.body") {
            lMd = lMd %>% filter(start >= minPos, end <= maxPos)
            windowSize = width(region)/20
            seqs = seq(minPos, maxPos, windowSize)

          } else if (stream_pos == "down.stream") {
            lMd = lMd %>% filter(start >= maxPos, end <= (maxPos+2000))
            windowSize = 100
            seqs = seq(maxPos+1, (maxPos+2000-windowSize)+1, windowSize)
          }

          ranges = GRanges(seqname, IRanges(seqs, seqs+windowSize-1))
          ranges$Proportion = rep(NA, 20)

          # Vectorized findOverlaps approach
          if (length(lMd) > 0) {
            hits.l = findOverlaps(lMd, ranges)

            if (length(hits.l) > 0) {
              hit_df = data.frame(
                subject_idx = subjectHits(hits.l),
                proportion = lMd[queryHits(hits.l)]$Proportion
              )

              bin_means = aggregate(proportion ~ subject_idx, hit_df, mean, na.rm = TRUE)
              ranges$Proportion[bin_means$subject_idx] = bin_means$proportion
            }
          }
          return(ranges)
        }

        CG_list = GRangesList("up.stream"=GRanges(),"gene.body"=GRanges(), "down.stream"=GRanges())
        CHG_list = GRangesList("up.stream"=GRanges(),"gene.body"=GRanges(), "down.stream"=GRanges())
        CHH_list = GRangesList("up.stream"=GRanges(),"gene.body"=GRanges(), "down.stream"=GRanges())
        
        for (stream_pos in c("up.stream","gene.body","down.stream")) {
          CG_list[[stream_pos]] = gene_2_bins(cntx.fun = "CG", stream_pos = stream_pos)
          CHG_list[[stream_pos]] = gene_2_bins(cntx.fun = "CHG", stream_pos = stream_pos)
          CHH_list[[stream_pos]] = gene_2_bins(cntx.fun = "CHH", stream_pos = stream_pos)
        }
        
        # if strand is minus, change order of granges to reverse
        if (region@strand@values == "-") {
          strand_minus <- function(x) {
            x[[1]] = rev(x[[1]])
            x[[2]] = rev(x[[2]])
            x[[3]] = rev(x[[3]])
            x = rev(x)
            return(x)
          }
          CG_list = strand_minus(CG_list)
          CHG_list = strand_minus(CHG_list)
          CHH_list = strand_minus(CHH_list)
        }
        
        #cat("\rprepare 20bp proportional bins:",gene.num,"/",length(ann.obj),"      ")
        cat("\rcaculate average methylation in 20bp proportional bins over", length(ann.obj), new_path.f, "body and 2Kb up-/down-stream regions", paste0("[", round((gene.num/length(ann.obj)*100), 0), "%]  "))

        return(list(CG_list=CG_list,
                    CHG_list=CHG_list,
                    CHH_list=CHH_list))
        
      }, error = function(cond) {
        return(NULL)
        })
    }

    cat("\n")
    cat("\rcaculate average methylation in 20bp proportional bins over", length(ann.obj), new_path.f, "body and 2Kb up-/down-stream regions", paste0("[", 0, "%]  "))
    #cat("\nprepare 20bp proportional bins to",length(ann.obj),new_path.f,"\n")
    results = mclapply(1:length(ann.obj) , gene_2_bins_run, mc.cores = n.cores.f)
    results = results[!sapply(results, is.null)]
    cat("\n")
    ###################################################
    # prepare grange list for uptream, body and downstream
    gr_obj_u = GRanges(rep("up.stream",20), IRanges(1:20,1:20))
    gr_obj_g = GRanges(rep("gene.body",20), IRanges(1:20,1:20))
    gr_obj_d = GRanges(rep("down.stream",20), IRanges(1:20,1:20))
    gr_obj_u$Proportion = rep(NA,20)
    gr_obj_g$Proportion = rep(NA,20)
    gr_obj_d$Proportion = rep(NA,20)
    gr_list_CG = GRangesList("up.stream"=gr_obj_u,"gene.body"=gr_obj_g, "down.stream"=gr_obj_d)
    gr_list_CHG = GRangesList("up.stream"=gr_obj_u,"gene.body"=gr_obj_g, "down.stream"=gr_obj_d)
    gr_list_CHH = GRangesList("up.stream"=gr_obj_u,"gene.body"=gr_obj_g, "down.stream"=gr_obj_d)
    
    # Define a function to apply for each row_num
    process_row <- function(row_num,results,string_loc) {
      
        CG = mean(sapply(results, function(x) x$CG_list[[string_loc]][row_num]$Proportion), na.rm = TRUE)
        CHG = mean(sapply(results, function(x) x$CHG_list[[string_loc]][row_num]$Proportion), na.rm = TRUE)
        CHH = mean(sapply(results, function(x) x$CHH_list[[string_loc]][row_num]$Proportion), na.rm = TRUE)
        
      cat(paste0("\rprocessing average of ", gsub("\\."," ",string_loc), " bins: ",(row_num/20)*100,"%  "))
      return(list(CG = CG,
             CHG = CHG,
             CHH = CHH))
    }
    
    for (string_loc in 1:3) {
      # Parallel processing using mclapply
      results_parallel = mclapply(1:20, function(row_num) process_row(row_num, results,string_loc), mc.cores = ifelse(n.cores.f >= 20, 20, n.cores.f))
      
      # Assigning results back to your lists
      for (row_num_l in 1:20) {
        gr_list_CG[[string_loc]][row_num_l]$Proportion <- results_parallel[[row_num_l]]$CG
        gr_list_CHG[[string_loc]][row_num_l]$Proportion <- results_parallel[[row_num_l]]$CHG
        gr_list_CHH[[string_loc]][row_num_l]$Proportion <- results_parallel[[row_num_l]]$CHH
      }
      cat("\n")
    }
    return(list(gr_list_CG=gr_list_CG,
                gr_list_CHG=gr_list_CHG,
                gr_list_CHH=gr_list_CHH))
    
  }
  
  var1_metaPlot = genes_metaPlot_fun(methylationPool_var1, coding.Genes)
  cat(paste0("finish calculating 'Genes_metaPlot' values to ",var1,"\n\n"))
  var2_metaPlot = genes_metaPlot_fun(methylationPool_var2, coding.Genes)
  cat(paste0("finish calculating 'Genes_metaPlot' values to ",var2,"\n\n"))

  v1.CG = var1_metaPlot$gr_list_CG
  v1.CHG = var1_metaPlot$gr_list_CHG
  v1.CHH = var1_metaPlot$gr_list_CHH
  
  v2.CG = var2_metaPlot$gr_list_CG
  v2.CHG = var2_metaPlot$gr_list_CHG
  v2.CHH = var2_metaPlot$gr_list_CHH
  
  dir.create("metaPlot_tables", showWarnings = F)
  for (name in c("up.stream","gene.body","down.stream")) {
    try({write.csv(v1.CG[[name]], paste0("metaPlot_tables/",var1,".CG.",name,".csv"), row.names = FALSE)})
    try({write.csv(v1.CHG[[name]], paste0("metaPlot_tables/",var1,".CHG.",name,".csv"), row.names = FALSE)})
    try({write.csv(v1.CHH[[name]], paste0("metaPlot_tables/",var1,".CHH.",name,".csv"), row.names = FALSE)})
    try({write.csv(v2.CG[[name]], paste0("metaPlot_tables/",var2,".CG.",name,".csv"), row.names = FALSE)})
    try({write.csv(v2.CHG[[name]], paste0("metaPlot_tables/",var2,".CHG.",name,".csv"), row.names = FALSE)})
    try({write.csv(v2.CHH[[name]], paste0("metaPlot_tables/",var2,".CHH.",name,".csv"), row.names = FALSE)})
  }
  
  # bind data frames for metaPlot
  v1.CG.bind = data.frame(pos = 1:60, Proportion = c(v1.CG[[1]],v1.CG[[2]],v1.CG[[3]])$Proportion)
  v1.CHG.bind = data.frame(pos = 1:60, Proportion = c(v1.CHG[[1]],v1.CHG[[2]],v1.CHG[[3]])$Proportion)
  v1.CHH.bind = data.frame(pos = 1:60, Proportion = c(v1.CHH[[1]],v1.CHH[[2]],v1.CHH[[3]])$Proportion)
  
  v2.CG.bind = data.frame(pos = 1:60, Proportion = c(v2.CG[[1]],v2.CG[[2]],v2.CG[[3]])$Proportion)
  v2.CHG.bind = data.frame(pos = 1:60, Proportion = c(v2.CHG[[1]],v2.CHG[[2]],v2.CHG[[3]])$Proportion)
  v2.CHH.bind = data.frame(pos = 1:60, Proportion = c(v2.CHH[[1]],v2.CHH[[2]],v2.CHH[[3]])$Proportion)
  
  binsPlot <- function(v1.cntx.stream,v2.cntx.stream,var1,var2,cntx.m) {
    
    v1.cntx.stream$V = "V1"
    v2.cntx.stream$V = "V2"
    v.cntx.stream = rbind(v1.cntx.stream,v2.cntx.stream)
    
    min_value = min(v.cntx.stream$Proportion)
    max_value = max(v.cntx.stream$Proportion)
    q1_value = min_value+((max_value-min_value)/3)
    q2_value = min_value+((max_value-min_value)/3)+((max_value-min_value)/3)
    #middle_value = round(mean(c(min_value, max_value)), 2)
    q_value = min_value+((max_value-min_value)/4)
    
    if (is_TE) {
      legend_labels = c(var1,paste0("\n",var2))
    } else {
      legend_labels = c(paste0("\n",var1),paste0("\n\n",var2))
    }

    main_title = ifelse(is_TE, "TE", "Gene body")
    breaks_and_labels <- list(breaks = c(1, 20, 40, 60), labels = c("    -2kb", "TSS", "TTS", "+2kb    "))
    
    plot_out = ggplot(data = v.cntx.stream, aes(x = pos, y = Proportion, color = V, group = V)) +
      geom_vline(xintercept = c(20, 40), colour = "gray", linetype = "solid", linewidth = 0.5) +
      geom_line(linewidth = 0.65) +#, aes(linetype = V)) +
      #scale_color_manual(values = c("V1" = "#e37320", "V2" = "#0072B2")) +
      scale_color_manual(values = c("V1" = "gray50", "V2" = "#bf6828")) +
      #scale_linetype_manual(values = c("V1" = "dashed", "V2" = "solid")) +
      theme_classic() +
      labs(title = main_title,
           x = "",
           y = paste0(cntx.m," methylation")) +
      theme(legend.position = "none",
            panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
            plot.title = element_text(hjust = 0.5),
            axis.text.y = element_text(size = 8),
            axis.text.x = element_text(size = 9)
      ) +
      scale_x_continuous(breaks = breaks_and_labels$breaks, labels = breaks_and_labels$labels, expand = expansion(add = c(0,0))) +
      scale_y_continuous(breaks = c(min_value, q1_value, q2_value, max_value),
                         labels = c(round(min_value, 3), round(q1_value, 3),
                                    round(q2_value, 3), round(max_value, 3))) +
      annotate("text",
               x = 3,
               y = ifelse(is_TE,max_value,q1_value),
               label = legend_labels,
               hjust = 0, vjust = 0.75, size = 3.25, 
               color = c("gray40","#bf6828"), fontface = "bold")
    
    
    if (is_TE) {
      svg(paste0("TEs_",cntx.m,"_metaPlot_",var2,"_vs_",var1,".svg"), width = 1.88, height = 1.94, family = "serif")
    } else {
      svg(paste0("Genes_",cntx.m,"_metaPlot_",var2,"_vs_",var1,".svg"), width = 1.88, height = 1.94, family = "serif")
    }
    #par(mar = c(2,2,1,2))
    print(plot_out)
    dev.off()
  }  
  
  binsPlot(v1.CG.bind,v2.CG.bind,var1,var2,"CG")
  binsPlot(v1.CHG.bind,v2.CHG.bind,var1,var2,"CHG")
  binsPlot(v1.CHH.bind,v2.CHH.bind,var1,var2,"CHH")
  
  message(paste("process average metaPlot to",
                 length(coding.Genes),
                 ifelse(!is_TE,"Protein Coding Genes","Transposable Elements\n")))
}
