Genes_features_metaPlot <- function(methylationPool_var1, methylationPool_var2, var1, var2, annotations_file, n.random, minReadsC, binSize, n.cores) {
  
  new_path.f = "Gene_features"
  dir.create(new_path.f, showWarnings = F)
  setwd(new_path.f)
  
  # Define regions as per your provided code
  Genes = annotations_file[which(annotations_file$type == "gene" & annotations_file$gene_model_type == "protein_coding")]
  Promoters = promoters(Genes, upstream=2000, downstream=0, use.names=TRUE)
  Promoters$type = "promoter"
  
  CDS = annotations_file[which(annotations_file$type == "CDS" & annotations_file$gene_model_type == "protein_coding")]
  fiveUTR = annotations_file[which(annotations_file$type == "five_prime_UTR" & annotations_file$gene_model_type == "protein_coding")]
  introns = annotations_file[which(annotations_file$type == "intron" & annotations_file$gene_model_type == "protein_coding")]
  threeUTR = annotations_file[which(annotations_file$type == "three_prime_UTR" & annotations_file$gene_model_type == "protein_coding")]
  
  # genes to analyze
  if (n.random[1] != "all") {
    if (length(n.random) == 1) {
    rndm_genes <- function(x) {return(sample(1:length(x), as.numeric(n.random)))} # random genes
    } else {
      rndm_genes <- function(x) {return(which(x$gene_id %in% n.random))} # list of TAIRs
    }
    
    Promoters = Promoters[rndm_genes(Promoters)]
    fiveUTR = fiveUTR[rndm_genes(fiveUTR)]
    CDS = CDS[rndm_genes(CDS)]
    introns = introns[rndm_genes(introns)]
    threeUTR = threeUTR[rndm_genes(threeUTR)]
  }
  
  regions_list = list("Promoters" = Promoters,
                      "fiveUTR" = fiveUTR,
                      "CDS" = CDS,
                      "introns" = introns,
                      "threeUTR" = threeUTR)
  
  region_names = names(regions_list)
  contexts = c("CG", "CHG", "CHH")
  
  # Function to process methylation data for each region and context
  genes_metaPlot_fun <- function(methylationData, regions_list, n.cores.f = n.cores) {
    
    # Filter positions below minReadsC
    methylationData = methylationData[which(methylationData$readsN >= minReadsC)]
    methylationData$Proportion = methylationData$readsM / methylationData$readsN
    
    # Filter methylationData by context before processing features
    methylationData_contexts = list()
    for (cntx in contexts) {
      methylationData_contexts[[cntx]] = methylationData[which(methylationData$context == cntx)]
    }
    
    
    result_list = list()
    
    for (region_name in region_names) {
      ann.obj = regions_list[[region_name]]
      
      #cat("\nProcessing region:", region_name, "with", length(ann.obj), "features\n")
      
      gene_2_bins_run <- function(feature.num) {
        tryCatch({
          region = ann.obj[feature.num][,0]
          seqname = seqnames(region)
          minPos = start(region)
          maxPos = end(region)
          
          context_results = list()
          
          # Create bins once
          windowSize = width(region) / binSize
          seqs = seq(minPos, maxPos, length.out = binSize+1)
          ranges = GRanges(seqname, IRanges(seqs[-length(seqs)], seqs[-1] - 1))
          
          for (cntx in contexts) {
            lMd_context = methylationData_contexts[[cntx]] %>% 
              subset(seqnames == as.character(seqname) & start >= minPos & end <= maxPos)
            
            if (length(lMd_context) >= 15) { # at least 15 sites...
              
              ranges$Proportion = rep(NA, length(ranges))
              
              # VECTORIZED APPROACH - Single findOverlaps call
              if (length(lMd_context) > 0) {
                hits.l = findOverlaps(lMd_context, ranges)
                
                if (length(hits.l) > 0) {
                  # Create data frame for aggregation
                  hit_df = data.frame(
                    subject_idx = subjectHits(hits.l),
                    proportion = lMd_context[queryHits(hits.l)]$Proportion
                  )
                  
                  # Calculate mean for each bin
                  bin_means = aggregate(proportion ~ subject_idx, hit_df, mean, na.rm = TRUE)
                  ranges$Proportion[bin_means$subject_idx] = bin_means$proportion
                }
              }
              
              # Handle reverse strand if necessary
              if (as.character(strand(region)) == "-") {
                ranges = rev(ranges)
              }
              
              context_results[[cntx]] = ranges
            } else {
              context_results[[cntx]] = NULL
            }
          }
          
          # print percentage every 100 genes
          if (gene.num %% 100 == 0) {
            cat("\rProcessed", length(ann.obj),  "features in", region_name, paste0("region [", round((feature.num/length(ann.obj))*100, 0), "%]  "))
          }

          return(context_results)
        }, error = function(cond) {
          return(NULL)
        })
      }

      cat("\n")
      cat("\rProcessed", length(ann.obj),  "features in", region_name, "region [0%]  ")
      #cat("\n\nPreparing", length(ann.obj), "bins for region:", region_name, "\n")
      results = mclapply(1:length(ann.obj), gene_2_bins_run, mc.cores = n.cores.f)
      results = results[!sapply(results, is.null)]
      cat("\n")

      # Average across features for each bin and context
      region_result = list()
      
      for (cntx in contexts) {
        gr_obj = GRanges(rep(region_name, binSize), IRanges(1:binSize,1:binSize))
        gr_obj$Proportion = rep(NA,binSize)
        
        # For each bin, average the Proportion values
        for (bin_num in 1:binSize) {
          bin_values = sapply(results, function(x) x[[cntx]]$Proportion[bin_num])
          gr_obj$Proportion[bin_num] = mean(unlist(bin_values), na.rm=TRUE)
        }
        
        region_result[[cntx]] = gr_obj
      }
      
      result_list[[region_name]] = region_result
      
    }
    
    return(result_list)
  }
  
  # Process methylation data for var1 and var2
  var1_metaPlot = genes_metaPlot_fun(methylationPool_var1, regions_list, n.cores.f = n.cores)
  cat(paste0("Finished calculating 'Genes_metaPlot' values for ", var1, "\n"))
  var2_metaPlot = genes_metaPlot_fun(methylationPool_var2, regions_list, n.cores.f = n.cores)
  cat(paste0("Finished calculating 'Genes_metaPlot' values for ", var2, "\n"))
  
  # Save the data
  dir.create("features_metaPlot_tables", showWarnings = F)
  for (region_name in region_names) {
    for (cntx in c("CG", "CHG", "CHH")) {
      try({write.csv(as.data.frame(var1_metaPlot[[region_name]][[cntx]]), paste0("features_metaPlot_tables/", var1, ".", cntx, ".", region_name, ".features.csv"), row.names = FALSE)})
      try({write.csv(as.data.frame(var2_metaPlot[[region_name]][[cntx]]), paste0("features_metaPlot_tables/", var2, ".", cntx, ".", region_name, ".features.csv"), row.names = FALSE)})
    }
  }
  
  # Plotting the data
  # bind data frames for metaPlot
  pos.end = binSize * length(region_names)
  
  #binsPlot <- function(v1.cntx.stream,v2.cntx.stream,var1,var2,cntx.m) {
  
  for (cntx in contexts) {
    
    # bind data frames by feature
    var1_proportions = c()
    var2_proportions = c()
    for (region_name in region_names) {
      var1_proportions = c(var1_proportions, var1_metaPlot[[region_name]][[cntx]]$Proportion)
      var2_proportions = c(var2_proportions, var2_metaPlot[[region_name]][[cntx]]$Proportion)
    }
    v.cntx = rbind(data.frame(pos = 1:pos.end, Proportion = var1_proportions, V = "V1"),
                   data.frame(pos = 1:pos.end, Proportion = var2_proportions, V = "V2"))
    
    # plot configuration
    min_value = min(v.cntx$Proportion)
    max_value = max(v.cntx$Proportion)
    q1_value = min_value+((max_value-min_value)/3)
    q2_value = min_value+((max_value-min_value)/3)+((max_value-min_value)/3)
    #middle_value = round(mean(c(min_value, max_value)), 2)
    q_value = min_value+((max_value-min_value)/4)
    
    legend_labels = c(paste0(" ",var1),paste0(" ",var2))
    
    region_names_plot = c("Promoter","5'UTR","CDS","Intron","3'UTR")
    main_title = "Protein Coding Genes"
    br.max = binSize * length(region_names)
    breaks_x = seq(0,br.max,by=binSize)
    breaks_vline = breaks_x[-c(1, length(breaks_x))]
    breaks_labels = breaks_x[1:length(region_names)] + binSize/2
    #breaks_and_labels <- list(breaks = seq(0,100,by=binSize), labels = region_names)
    
    plot_out = ggplot(data = v.cntx, aes(x = pos, y = Proportion, color = V, group = V)) +
      geom_vline(xintercept = breaks_vline, colour = "gray", linetype = "solid", linewidth = 0.6) +
      geom_line(linewidth = 0.65) +#, aes(linetype = V)) +
      scale_color_manual(values = c("V1" = "gray50", "V2" = "#bf6828")) +
      theme_classic() +
      labs(title = main_title, x = "", y = paste0(cntx," Methylation")) +
      theme(legend.position = "none",
            panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
            plot.title = element_text(hjust = 0.5),
            axis.text.y = element_text(size = 8),
            axis.text.x = element_text(size = 8, face = "bold")
      ) +
      scale_x_continuous(breaks = breaks_labels, labels = region_names_plot, minor_breaks = breaks_x, expand = expansion(add = c(0, 0))) +
      scale_y_continuous(breaks = c(min_value, q1_value, q2_value, max_value),
                         labels = c(round(min_value, 3), round(q1_value, 3),
                                    round(q2_value, 3), round(max_value, 3))) +
      annotate("text",
               x = max(breaks_vline),
               y = c(max_value*0.95, max_value*0.85),
               label = legend_labels,
               hjust = 0, vjust = 0.75, size = 3.25, 
               color = c("gray40","#bf6828"), fontface = "bold")
    
    
    svg(paste0("Genes_features_",cntx,"_metaPlot_",var2,"_vs_",var1,".svg"), width = 3.5, height = 2, family = "serif")
    #par(mar = c(2,2,1,2))
    print(plot_out)
    dev.off()
  } 
  
  message(paste("Processed average metaPlot for", length(n.random), "Protein Coding Genes\n"))
}
