DMRs_ann <- function(annotation_vec, DMRsReplicates, context, description_file) {
  
  ann_count_df = data.frame(type = c(names(annotation_vec)),
                            x = rep(NA, length(annotation_vec))) # counts DMRs for log file
  names(ann_count_df)[2] = context
  
  dir.create(context)
  setwd(context)
  
  ######################## main annotation loop
  for (i in 1:length(annotation_vec)) {
    # find overlaps for DMRs with annotation file
    m <- findOverlaps(DMRsReplicates, annotation_vec[[i]])
    DMRs_annotation <- DMRsReplicates[queryHits(m)]
    mcols(DMRs_annotation) <- cbind.data.frame(
      mcols(DMRs_annotation),
      mcols(annotation_vec[[i]][subjectHits(m)]))
    
    if (length(DMRs_annotation) != 0) {
      
      type_name = names(annotation_vec[i])
      
      DMRs_annotation_df = sort(DMRs_annotation) %>% as.data.frame()
      
      # edit to merge with 'description_file' (if not TE annotations)
      if (type_name != "Transposable_Elements") {
        DMRs_annotation_df = apply(DMRs_annotation_df,2,as.character)
        
        # keep the columns from 'DMRsReplicates' and 'gene_id' from 
        col_keep = c("gene_id","type", names(as.data.frame(DMRsReplicates)))
        DMRs_annotation_df = DMRs_annotation_df[,which(colnames(DMRs_annotation_df) %in% col_keep)]
        
        # merge with araprot11 and uniprot databases description files
        DMRs_annotation_df = merge.data.frame(DMRs_annotation_df,description_file, by = "gene_id", all.x = T)
      }
      
      ### edit columns positions
      end_columns = grep("sumReads|proportion|cytosinesCount", names(DMRs_annotation_df))
      DMRs_annotation_df = DMRs_annotation_df %>% 
        relocate(names(DMRs_annotation_df)[end_columns], .after = last_col()) 
      
      # make unique by both TAIR and position
      DMRs_annotation_df$tmp_pos = paste(DMRs_annotation_df$gene_id, 
                                         DMRs_annotation_df$seqnames,
                                         DMRs_annotation_df$start,
                                         DMRs_annotation_df$end , 
                                         sep = "_DEL_")
      DMRs_annotation_df = DMRs_annotation_df %>% distinct(tmp_pos, .keep_all = T) %>% select(-tmp_pos)
      ann_count_df[i,context] = nrow(DMRs_annotation_df) # counts DMRs for log file
      
      ###########################################
      ### annotate 'TEG' and 'pseudogene' if needed
      if (!contains_TEG_n_pseudogene) {
        
        if (names(annotation_vec[i]) == "Genes") {
          TEG = DMRs_annotation_df[DMRs_annotation_df$gene_model_type == "transposable_element_gene",]
          write.csv(TEG,
                    paste0("TEG_",context,"_genom_annotations.csv"), 
                    row.names = F, na = "")
          
          pseudogene = DMRs_annotation_df[DMRs_annotation_df$gene_model_type == "pseudogene",]
          write.csv(pseudogene,
                    paste0("pseudogene_",context,"_genom_annotations.csv"), 
                    row.names = F, na = "")
          
          message(paste0("\t\tdone: TEG (",nrow(TEG),")"))
        }
        
        # remove TEG if there are (in non-TEGs annotations)
        #if (type_name != "Transposable_Elements") {
        #  DMRs_annotation_df = DMRs_annotation_df[!DMRs_annotation_df$gene_model_type == "transposable_element_gene",]
        #}
      }
      
      ###########################################
      
      tryCatch(
        {write.csv(DMRs_annotation_df,paste0(type_name,"_",context,"_genom_annotations.csv"),
                   row.names = F, na = "")
          message(paste0("\t\tdone: ",ann_count_df[i,1]," (",ann_count_df[i,context],")"))
        }, error = function(cond) {
          message(paste0("\t\tfail: ",ann_count_df[i,1]))}
      )
    } else {ann_count_df[i,context] = 0}
  }
  
  setwd("../")
}