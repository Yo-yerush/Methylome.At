Genes_type_plot <- function(context) {
  
  #dir.create("TEs_addiotionnal_results", showWarnings = F)
  #setwd("TEs_addiotionnal_results")
  
  tmp_path = "C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_pipeline/DMRcaller/Methylome.At/check_ann_files_res_CG"
  Genes_path = paste0(tmp_path,"/Genes_",context,"_genom_annotations.csv")
  TEG_path = paste0(tmp_path,"/TEG_",context,"_genom_annotations.csv")
  pseudogene_path = paste0(tmp_path,"/pseudogene_",context,"_genom_annotations.csv")

  #DMRsReplicates_Genes_path = paste0("../",context,"/Genes_",context,"_genom_annotations.csv")
  #DMRsReplicates_TEG_path = paste0("../",context,"/Genes_",context,"_genom_annotations.csv")
  #DMRsReplicates_pseudogene_path = paste0("../",context,"/Genes_",context,"_genom_annotations.csv")
    
  read.dmr.file <- function(DMRs_ann_path, is.genes=F) {
    if (file.exists(DMRs_ann_path)) {
      x = read.csv(DMRs_ann_path)
      if (is.genes) {x = x[x$gene_model_type == "protein_coding",]}
      x = x %>% 
        distinct(gene_id, .keep_all = T) %>% 
        select(gene_id, regionType)
      return(x)
    }
  }
    Genes = read.dmr.file(Genes_path,T)
    TEG = read.dmr.file(TEG_path)
    pseudogene = read.dmr.file(pseudogene_path)
    
    ######################### TE vs Genes plot
    GvsT_df = data.frame(x=c("Coding-Genes","TEGs", "Pseudogenes"),
                         y=c(nrow(Genes),nrow(TEG),nrow(pseudogene)))
    GvsT_plot = ggplot(GvsT_df,aes(x=x, y=y)) + 
      geom_bar(stat = 'identity', position= "stack")+
      geom_col()+ 
      xlab("")+
      ylab("Number of DMRs")+
      ggtitle(context)+
      theme_classic() + theme(axis.text.x =  element_text(face="bold",size=11,angle = 45))
    svg(file = paste0(context,"_TE.vs.ProteinCodingGenes.svg"), width = 1.6, height = 2.1, family = "serif")
    plot(GvsT_plot)
    dev.off()
  }
}