run_GO <- function(comparison_name, genome_ann_path, GO_path) {
  
  for (ann_loop in c("Genes","Promoters","CDS","Introns","fiveUTRs","threeUTRs")) {
    for (contx_loop in c("CG","CHG","CHH")) {
      for (gain_loss_loop in c("gain","loss")) {
        for (GO_type_loop in c("BP","MF","CC")) {
          tryCatch({
            suppressMessages(top.GO.fun(treatment = comparison_name,
                                        gain_OR_loss = gain_loss_loop,
                                        GO_Ontology_type = GO_type_loop,
                                        context = contx_loop,
                                        annotation = ann_loop,
                                        genome_ann_path = genome_ann_path,
                                        path_for_results = GO_path))}, 
            error = function(cond) {
              message(paste0("Error in 'top.GO.fun': ",paste(ann_loop,contx_loop,gain_loss_loop,GO_type_loop,sep = "-")))}
          )
        }
      }
      tryCatch(
        {
          GO_one_plot(comparison_name,contx_loop,ann_loop,GO_path)
          #message(paste0("\tdone: ",ann_loop," in ",contx_loop," context\t"))
        }, error = function(cond) {
          message(paste0("Error in GO plot: ",paste(ann_loop,contx_loop,gain_loss_loop,GO_type_loop,sep = "-")))}
      )
    }
    message(paste0("\tdone: ",ann_loop))
  }
}

########################################################

top.GO.fun = function(treatment,
                      gain_OR_loss, # "gain" OR "loss"
                      GO_Ontology_type, # BP, CC, MF
                      context, # "CG", "CHG", "CHH" or "all"
                      annotation, # "Genes" for example
                      pcutoff = 0.01,
                      GO_algorithm = "weight01", # classic, weight, weight01
                      genome_ann_path,
                      path_for_results,
                      n.nodes = NULL,
                      get_GO_term_genes = c(F, "GO_term"),
                      save.files = T) {
  
  ##########################################
  start_path = getwd()
  
  #####
  DMR_file = read.csv(paste0(genome_ann_path,"/",context,"/",annotation,"_",context,"_genom_annotations.csv"))
  DMR_file = DMR_file[DMR_file$regionType == gain_OR_loss,
                      c("gene_id","pValue")]
  
  
  tair_ids <- data.frame(gene_id = keys(org.At.tair.db))
  
  all_genes = merge.data.frame(DMR_file, tair_ids, by = "gene_id", all.y = T)
  all_genes$pValue[is.na(all_genes$pValue)] <- 0.999
  all_genes$pValue[all_genes$pValue == 0] <- 1e-300
  
  tmp = ifelse(all_genes$pValue < pcutoff, 1, 0)
  geneList = tmp
  names(geneList) = all_genes$gene_id
  
  ############################################################
  #### type Ontology  ####
  myGOdata <- new("topGOdata",
                  ontology = GO_Ontology_type,
                  allGenes = geneList,
                  geneSelectionFun = function(x)(x == 1),
                  #description = "Test",
                  annot = annFUN.org,
                  #nodeSize = 5,
                  mapping="org.At.tair.db") 
  
  sg <- sigGenes(myGOdata)
  str(sg)
  numSigGenes(myGOdata)
  
  resultFisher <- runTest(myGOdata, algorithm = GO_algorithm, statistic = "fisher")
  
  
  allRes <- GenTable(myGOdata, Fisher = resultFisher,
                     orderBy= "Fisher", ranksOf= "Fisher", topNodes=length(resultFisher@score))
  allRes$Fisher = as.numeric(allRes$Fisher)
  allRes$Term = gsub(",", ";", allRes$Term)
  
  allRes_view = allRes[allRes$Fisher <= pcutoff,]
  
  ### save
  new_path_1 = paste(path_for_results,context, sep = "/")
  new_path = paste(path_for_results,context,annotation, sep = "/")
  if (save.files == T) {
    if (dir.exists(new_path_1) == F) {
      dir.create(new_path_1)
    }
    if (dir.exists(new_path) == F) {
      dir.create(new_path)
    }
    write.csv(allRes_view, paste0(new_path,"/",GO_Ontology_type ,"." ,annotation,
                                  ".",context,".",gain_OR_loss,".",treatment,".topGO.csv"),
              quote = F, row.names = F)
  }
  
  
  if (is.null(n.nodes)) {
    if (nrow(allRes_view) > 10) {
      n.nodes = 10
    } else {
      n.nodes = nrow(allRes_view)
    }
  }
  
  setwd(new_path)
  printGraph(myGOdata, resultFisher, 
             firstSigNodes = n.nodes, fn.prefix= paste(GO_Ontology_type,annotation,context,gain_OR_loss,treatment, sep = "_"),
             useInfo="all", pdfSW= T)
  
  setwd(start_path)
}

########################################################

GO_one_plot <- function(treatment,
                        context,
                        annotation,
                        path_for__results,
                        breaks_gain = NULL, # c(2,4,6),
                        limits_gain = NULL, # c(0.5,6.5)
                        breaks_loss = NULL,
                        limits_loss = NULL
) {
  
  new_path = paste(path_for__results,context,annotation, sep = "/")
  
  bp_gain = read.csv(paste0(new_path,"/BP.",annotation,".",context,".gain.",treatment,".topGO.csv"))
  bp_loss = read.csv(paste0(new_path,"/BP.",annotation,".",context,".loss.",treatment,".topGO.csv"))
  cc_gain = read.csv(paste0(new_path,"/CC.",annotation,".",context,".gain.",treatment,".topGO.csv"))
  cc_loss = read.csv(paste0(new_path,"/CC.",annotation,".",context,".loss.",treatment,".topGO.csv"))
  mf_gain = read.csv(paste0(new_path,"/MF.",annotation,".",context,".gain.",treatment,".topGO.csv"))
  mf_loss = read.csv(paste0(new_path,"/MF.",annotation,".",context,".loss.",treatment,".topGO.csv"))
  
  if (length(bp_gain$GO.ID) != 0) {bp_gain$type = "BP"}
  if (length(bp_loss$GO.ID) != 0) {bp_loss$type = "BP"}
  if (length(cc_gain$GO.ID) != 0) {cc_gain$type = "CC"}
  if (length(cc_loss$GO.ID) != 0) {cc_loss$type = "CC"}
  if (length(mf_gain$GO.ID) != 0) {mf_gain$type = "MF"}
  if (length(mf_loss$GO.ID) != 0) {mf_loss$type = "MF"}
  
  gain_bind = rbind(bp_gain, cc_gain, mf_gain)
  loss_bind = rbind(bp_loss, cc_loss, mf_loss)
  
  gain_bind = gain_bind[!grepl("cellular_component|biological_process|molecular_function", gain_bind$Term),]
  loss_bind = loss_bind[!grepl("cellular_component|biological_process|molecular_function", loss_bind$Term),]
  
  gain_col = "#cf534c"
  loss_col = "#6397eb"
  
  bubble_gain = gain_bind %>% 
    ggplot(aes(Significant, reorder(Term,-Fisher), size = Annotated, color = Fisher)) + 
    scale_color_gradient("p.value", low = gain_col, high = "black") + #theme_classic() +
    labs(#title = paste0(type_name," - ",gain_loss,"regulated transcripts"),
      x = "Significant", y = "") + theme_bw() + 
    theme(
      #plot.title=element_text(hjust=0.5),
      legend.key.size = unit(0.25, 'cm'),
      legend.title = element_text(size=9.5),
      legend.position = "right",
      text = element_text(family = "serif")) + 
    geom_point() +
    facet_grid(rows = vars(type), scales = "free_y", space = "free_y") + 
    guides(color = guide_colorbar(order = 1, barheight = 4))
  
  if (is.null(breaks_gain) == F & is.null(limits_gain) == F) {
    bubble_gain = bubble_gain + scale_x_continuous(breaks=breaks_gain,
                                                   limits=limits_gain)
  }
  
  
  
  bubble_loss = loss_bind %>% 
    ggplot(aes(Significant, reorder(Term,-Fisher), size = Annotated, color = Fisher)) + 
    scale_color_gradient("p.value", low = loss_col, high = "black") + #theme_classic() +
    labs(#title = paste0(type_name," - ",gain_loss,"regulated transcripts"),
      x = "Significant", y = "") + theme_bw() + 
    theme(
      #plot.title=element_text(hjust=0.5),
      legend.key.size = unit(0.25, 'cm'),
      legend.title = element_text(size=9.5),
      legend.position = "right",
      text = element_text(family = "serif")) + 
    geom_point() +
    facet_grid(rows = vars(type), scales = "free_y", space = "free_y") + 
    guides(color = guide_colorbar(order = 1, barheight = 4))
  
  if (is.null(breaks_loss) == F & is.null(limits_loss) == F) {
    bubble_loss = bubble_loss + scale_x_continuous(breaks=breaks_loss,
                                                   limits=limits_loss)
  }
  
  
  Height = max(c(nrow(gain_bind),nrow(loss_bind)))/6.25
  if (Height < 3) {Height = 3}
  #print(paste(treatment,context,annotation, sep = "_"))
  
  svg(paste0(path_for__results,"/",treatment,"_",context,"_",annotation,"_GO.svg"),
      width = 9.90, height = Height, family = "serif")
  multiplot(bubble_gain, bubble_loss, cols=2)
  dev.off()
}
