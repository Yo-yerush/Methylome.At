DMRs_circular_plot <- function(ann.file, TE_4_dens, comparison_name) {
chr_amount = length(seqnames(ann.file)@values)
  
genes_type = ann.file[which(ann.file$type == "gene")]

  cntx_file <- function(context) {
    ############# read DMRs file
    dmrs_file = read.csv(paste0("DMRs_",context,"_",comparison_name,".csv"))
    dmrs_file = dmrs_file[,c("seqnames","start","end","log2FC")]
    return(dmrs_file)
  }
  
  CG_file = cntx_file("CG")
  CHG_file = cntx_file("CHG")
  CHH_file = cntx_file("CHH")
  
  #####################################
  ############# the plot #############
  img_device(paste0("DMRs_Density_",comparison_name), w = 3.25, h = 3.25, family = "serif")
  
  circos.par(start.degree = 90)
  circos.genomicInitialize(as.data.frame(ann.file)[,1:3], sector.names = paste0("Chr ",seq(chr_amount)), axis.labels.cex = 0.325, labels.cex = 1.35)
  
  circos.genomicDensity(list(CG_file[CG_file$log2FC > 0, 1:3],
                             CG_file[CG_file$log2FC < 0, 1:3]),
                        bg.col = "#fafcff", bg.border = NA, count_by = "number",
                        col = c("#FF000080","#304ed180"), border = T, track.height = 0.165, track.margin=c(0, 0))
  
  circos.genomicDensity(list(CHG_file[CHG_file$log2FC > 0, 1:3],
                             CHG_file[CHG_file$log2FC < 0, 1:3]),
                        bg.col = "#fafcff", bg.border = NA, count_by = "number",
                        col = c("#FF000080","#304ed180"), border = T, track.height = 0.165, track.margin=c(0, 0))
  
  circos.genomicDensity(list(CHH_file[CHH_file$log2FC > 0, 1:3],
                             CHH_file[CHH_file$log2FC < 0, 1:3]),
                        bg.col = "#fafcff", bg.border = NA, count_by = "number",
                        col = c("#FF000080","#304ed180"), border = T, track.height = 0.165, track.margin=c(0, 0))
  
  circos.genomicDensity(list(as.data.frame(genes_type)[1:3],
                             as.data.frame(TE_4_dens)[1:3]),
                        bg.col = "#fafcff", bg.border = NA, count_by = "number",
                        col =c("gray80","#fcba0320"), border = T, track.height = 0.165, track.margin=c(0, 0))
  
  circos.clear()
  dev.off()
}
############################################################

# legends
DMRs_circular_plot_legends <- function() {
  
  rndm_dis <- function() {
    counts_norm_hyper = hist(rnorm(3000, mean=0, sd=1), plot = F, breaks = 30)$counts
    counts_norm_hypo = hist(rnorm(3000, mean=0, sd=1), plot = F, breaks = 30)$counts
    rndm_norm_hyper = c((counts_norm_hyper - min(counts_norm_hyper)) / (max(counts_norm_hyper) - min(counts_norm_hyper)), rep(0,10))
    rndm_norm_hypo = c(rep(0,10), (counts_norm_hypo - min(counts_norm_hypo)) / (max(counts_norm_hypo) - min(counts_norm_hypo)))
    #rndm_norm_TE = c((counts_norm_hyper - min(counts_norm_hyper)) / (max(counts_norm_hyper) - min(counts_norm_hyper)), (counts_norm_hypo - min(counts_norm_hypo)) / (max(counts_norm_hypo) - min(counts_norm_hypo)))
    #rndm_norm_TE[rndm_norm_TE < 0.05] = rnorm(length(rndm_norm_TE[rndm_norm_TE < 0.05]), mean = 0.05, sd = 0.02)
    #rndm_norm_TE = rndm_norm_TE[!rndm_norm_TE == 0]
    return(list(hyper = rndm_norm_hyper, hypo = rndm_norm_hypo))
  }
  hyper.1 = rndm_dis()$hyper
  hyper.2 = rndm_dis()$hyper
  hyper.3 = rndm_dis()$hyper
  hypo.1 = rndm_dis()$hypo
  hypo.2 = rndm_dis()$hypo
  hypo.3 = rndm_dis()$hypo
  
  rndm_norm_TE = c(0, rndm_dis()$hypo[-c(1:10)], 0)
  rndm_norm_genes = (1-rndm_norm_TE)/1.8
  
  rndm = runif(12, min = 0, max = 1)
  legend_titles_circular = c("DMRs","  Group sig. genes","Methylation values","Expression values","  Overlapping TEs","Genome annotation")
  tracks_total_size = 16
  legend_tracks_pos = as.character((tracks_total_size/4)+1)
  text_tracks_pos = as.character((tracks_total_size/4))
  
  colfunc <- colorRampPalette(c("red", "white", "blue"))
  colfunc_vec_1 = colfunc(16)[-c((round(16/2)-1):(16/2), ((16/2)+1):((16/2)+1+1))]
  colfunc_vec_2 = colfunc(100)[-c((round(100/2)-5):(100/2)-1, ((100/2)+1):((100/2)+1+5))]
  
  
  img_device("legends", w = 2.5, h = 3.35, family = "serif")
  par(mar = c(0, 0, 0, 0) )
  circos.par("track.height" = 0.6, "canvas.xlim" = c(-0.2, 0.3), "canvas.ylim" = c(-0.25, 1), "gap.degree" = 0, "clock.wise" = FALSE)
  circos.initialize(factors = as.character(1:tracks_total_size), xlim = c(0, 1)) 
  
  circos.trackPlotRegion(factors = 1:tracks_total_size, ylim = c(0, 1), bg.border = NA, track.height = 0.175) 
  circos.updatePlotRegion(sector.index = legend_tracks_pos, bg.border = NA , bg.col = "#fafcfc", bg.lwd=0.2)
  circos.lines(seq(0,1, 1/(length(hyper.1)-1)), hyper.1, area = T, border = T, col = "#FF000080")
  circos.lines(seq(0,1, 1/(length(hypo.1)-1)), hypo.1, area = T, border = T, col = "#304ed180")
  circos.lines(seq(0,1, 1/(length(hyper.1)-1)), hyper.1, border = T)
  circos.text(text_tracks_pos, x = 0.1, y = 0.5, labels = "CG context", facing = "downward", ad = c(0, 0.5))
  
  circos.trackPlotRegion(factors = 1:tracks_total_size, ylim = c(0, 1), bg.border = NA, track.height = 0.175) 
  circos.updatePlotRegion(sector.index = legend_tracks_pos, bg.border = NA , bg.col = "#fafcfc", bg.lwd=0.2)
  circos.lines(seq(0,1, 1/(length(hyper.2)-1)), hyper.2, area = T, border = T, col = "#FF000080")
  circos.lines(seq(0,1, 1/(length(hypo.2)-1)), hypo.2, area = T, border = T, col = "#304ed180")
  circos.lines(seq(0,1, 1/(length(hyper.2)-1)), hyper.2, border = T)
  circos.text(text_tracks_pos, x = 0.14, y = 0.5, labels = "CHG context", facing = "downward", adj = c(0, 0.5))
  
  circos.trackPlotRegion(factors = 1:tracks_total_size, ylim = c(0, 1), bg.border = NA, track.height = 0.175) 
  circos.updatePlotRegion(sector.index = legend_tracks_pos, bg.border = NA , bg.col = "#fafcfc", bg.lwd=0.2)
  circos.lines(seq(0,1, 1/(length(hyper.3)-1)), hyper.3, area = T, border = T, col = "#FF000080")
  circos.lines(seq(0,1, 1/(length(hypo.3)-1)), hypo.3, area = T, border = T, col = "#304ed180")
  circos.lines(seq(0,1, 1/(length(hyper.3)-1)), hyper.3, border = T)
  circos.text(text_tracks_pos, x = 0.22, y = 0.5, labels = "CHH context", facing = "downward", adj = c(0, 0.5))
  
  circos.trackPlotRegion(factors = 1:tracks_total_size, ylim = c(0, 1), bg.border = NA, track.height = 0.15) 
  circos.updatePlotRegion(sector.index = legend_tracks_pos, bg.border = NA , bg.col = "#fafcfc", bg.lwd=0.2)
  circos.lines(seq(0,1, 1/(length(rndm_norm_genes)-1)), rndm_norm_genes, area = T, border = T, col = "gray80")
  circos.lines(seq(0,1, 1/(length(rndm_norm_TE)-1)), rndm_norm_TE, area = T, border = T, col = "#fcba0330")
  circos.text(text_tracks_pos, x = 0.37, y = 0.5, labels = "", facing = "downward", adj = c(0, 0.5))
  
  circos.clear()
  
  # genes/TEs
  legend(0,0.425,
         legend=c("Genes","TEs"),
         fill = c("#b2b2b2","#fcba0360"),
         bty="n")
  
  # DMRs
  legend(-0.12,0.125,
         legend=c(substitute(paste(bold("Hyper"),"-DMRs")),
                  substitute(paste(bold("Hypo"),"-DMRs")),
                  "Shared"),
         fill = c("#FF000095","#304ed195","#5e0d3d99"),
         bty="n")

  dev.off()
}
