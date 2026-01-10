DMRs_circular_plot <- function(ann.file, TE_4_dens, comparison_name) {
  chr_amount <- length(seqnames(ann.file)@values)

  genes_type <- ann.file[which(ann.file$type == "gene")]

  cntx_file <- function(context) {
    ############# read DMRs file
    dmrs_file <- read.csv(paste0("DMRs_", context, "_", comparison_name, ".csv"))
    dmrs_file <- dmrs_file[, c("seqnames", "start", "end", "log2FC")]
    return(dmrs_file)
  }

  CG_file <- cntx_file("CG")
  CHG_file <- cntx_file("CHG")
  CHH_file <- cntx_file("CHH")

  #####################################
  ############# the plot #############
  img_device(paste0("DMRs_Density_", comparison_name), w = 3, h = 3)

  par(mar = c(0, 0, 0, 0))

  circos.par(gap.degree = c(rep(1, chr_amount - 1), 25), start.degree = 90, points.overflow.warning = FALSE)
  circos.genomicInitialize(as.data.frame(ann.file)[, 1:3], sector.names = paste0("Chr ", seq(chr_amount)), axis.labels.cex = 0.325, labels.cex = 1.25)

  circos.genomicDensity(
    list(
      CG_file[CG_file$log2FC > 0, 1:3],
      CG_file[CG_file$log2FC < 0, 1:3]
    ),
    bg.col = "#fafcff", bg.border = NA, count_by = "number",
    col = c("#FF000080", "#304ed180"), border = T, track.height = 0.165, track.margin = c(0, 0)
  )
  circos.text("Chr1", x = 0, y = 1, labels = "CG", facing = "downward", cex = 0.6, adj = c(1, -0.4))

  circos.genomicDensity(
    list(
      CHG_file[CHG_file$log2FC > 0, 1:3],
      CHG_file[CHG_file$log2FC < 0, 1:3]
    ),
    bg.col = "#fafcff", bg.border = NA, count_by = "number",
    col = c("#FF000080", "#304ed180"), border = T, track.height = 0.165, track.margin = c(0, 0)
  )
  circos.text("Chr1", x = 0, y = 1, labels = "CHG", facing = "downward", cex = 0.6, adj = c(1, -0.4))

  circos.genomicDensity(
    list(
      CHH_file[CHH_file$log2FC > 0, 1:3],
      CHH_file[CHH_file$log2FC < 0, 1:3]
    ),
    bg.col = "#fafcff", bg.border = NA, count_by = "number",
    col = c("#FF000080", "#304ed180"), border = T, track.height = 0.165, track.margin = c(0, 0)
  )
  circos.text("Chr1", x = 0, y = 1, labels = "CHH", facing = "downward", cex = 0.6, adj = c(1, -0.4))

  circos.genomicDensity(
    list(
      as.data.frame(genes_type)[1:3],
      as.data.frame(TE_4_dens)[1:3]
    ),
    bg.col = "#fafcff", bg.border = NA, count_by = "number",
    col = c("gray80", "#fcba0320"), border = T, track.height = 0.165, track.margin = c(0, 0)
  )

  circos.clear()

  legend("topleft",
    legend = c("Hyper-DMRs", "Hypo-DMRs", "Overlay", "Genes", "TEs"),
    fill = c("#FF000095", "#304ed195", "#8208b695", "#9c9c9c", "#fcba0360"),
    cex = 0.5,
    bty = "n"
  )

  dev.off()
}
