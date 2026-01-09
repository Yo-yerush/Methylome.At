TEs_superfamily_circular_plot <- function(ann.file) {
    chr_amount <- length(seqnames(ann.file)@values)
    genes_type <- ann.file[which(ann.file$type == "gene")]

    ############# Heterochromatin positions #############
    heteroChr <- data.frame(
        Chr = paste0("Chr", c(1:chr_amount)),
        start = c(12500000, 1250000, 11000000, 1666667, 9444444),
        end = c(17500000, 7500000, 16250000, 7000000, 15000000)
    )

    ############# Centromere positions #############
    cenChr <- data.frame(
        Chr = paste0("Chr", c(1:chr_amount)),
        start = c(14476796, 3462971, 13780083, 3177188, 11207348),
        end = c(15081019, 3650512, 14388500, 3248799, 11555278)
    )

    #############

    cntx_file <- function(x) read.csv(paste0(x, "/DMRs_Transposable_Elements_", x, "_genom_annotations.csv"))

    all_cntx_gr <- rbind(
        cntx_file("CG"),
        cntx_file("CHG"),
        cntx_file("CHH")
    )

    ##########

    family_res_list <- split(all_cntx_gr[, 1:3], all_cntx_gr$Transposon_Super_Family)[order(sapply(split(all_cntx_gr[, 1:3], all_cntx_gr$Transposon_Super_Family), nrow), decreasing = TRUE)][1:6]
    names_2_keep <- names(family_res_list)
    names(family_res_list) <- gsub("LTR/|RC/|^DNA/|/L1$", "", names(family_res_list))

    ###############################################
    ############# circular lines plot #############
    img_device("TEs_addiotionnal_results/TEs_lines_density_superfamilies",
        w = 4.25, h = 4.25
    )
    circos.par(gap.degree = c(rep(4, 4), 35), start.degree = 90)
    circos.genomicInitialize(as.data.frame(ann.file)[, 1:3], sector.names = paste0("Chr ", 1:chr_amount), axis.labels.cex = 0.4, labels.cex = 1.35)
    for (family.i in 1:length(family_res_list)) {
        suppressMessages({
            density_data <- genomicDensity(family_res_list[[family.i]], window.size = 1e5, count_by = "number")
            ylims <- range(density_data$value) * 1.435
            circos.genomicTrackPlotRegion(density_data, ylim = range(density_data$value), bg.border = NA, track.height = 0.105, track.margin = c(0, 0), panel.fun = function(region, value, ...) {
                chr.n <- gsub("Chr", "", get.cell.meta.data("sector.index"))

                ### heterocromatin
                circos.rect(heteroChr[chr.n, 2], ylims[1], heteroChr[chr.n, 3], ylims[2], # xleft, ybottom, xright, ytop
                    col = "#fcba0320",
                    border = NA
                )
                ### centromere
                circos.rect(cenChr[chr.n, 2], ylims[1], cenChr[chr.n, 3], ylims[2], # xleft, ybottom, xright, ytop
                    col = "#fcba0360",
                    border = NA
                )

                ### density lines
                colors <- ifelse(value >= 3, "green4", "gray15")
                circos.genomicLines(region, value,
                    col = colors,
                    border = TRUE, lty = 1, lwd = 0.5, type = "h"
                )
            })
            ### y-axis labels
            circos.text("Chr1", x = 0, y = 0.5, labels = paste0(names(family_res_list)[family.i], "  "), facing = "downward", cex = 0.6, adj = c(1, 0))
        })
    }
    circos.clear()
    dev.off()





    ###############################################
    ############# circular dencsity plot ##########

    gain_gr <- all_cntx_gr[all_cntx_gr$regionType == "gain", ]
    loss_gr <- all_cntx_gr[all_cntx_gr$regionType == "loss", ]

    family_list_gain <- split(gain_gr[, 1:3], gain_gr$Transposon_Super_Family)[order(sapply(split(gain_gr[, 1:3], gain_gr$Transposon_Super_Family), nrow), decreasing = TRUE)][names_2_keep]
    family_list_loss <- split(loss_gr[, 1:3], loss_gr$Transposon_Super_Family)[order(sapply(split(loss_gr[, 1:3], loss_gr$Transposon_Super_Family), nrow), decreasing = TRUE)][names_2_keep]
    names(family_list_gain) <- gsub("LTR/|RC/|^DNA/|/L1$", "", names(family_list_gain))
    names(family_list_loss) <- gsub("LTR/|RC/|^DNA/|/L1$", "", names(family_list_loss))

    ##########

    img_device("TEs_addiotionnal_results/TEs_direction_density_superfamilies",
        w = 4.25, h = 4.25
    )
    circos.par(gap.degree = c(rep(4, 4), 35), start.degree = 90)
    circos.genomicInitialize(as.data.frame(ann.file)[, 1:3], sector.names = paste0("Chr ", 1:chr_amount), axis.labels.cex = 0.4, labels.cex = 1.35)
    for (family.i in 1:length(family_list_gain)) {
        suppressMessages({
            density_gain <- genomicDensity(family_list_gain[[family.i]], window.size = 1e6, count_by = "number")
            density_loss <- genomicDensity(family_list_loss[[family.i]], window.size = 1e6, count_by = "number")
            y_max <- max(c(density_gain$value, density_loss$value))
            ylims <- c(0, y_max) * 1.435

            circos.genomicTrackPlotRegion(cbind(density_gain, density_loss$value),
                ylim = c(0, y_max),
                bg.border = NA, track.height = 0.105, track.margin = c(0, 0),
                panel.fun = function(region, value, ...) {
                    chr.n <- gsub("Chr", "", get.cell.meta.data("sector.index"))

                    ### heterocromatin
                    circos.rect(heteroChr[chr.n, 2], ylims[1], heteroChr[chr.n, 3], ylims[2], # xleft, ybottom, xright, ytop
                        col = "#fcba0320",
                        border = NA
                    )
                    ### centromere
                    circos.rect(cenChr[chr.n, 2], ylims[1], cenChr[chr.n, 3], ylims[2], # xleft, ybottom, xright, ytop
                        col = "#fcba0360",
                        border = NA
                    )

                    ### density lines
                    circos.genomicLines(region, value[, 1], col = "#FF000080", border = TRUE, lty = 1, lwd = 0.5, type = "l", area = T)
                    circos.genomicLines(region, value[, 2], col = "#304ed180", border = TRUE, lty = 1, lwd = 0.5, type = "l", area = T)
                }
            )
            ### y-axis labels
            circos.text("Chr1", x = 0, y = 0.5, labels = paste0(names(family_list_gain)[family.i], "  "), facing = "downward", cex = 0.6, adj = c(1, 0))
        })
    }
    circos.clear()
    dev.off()
}
