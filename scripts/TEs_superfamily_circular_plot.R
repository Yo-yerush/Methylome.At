TEs_superfamily_circular_plot <- function(ann.file, TE_4_dens) {
    ############# Heterochromatin positions #############
    heteroChr <- data.frame(
        Chr = paste0("Chr", c(1:5)),
        start = c(12500000, 1250000, 11000000, 1666667, 9444444),
        end = c(17500000, 7500000, 16250000, 7000000, 15000000)
    )

    ############# Centromere positions #############
    cenChr <- data.frame(
        Chr = paste0("Chr", c(1:5)),
        start = c(14476796, 3462971, 13780083, 3177188, 11207348),
        end = c(15081019, 3650512, 14388500, 3248799, 11555278)
    )

    #############

    cntx_file <- function(x) read.csv(paste0(x, "/Transposable_Elements_", x, "_genom_annotations.csv"))

    all_cntx_gr <- rbind(
        cntx_file("CG"),
        cntx_file("CHG"),
        cntx_file("CHH")
    ) %>%
        select(seqnames, start, end, Transposon_Super_Family)

    ##########
    family_res_list <- list(
        SINE = all_cntx_gr[grep("^SINE|RathE", all_cntx_gr$Transposon_Super_Family), 1:3],
        LINE = all_cntx_gr[grep("^LINE", all_cntx_gr$Transposon_Super_Family), 1:3],
        Copia = all_cntx_gr[grep("LTR/Copia", all_cntx_gr$Transposon_Super_Family), 1:3],
        Helitron = all_cntx_gr[grep("RC/Helitron", all_cntx_gr$Transposon_Super_Family), 1:3],
        TIR = all_cntx_gr[grep("^DNA", all_cntx_gr$Transposon_Super_Family), 1:3],
        Gypsy = all_cntx_gr[grep("LTR/Gypsy", all_cntx_gr$Transposon_Super_Family), 1:3]
    )

    ### Calculate density
    density_TE <- genomicDensity(as.data.frame(TE_4_dens)[, 1:3], window.size = 1e6, count_by = "number")
    # new scaling
    density_TE$value <- density_TE$value / max(density_TE$value)


    #####################################
    ############# the plot #############
    img_device(
        paste0("genome_annotation/TEs_addiotionnal_results/TEs_superfamilies_circular_plot"),
        w = 4.25, h = 4.25
    )
    circos.par(gap.degree = c(rep(4, 4), 35), start.degree = 90)
    circos.genomicInitialize(as.data.frame(ann.file)[, 1:3], sector.names = paste0("Chr ", 1:5), axis.labels.cex = 0.4, labels.cex = 1.35)
    for (family.i in 1:length(family_res_list)) {
        density_data <- genomicDensity(family_res_list[[family.i]], window.size = 1e5, count_by = "number")
        ylims <- range(density_data$value) * 1.435
        circos.genomicTrackPlotRegion(density_data, ylim = range(density_data$value), bg.border = NA, track.height = 0.085, track.margin = c(0, 0), panel.fun = function(region, value, ...) {
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
            colors <- ifelse(value > 5, "red3", "gray15")
            circos.genomicLines(region, value,
                col = colors,
                border = TRUE, lty = 1, lwd = 0.5, type = "h"
            )
        })
        ### y-axis labels
        circos.text("Chr1", x = 0, y = 0.5, labels = paste0(names(family_res_list)[family.i], "  "), facing = "downward", cex = 0.6, adj = c(1, 0))
    }
    ### TEs
    circos.genomicTrackPlotRegion(density_TE, bg.border = NA, track.height = 0.05, track.margin = c(0, 0), panel.fun = function(region, value, ...) {
        circos.genomicLines(region, value,
            col = "#fcba0320", border = TRUE, type = "l", area = T,
            track.index = get.cell.meta.data("track.index") - 1
        )
    })
    circos.clear()
    dev.off()
}