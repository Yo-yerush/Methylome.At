library(data.table)
library(dplyr)
library(GenomicRanges)
library(circlize)



jointed_gr <- methylationDataReplicates_joints
windowSize <- 1000

# read TFBS file
tfbs_data <- fread("/home/yoyerush/yo/methylome_pipeline/motifs_db_files/unibind/TAIR10_compressed_TFBSs.bed.gz")

colnames(tfbs_data) <- c("seqnames", "start", "end", "motif", "phase", "strand", "thickStart", "thickEnd", "itemRgb")

# extract moif column information
tfbs_data$name_parsed <- strsplit(as.character(tfbs_data$motif), "_")
tfbs_data$Experiment <- sapply(tfbs_data$name_parsed, `[`, 1)
tfbs_data$Condition <- sapply(tfbs_data$name_parsed, `[`, 2)
tfbs_data$TF <- sapply(tfbs_data$name_parsed, `[`, 3)
tfbs_data$MotifID <- sapply(tfbs_data$name_parsed, `[`, 4)

# group TF families
tf_family <- c(
    AG = "MADS-box", AP1 = "MADS-box", AP3 = "MADS-box", PI = "MADS-box",
    SEP3 = "MADS-box", SOC1 = "MADS-box", SVP = "MADS-box", AGL27 = "MADS-box",
    FLC = "MADS-box", LFY = "MADS-box",
    PIF1 = "bHLH", PIF3 = "bHLH", PIF4 = "bHLH", PIF5 = "bHLH",
    WRKY18 = "WRKY", WRKY33 = "WRKY", WRKY40 = "WRKY",
    ABF1 = "bZIP", ABF3 = "bZIP", ABF4 = "bZIP", GBF2 = "bZIP", GBF3 = "bZIP", HY5 = "bZIP",
    ARR1 = "Type-B ARR", ARR10 = "Type-B ARR", ARR14 = "Type-B ARR",
    `ATHB-5` = "HD-ZIP", `ATHB-6` = "HD-ZIP", `ATHB-7` = "HD-ZIP", HAT22 = "HD-ZIP", JKD = "HD-ZIP",
    TCP4 = "TCP",
    MYB3 = "MYB",
    E2FA = "E2F",
    REF6 = "Chromatin regulator"
)

tfbs_data$TF_family <- tf_family[tfbs_data$TF]


tfbs_gr <- tfbs_data %>%
    mutate(
        seqnames = gsub("chr", "Chr", tfbs_data$seqnames),
        hex = sapply(strsplit(as.character(tfbs_data$itemRgb), ","), function(x) {
            rgb(as.numeric(x[1]), as.numeric(x[2]), as.numeric(x[3]), maxColorValue = 255)
        })
    ) %>%
    filter(seqnames != "ChrMt") %>%
    filter(seqnames != "ChrPt") %>%
    dplyr::select(-c(motif, name_parsed, itemRgb, thickStart, thickEnd, phase)) %>%
    dplyr::relocate(strand, TF, TF_family, .after = end) %>%
    makeGRangesFromDataFrame(., keep.extra.columns = T)


# m1 <- findOverlaps(meth_var1, tfbs_gr)
# var1_meth_tfbs <- meth_var1[queryHits(m1)]
# mcols(var1_meth_tfbs) <- cbind.data.frame(
#     mcols(var1_meth_tfbs),
#     mcols(tfbs_gr[subjectHits(m1)])
# )

# m2 <- findOverlaps(meth_var2, tfbs_gr)
# var2_meth_tfbs <- meth_var2[queryHits(m2)]
# mcols(var2_meth_tfbs) <- cbind.data.frame(
#     mcols(var2_meth_tfbs),
#     mcols(tfbs_gr[subjectHits(m2)])
# )


m <- findOverlaps(jointed_gr, tfbs_gr)
joined_meth_tfbs <- jointed_gr[queryHits(m)]
mcols(joined_meth_tfbs) <- cbind.data.frame(
    mcols(joined_meth_tfbs),
    mcols(tfbs_gr[subjectHits(m)])
)



dmp_df <- dmp_def(joined_meth_tfbs)



chr_amount <- 5 # length(seqnames(ann.file)@values)

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

img_device("TFs_lines_density_superfamilies",
    w = 4.25, h = 4.25
)


par(mar = c(0, 0, 0, 0))

circos.par(gap.degree = c(rep(1, chr_amount - 1), 35), start.degree = 90, points.overflow.warning = FALSE)
circos.genomicInitialize(as.data.frame(ann.file)[, 1:3], sector.names = paste0("Chr ", 1:chr_amount), axis.labels.cex = 0.4, labels.cex = 1.25)
for (family.i in unique(as.character(tf_family))) {
    cat(family.i, "\n")

    suppressMessages({
        density_data <- genomicDensity(dmp_df[which(dmp_df$TF_family == family.i)], window.size = windowSize, count_by = "number")
        ylims <- range(density_data$value) * 1.435
        circos.genomicTrackPlotRegion(density_data, ylim = range(density_data$value), bg.border = NA, track.height = 0.065, track.margin = c(0, 0), panel.fun = function(region, value, ...) {
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
            colors <- ifelse(value >= 15, "#440154",
                ifelse(value >= 10, "#31688e",
                    ifelse(value >= 7, "#21918c",
                        ifelse(value >= 5, "#35b779",
                            ifelse(value >= 3, "#90d743",
                                "#d9d9d9"
                            )
                        )
                    )
                )
            )
            circos.genomicLines(region, value,
                col = colors,
                border = TRUE, lty = 1, lwd = 0.5, type = "h"
            )
        })
        ### y-axis labels
        circos.text("Chr1", x = 0, y = 0.5, labels = paste0(family.i, "  "), facing = "downward", cex = 0.6, adj = c(0.85, -0.15))
    })
}
circos.clear()

legend("topleft",
    legend = c("≥15", "≥10", "≥7", "≥5", "≥3", "<3"),
    fill = c("#440154", "#31688e", "#21918c", "#35b779", "#90d743", "#d9d9d9"),
    title = "DMR count",
    cex = 0.6,
    bty = "n"
)

dev.off()
