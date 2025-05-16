DMRs_2_bigWig <- function(var1, var2, context, ann.res) {
    if (ann.res == "all") {
        DMRsReplicates <- read.csv(paste0("DMRs_", context, "_", var2, "_vs_", var1, ".csv"))
        file_name <- paste(context, "all", var2, "vs", var1, sep = "_")
    } else {
        DMRsReplicates <- read.csv(paste0("genome_annotation/", context, "/", ann.res))
        ann.name <- sub(paste0("_", context, ".*$"), "", ann.res)
        file_name <- paste(context, ann.name, var2, "vs", var1, sep = "_")
    }

    # keep genome location and methylation proportion
    DMRsReplicates_wig <- DMRsReplicates[, c("seqnames", "start", "end")]
    DMRsReplicates_wig$score <- round(log2(DMRsReplicates$proportion2 / DMRsReplicates$proportion1), 2)
    DMRsReplicates_wig <- makeGRangesFromDataFrame(DMRsReplicates_wig, keep.extra.columns = T)

    # get seqlength for arabidopsis
    chrInfo <- as.list(org.At.tairCHRLENGTHS)
    seqlengths <- unlist(chrInfo)[1:5]
    names(seqlengths) <- paste0("Chr", 1:5)

    # save as BigWig file
    seqlengths(DMRsReplicates_wig) <- seqlengths
    export(DMRsReplicates_wig, paste0("DMRs_bigWig/", file_name, ".bw"), format = "bigWig")
}
