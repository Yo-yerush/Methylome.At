gr_2_bigWig <- function(gr, out_file_name = "gr_2_bigWig.bw", out_format = "bigWig") { # scaling=TRUE
    
    # keep genome location and methylation proportion
    gr_out <- gr
    mcols(gr_out) <- NULL
    gr_out$score <- round(gr$log2FC, 2)

    # get seqlength for arabidopsis
    chrInfo <- as.list(org.At.tairCHRLENGTHS)
    seqlengths <- unlist(chrInfo)[1:5]
    names(seqlengths) <- paste0("Chr", 1:5)

    # save as BigWig file
    seqlengths(gr_out) <- seqlengths
    export(gr_out, pout_file_name, format = out_format)
}
