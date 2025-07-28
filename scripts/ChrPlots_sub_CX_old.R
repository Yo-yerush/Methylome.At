# ############
# ctrl_name <- var1
# trnt_name <- var2
# ctrl_pool <- meth_var1
# trnt_pool <- meth_var2
# num_cores <- 30
# chr_names <- unique(as.character(seqnames(annotation.gr)))
# ###########

run_ChrPlots_sub_CX <- function(ctrl_name, trnt_name, ctrl_pool, trnt_pool, num_cores, chr_names) {
    context_list <- list(
        CG = c("CGA", "CGT", "CGC", "CGG"),
        CHG = c("CAG", "CTG", "CCG"),
        CHH = c("CAA", "CAT", "CAC", "CTA", "CTT", "CTC", "CCA", "CCT", "CCC")
    )

    max_cores <- 9 * length(chr_names)
    if (num_cores >= max_cores) {
       num_cores_1 <- 9
       num_cores_2 <- length(chr_names)
    } else if (num_cores >= 18) {
       num_cores_1 <- 9
       num_cores_2 <- 2
    } else if (num_cores >= 9) {
       num_cores_1 <- 9
       num_cores_2 <- 1
    } else if (num_cores >= 3) {
       num_cores_1 <- 3
       num_cores_2 <- 1
    } else {
       num_cores_1 <- 1
       num_cores_2 <- 1
    }

    for (context_name in names(context_list)) {
        cat(context_name, "sub-contexts ... \n")

        subContext_list_ctrl <- GRangesList()
        subContext_list_trnt <- GRangesList()
        subContext_list_delta <- GRangesList()

        # parallel processing of subContexts
        subContext_par <- mclapply(context_list[[context_name]], function(sub_context) {
            profile_ctrl <- GRanges()
            profile_trnt <- GRanges()

            # parallel processing of chromosomes within each sub-context
            chr_results <- mclapply(as.character(ctrl_pool@seqnames@values), function(chr.n) {
                chr.vec <- ctrl_pool[ctrl_pool@seqnames == chr.n]
                regions <- GRanges(seqnames = Rle(chr.n), ranges = IRanges(chr.vec@ranges@start[1], chr.vec@ranges@start[length(chr.vec@ranges@start)]))

                ctrl_profile <- MethylationProfile(ctrl_pool, regions, windowSize = 500000, context = sub_context)
                trnt_profile <- MethylationProfile(trnt_pool, regions, windowSize = 500000, context = sub_context)

                return(list(ctrl = ctrl_profile, trnt = trnt_profile))
            }, mc.cores = num_cores)

            # combine all chromosomes
            for (chr_result in chr_results) {
                profile_ctrl <- c(profile_ctrl, chr_result$ctrl)
                profile_trnt <- c(profile_trnt, chr_result$trnt)
            }

            # delta
            profile_delta <- profile_ctrl
            profile_delta$Proportion <- (profile_trnt$Proportion - profile_ctrl$Proportion)

            return(list(
                ctrl = profile_ctrl,
                trnt = profile_trnt,
                delta = profile_delta
            ))
        }, mc.cores = num_cores_1)

        for (isc in seq_along(context_list[[context_name]])) {
            sub_context <- context_list[[context_name]][isc]
            subContext_list_ctrl[[sub_context]] <- subContext_par[[isc]]$ctrl
            subContext_list_trnt[[sub_context]] <- subContext_par[[isc]]$trnt
            subContext_list_delta[[sub_context]] <- subContext_par[[isc]]$delta
        }

        ################################################################
        ############# the plot #############

        max_proportions_fun <- function(x) {
            xx <- sapply(x, function(gr) {
                max(gr$Proportion)
            }) %>% max()
            ifelse(xx < 0, 0, xx)
        }
        min_proportions_fun <- function(x) {
            xx <- sapply(x, function(gr) {
                min(gr$Proportion)
            }) %>% min()
            ifelse(xx > 0, 0, xx)
        }

        max_proportions_ctrl <- max_proportions_fun(subContext_list_ctrl)
        max_proportions_trnt <- max_proportions_fun(subContext_list_trnt)
        max_proportions_delta <- max_proportions_fun(subContext_list_delta)

        min_proportions_ctrl <- min_proportions_fun(subContext_list_ctrl)
        min_proportions_trnt <- min_proportions_fun(subContext_list_trnt)
        min_proportions_delta <- min_proportions_fun(subContext_list_delta)

        max_comb <- max(c(max_proportions_ctrl, max_proportions_trnt))

        dir.create(paste0("sub", context_name), showWarnings = F)

        subC_plot(ctrl_name, subContext_list_ctrl, max_comb, 0, context_name, context_list)
        subC_plot(trnt_name, subContext_list_trnt, max_comb, 0, context_name, context_list)
        subC_plot("delta", subContext_list_delta, max_proportions_delta, min_proportions_delta, context_name, context_list)
    }
}

###########################################################################


subC_plot <- function(subC_var, subC_list, subC_max, subC_min, subC_CNTX, CNTX_list) {
    if (subC_var == "delta") {
        y_lab_fun <- paste(subC_CNTX, " methylation (Δ)")
    } else {
        y_lab_fun <- paste(subC_CNTX, " methylation")
    }

    svg(paste0("sub", subC_CNTX, "/", subC_var, "_sub", subC_CNTX, "_ChrPlot.svg"),
        width = 8, height = 2, family = "serif"
    )

    par(mar = c(1, 4, 2, 0))
    par(fig = c(0, 2, 0, 10) / 10)
    plot(runif(10), runif(10),
        xlim = c(0, 0.01),
        ylim = c(subC_min, subC_max), axes = FALSE, type = "n", ylab = y_lab_fun, xlab = ""
    )
    axis(2, c(subC_min, subC_max),
        lty = 1,
        labels = c(round(subC_min, 4), round(subC_max, 4))
    )
    par(new = T)

    i <- 1
    u <- (10 - i) / 6
    for (chr_number in 1:length(chr_names)) {
        par(mar = c(1, 0, 2, 0))
        par(fig = c(i, i + u, 0, 10) / 10)
        chromosome_plot(subC_list, # subContext_list_ctrl, subContext_list_trnt,
            chr_number, chr_names[chr_number], subC_CNTX, subC_var,
            difference = F, y.new.scale = T, y_cntx_max = subC_max, y_cntx_min = subC_min
        )
        # chromosome_plot(subC_list, chr_names[chr_number], subC_max, (subC_max + subC_nub) / 2, subC_min, "red", subC_CNTX)
        par(new = T)
        i <- i + u
    }

    ### legend
    par(mar = c(1, 0, 1, 0))
    par(fig = c(8, 10, 0, 10) / 10)
    legend("top",
        legend = CNTX_list[[subC_CNTX]],
        lty = 0, col = brewer.pal(n = length(CNTX_list[[subC_CNTX]]), name = "Set1"),
        pch = 15, bty = "n", cex = 0.75
    )

    dev.off()
}

########################################

chromosome_plot <- function(delta_profile, # var1_profile, var2_profile,
                            chr.n, chr.name, cntx, trnt,
                            difference = F, y.new.scale = F, y_cntx_max = NULL, y_cntx_min = NULL) {
    delta_profile <- delta_profile[seqnames(delta_profile) %in% chr.name]
    if (chr.n == 1) {
        ylab <- "methylation"
    } else {
        ylab <- ""
    }

    col <- brewer.pal(n = length(delta_profile), name = "Set1")
    pch <- rep(26, length(delta_profile))
    lty <- rep(1, length(delta_profile))
    lwd <- rep(2, length(delta_profile))

    if (y.new.scale) {
        ymax <- y_cntx_max
        ymin <- y_cntx_min
    } else {
        ymax <- 1
        ymin <- 0
    }

    pos <- (start(delta_profile[[1]]) + end(delta_profile[[1]])) / 2
    # lines for treatment var
    plot(pos, delta_profile[[1]]$Proportion,
        type = "o",
        ylim = c(ymin, ymax), xlab = "genomic coordinate", ylab = ylab,
        col = col[1], pch = pch[1], lty = lty[1],
        yaxt = "n", xaxt = "n",
        main = ""
    )


    # lines for the rest of treatment vars
    for (i in 2:length(delta_profile)) {
        lines(pos, delta_profile[[i]]$Proportion,
            type = "o", col = col[i], lty = lty[i], pch = pch[i]
        )
    }

    # line at 0 (yaxis)
    lines(pos, rep(0, length(delta_profile[[1]])),
        type = "o", col = "gray30", lty = 2, pch = 26
    )

    mtext(paste0("Chr ", chr.n), side = 1, line = 0, cex = 1) # adj = 0,

    # Adding text at the top left corner
    if (chr.n == 1) {
        text(
            x = par("usr")[1] + 0.125 * diff(par("usr")[1:2]), # Move a small fraction of the plot width to the right,
            y = par("usr")[4],
            labels = trnt,
            pos = 1, # Position text to the right of the specified coordinates
            # offset = 1,  # Add some space between the text and the top left corner
            cex = 0.9,
            font = 2,
            adj = c(0, 1)
        )
    }
}

########################################

####### 'DMRcaller' functions for manual editing
.sumReadsM <- function(methylationData) {
    return(sum(methylationData$readsM))
}

.sumReadsN <- function(methylationData) {
    return(sum(methylationData$readsN))
}

.analyseReadsInsideRegionsOneSample <- function(methylationData, regions) {
    overlaps <- findOverlaps(methylationData, regions, ignore.strand = TRUE)
    methylationDataContextList <- S4Vectors::splitAsList(methylationData[queryHits(overlaps)], subjectHits(overlaps))
    regionsIndexes <- as.integer(names(methylationDataContextList))

    regions$sumReadsM <- rep(0, times = length(regions))
    regions$sumReadsN <- rep(0, times = length(regions))
    regions$Proportion <- rep(0, times = length(regions))


    regions$sumReadsM[regionsIndexes] <- sapply(methylationDataContextList, .sumReadsM)
    regions$sumReadsN[regionsIndexes] <- sapply(methylationDataContextList, .sumReadsN)

    regions$Proportion[regionsIndexes] <- regions$sumReadsM[regionsIndexes] / regions$sumReadsN[regionsIndexes]
    return(regions)
}

####### yo modified 'DMRcaller' scripts
MethylationProfile <- function(methylationData, region, windowSize, context) {
    seqname <- seqnames(region)
    minPos <- start(region)
    maxPos <- end(region)
    hits <- findOverlaps(methylationData, region)
    localMethylationData <- methylationData[queryHits(hits)]
    rm(methylationData)

    if (context == "CG" | context == "CHG" | context == "CHH") {
        contextMethylationData <- localMethylationData[localMethylationData$context %in%
            context]
    } else {
        contextMethylationData <- localMethylationData[localMethylationData$trinucleotide_context %in%
            context]
    }

    rm(localMethylationData)

    seqs <- seq(minPos, maxPos - windowSize, windowSize)
    ranges <- GRanges(seqname, IRanges(seqs, seqs + windowSize -
        1))


    ranges <- .analyseReadsInsideRegionsOneSample(
        contextMethylationData,
        ranges
    )
    ranges$context <- paste(context, collapse = "_")
    return(ranges)
}
