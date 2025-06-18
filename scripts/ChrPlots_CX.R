chromosome_plot <- function(profile_vars, chr.n, y_max, y_mid, y_min, col, cntx) {
    ### filter by chromosome
    profile_vars <- lapply(profile_vars, function(x) {
        x[x@seqnames == chr.n]
    })

    ### if the seqlevels are not the same, align them
    ### by taking the max length of each seqlevel
    maxLens <- Reduce(
        function(a, b) pmax(a, b, na.rm = TRUE),
        lapply(profile_vars, seqlengths)
    )
    common_si <- Seqinfo(names(maxLens), maxLens)
    profile_vars <- lapply(profile_vars, function(gr) {
        seqlevels(gr) <- seqlevels(common_si) # ensure same ordering
        seqinfo(gr) <- common_si # overwrite lengths
        gr
    })

    ### make it as gr list
    methylationProfiles <- GRangesList()
    for (i in seq_along(profile_vars)) {
        methylationProfiles[[i]] <- profile_vars[[i]]
    }

    pos <- (start(methylationProfiles[[1]]) + end(methylationProfiles[[1]])) / 2
    plot(pos, methylationProfiles[[1]]$Proportion,
        type = "o",
        ylim = c(y_min, y_max), xlab = "", ylab = "",
        col = col[1], pch = 26, lty = 1,
        yaxt = "n", xaxt = "n",
        main = "", frame.plot = FALSE
    )

    if (length(profile_vars) == 2) {
        lines(pos, methylationProfiles[[2]]$Proportion,
            type = "o", col = col[2], lty = 1, pch = 26
        )
    }

    if (cntx != "TE") {
        axis(1, at = c(min(pos), max(pos)), labels = FALSE, col = "gray20", tck = 0.025, pos = y_min) ###

    } else {
        # color bellow the lines
        polygon(c(pos, rev(pos)),
            c(methylationProfiles[[1]]$Proportion, rep(y_min, length(pos))),
            col = adjustcolor(col[1], alpha.f = 0.25), border = NA
        )
        axis(1, at = c(min(pos), max(pos)), labels = FALSE, col = "gray10", tck = 0.1, pos = 0) ###
    }

    if (cntx != "TE" & y_min != 0) {
        lines(pos, rep(0, length(pos)), type = "o", col = "gray50", lty = 2, pch = 26)
    }
}

###################################################################

ChrPlots_CX_all <- function(
    meth_var_list,
    meth_names,
    y_max_cg = 1,
    y_max_chg = 0.5,
    y_max_chh = 0.2,
    y_mid_cg = NULL,
    y_mid_chg = NULL,
    y_mid_chh = NULL,
    y_min_cg = 0,
    y_min_chg = 0,
    y_min_chh = 0,
    italic_legend_names = TRUE,
    ylab_suffix = NULL,
    y_title_cex = 1,
    TE_as_gr = NULL) {

    ### color palette
    if (length(meth_var_list) == 2) {
       col_vec <- c("gray40", "#bf6828")
    } else {
       col_vec <- "gray20"
    }

    ### y-mid edit
    y_mid_cg <- ifelse(is.null(y_mid_cg), (y_max_cg + y_min_cg) / 2, y_mid_cg)
    y_mid_chg <- ifelse(is.null(y_mid_chg), (y_max_chg + y_min_chg) / 2, y_mid_chg)
    y_mid_chh <- ifelse(is.null(y_mid_chh), (y_max_chh + y_min_chh) / 2, y_mid_chh)

    ### ylab suffix - in addition to 'CNTX methylation'
    ### add 'ylab_suffix=(delta)' to get 'CG methylation (delta)'
    if (!is.null(ylab_suffix)) {
        ylab_suffix <- paste0(" ", ylab_suffix)
    }

    ### change seqnames and value column name for plot
    meth_vars_trimmed <- lapply(meth_var_list, function(inner_list) {
        lapply(inner_list, function(m) {
            m <- renameSeqlevels(m, gsub("Chr", "", seqlevels(m)))
            names(mcols(m)) <- "Proportion"
            m
        })
    })

    ### normelize cheomosome panel to its size
    seq_len_vec <- as.numeric(seqlengths(meth_vars_trimmed[[1]][["chh"]]))
    chr_len <- seq_len_vec / max(seq_len_vec)

    ### Low resolution profiles plot ###
    ## column 1: y-axis strip
    ## columns 2-6: chromosomes 1-5
    ## column 7: legend

    if (is.null(TE_as_gr)) {
        lay <- cbind(
            matrix(1:24, nrow = 4, byrow = TRUE),
            rep(25, 4) # legend in column 7, spans all rows
        )

        layout(lay,
            widths  = c(0.375, chr_len, 0.5),
            heights = c(1, 1, 1, 0.25) # last row for Chrs
        )
    } else {
        lay <- cbind(
            matrix(1:30, nrow = 5, byrow = TRUE),
            rep(31, 5)
        )

        layout(lay,
            widths  = c(0.375, chr_len, 0.5),
            heights = c(1, 1, 1, 0.35, 0.30)
        )
    }


    for (cntx in c("CG", "CHG", "CHH", "TE")) {
        meth_vars_context <- lapply(meth_vars_trimmed, function(inner_list) inner_list[[tolower(cntx)]])

        if (cntx != "TE") {
            ## y-axis
            y_title <- paste0(cntx, " methylation", ylab_suffix)
            y_max_cntx <- ifelse(cntx == "CG", y_max_cg, ifelse(cntx == "CHG", y_max_chg, y_max_chh))
            y_mid_cntx <- ifelse(cntx == "CG", y_mid_cg, ifelse(cntx == "CHG", y_mid_chg, y_mid_chh))
            y_min_cntx <- ifelse(cntx == "CG", y_min_cg, ifelse(cntx == "CHG", y_min_chg, y_min_chh))

            par(mar = c(0, 4, 0, 0)) # c(1, 4, 2, 0))
            plot(NA, NA,
                xlim = c(0, 1), ylim = c(y_min_cntx, y_max_cntx),
                axes = FALSE, xlab = "",
                ylab = y_title,
                cex.lab = y_title_cex
            )
            axis(2,
                at = c(y_min_cntx, y_mid_cntx, y_max_cntx),
                labels = FALSE # ,  c(paste0("       ", y_min_cntx), y_mid_cntx, paste0(y_max_cntx, "      ")),
                # col = "gray35", cex.axis = 1
            )
            mtext(
                side = 2,
                text = y_min_cntx,
                at = y_min_cntx,
                adj = ifelse(y_min_cntx == 0, 0.5, 0), #
                line = 0.65,
                col = "gray25",
                cex = 0.55
            )
            mtext(
                side = 2,
                text = y_mid_cntx,
                at = y_mid_cntx,
                adj = 0.5, #
                line = 0.65,
                col = "gray25",
                cex = 0.55
            )
            mtext(
                side = 2,
                text = y_max_cntx,
                at = y_max_cntx,
                adj = ifelse(y_max_cntx == 0 | y_max_cntx == 1, 0.5, 0.9), #
                line = 0.65,
                col = "gray25",
                cex = 0.55
            )

            ## chromosome
            par(mar = c(0, 0, 0, 0)) # c(1, 0, 2, 0))
            for (chr in 1:5) {
                chromosome_plot(meth_vars_context, chr, y_max_cntx, y_mid_cntx, y_min_cntx, col_vec, cntx)
            }
        } else if (!is.null(TE_as_gr)) {
            te_plot_conf(TE_as_gr)
        }
    }

    ## x-axis - chromosome
    par(mar = c(1, 0, 0, 0))
    plot.new()
    for (chr in 1:5) {
        plot.new()
        mtext(
            side = 1,
            text = paste0("Chr ", chr),
            line = -0.5,
            at = 0.5,
            adj = 0.5,
            col = "gray25"
        )
    }

    ## legend
    par(mar = c(0, 0, 0, 0))
    plot.new()
    legend("top",
        legend = meth_names,
        text.font = ifelse(italic_legend_names, 3, 1),
        col = col_vec,
        lty = 1, bty = "n", cex = 1.2, lwd = 2
    )
}

###################################################################

te_plot_conf <- function(x) {
    te_gr <- x %>%
        circlize::genomicDensity(window.size = 0.5e6) %>%
        makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
        renameSeqlevels(gsub("Chr|chr|chromosome", "", seqlevels(.)))
    names(mcols(te_gr)) <- "Proportion"

    y_max_te <- 1 # max(te_gr$Proportion)
    y_mid_te <- 0.5 # y_max_te / 2
    y_min_te <- 0

    par(mar = c(0, 4, 0, 0))
    plot(NA, NA,
        xlim = c(0, 1), ylim = c(y_min_te, y_max_te),
        axes = FALSE, xlab = "",
        ylab = "TEs", cex.lab = 1.2
    )

    axis(2, at = c(y_min_te, y_max_te), labels = FALSE)
    mtext(
        side = 2, text = c("", ""), # c(y_min_te, y_mid_te, round(y_max_te, 2)),
        at = c(y_min_te, y_max_te) # ,
        # line = 0.65, col = "gray25", cex = 0.7, adj = c(0, .5, .9)
    )

    par(mar = c(0, 0, 0, 0))
    te_list <- list(te_gr)
    te_names <- "TE"

    for (chr in 1:5) {
        cat(".")
        chromosome_plot(te_list, chr,
            y_max_te, y_mid_te, y_min_te,
            col = "#55555590", cntx = "TE"
        ) ###
    }
}

###################################################################

windowSize_mcol <- function(x, mcol_name, windowSize = 1.5e5) {
    # chromosome lengths
    seqlens <- vapply(split(end(x), seqnames(x)), max, numeric(1))
    seqlengths(x) <- seqlens[seqlevels(x)]

    # windows by window size
    windows <- tileGenome(seqlens,
        tilewidth = windowSize,
        cut.last.tile.in.chrom = TRUE
    )

    # map to windows and calculate mean value
    hits <- findOverlaps(windows, x, ignore.strand = TRUE)
    mValue <- tapply(mcols(x)[[mcol_name]][subjectHits(hits)],
        queryHits(hits),
        mean,
        na.rm = TRUE
    )
    mcols(windows)$mean_value <- NA_real_
    mcols(windows)$mean_value[as.integer(names(mValue))] <- mValue

    cat(".")
    windows
}

###################################################################

var_sep <- function(a) {
    list(
        cg = a[a$context == "CG", ] %>% makeGRangesFromDataFrame(keep.extra.columns = T) %>% windowSize_mcol("Proportion"),
        chg = a[a$context == "CHG", ] %>% makeGRangesFromDataFrame(keep.extra.columns = T) %>% windowSize_mcol("Proportion"),
        chh = a[a$context == "CHH", ] %>% makeGRangesFromDataFrame(keep.extra.columns = T) %>% windowSize_mcol("Proportion")
    )
}

###################################################################

run_ChrPlots_CX <- function(comparison_name, meth_var1, meth_var2, var1, var2, TE.gr) {
    cat("\nCalculate methylated/unmethylated C's ratio...")
    meth_var1 <- as.data.frame(meth_var1) %>%
        mutate(Proportion = readsM / readsN) %>%
        select(seqnames, start, end, Proportion, context) %>%
        mutate(Proportion = ifelse(is.nan(Proportion), 0, Proportion)) # %>%
    # filter(!is.nan(Proportion)) ###
    cat(".")
    meth_var2 <- as.data.frame(meth_var2) %>%
        mutate(Proportion = readsM / readsN) %>%
        select(seqnames, start, end, Proportion, context) %>%
        mutate(Proportion = ifelse(is.nan(Proportion), 0, Proportion)) # %>%
    # filter(!is.nan(Proportion)) ###
    cat(".")
    meth_delta <- meth_var1 %>%
        mutate(Proportion = meth_var2$Proportion - meth_var1$Proportion) %>%
        select(seqnames, start, end, Proportion, context) %>%
        filter(!is.nan(Proportion))
    cat(" done\n")

    cat("\nChrPlots...")
    svg(paste0("ChrPlot_", comparison_name, ".svg"), width = 7, height = 4, family = "serif")
    try({
        ChrPlots_CX_all(
            meth_var_list = list(var_sep(meth_var1), var_sep(meth_var2)),
            meth_names = c(var1, var2),
            y_max_cg = 1,
            y_max_chg = 0.6,
            y_max_chh = 0.25,
            y_mid_cg = NULL,
            y_mid_chg = NULL,
            y_mid_chh = NULL,
            y_min_cg = 0,
            y_min_chg = 0,
            y_min_chh = 0,
            italic_legend_names = FALSE,
            ylab_suffix = NULL,
            y_title_cex = 1,
            TE_as_gr = TE.gr
        )
    })
    dev.off()
    cat(" done\n")

    cat("ChrPlots (difference)...")
    svg(paste0("ChrPlot_difference_", comparison_name, ".svg"), width = 7, height = 4, family = "serif")
    try({
        ChrPlots_CX_all(
            meth_var_list = list(var_sep(meth_delta)),
            meth_names = comparison_name,
            y_max_cg = 0.1,
            y_max_chg = 0.1,
            y_max_chh = 0.1,
            y_mid_cg = 0,
            y_mid_chg = 0,
            y_mid_chh = 0,
            y_min_cg = -0.1,
            y_min_chg = -0.1,
            y_min_chh = -0.1,
            italic_legend_names = FALSE,
            ylab_suffix = "(Δ)",
            y_title_cex = 1,
            TE_as_gr = TE.gr
        )
    })
    dev.off()
    cat(" done\n")
}