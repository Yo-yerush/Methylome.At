run_sum_deltaH_CX <- function(ctrl_name, trnt_name, ctrl_pool, trnt_pool, fdr = 0.95) {
    ctrl_pool <- as.data.frame(ctrl_pool) %>%
        # filter(readsN > 6) %>%
        select(seqnames, start, context, m_ctrl = readsM, n_ctrl = readsN)
    trnt_pool <- as.data.frame(trnt_pool) %>%
        # filter(readsN > 6) %>%
        select(seqnames, start, context, m_trnt = readsM, n_trnt = readsN)

    merged_pool <- inner_join(
        ctrl_pool, trnt_pool,
        by = c("seqnames", "start", "context")
    )

    # main loop
    cat("sum dH in 100bp windowsize\n")
    hd_list <- mclapply(c("CG", "CHG", "CHH"), function(cntx) {
        dH_bins(merged_pool, cntx)
    }, mc.cores = 3)
    cg_hd <- hd_list[[1]]
    chg_hd <- hd_list[[2]]
    chh_hd <- hd_list[[3]]

    # FDR cut-off
    q_cg <- quantile(cg_hd$sum_surprisal, fdr) %>% as.numeric()
    q_chg <- quantile(chg_hd$sum_surprisal, fdr) %>% as.numeric()
    q_chh <- quantile(chh_hd$sum_surprisal, fdr) %>% as.numeric()

    # filter by FDR threshold
    cg_filtered <- cg_hd %>% filter(sum_surprisal > q_cg)
    chg_filtered <- chg_hd %>% filter(sum_surprisal > q_chg)
    chh_filtered <- chh_hd %>% filter(sum_surprisal > q_chh)

    cat("\nsignificants dH:\n")
    print(knitr::kable(
        data.frame(
            context = c("CG", "CHG", "CHH"),
            total_regions = c(nrow(cg_hd), nrow(chg_hd), nrow(chh_hd)),
            "5_percent_top_regions" = c(nrow(cg_filtered), nrow(chg_filtered), nrow(chh_filtered)),
            sValue_threshold = paste("S >", c(round(q_cg, 1), round(q_chg, 1), round(q_chh, 1)))
        )
    ))

    # genome density plot
    dHRs_circular_plot(cg_filtered, chg_filtered, chh_filtered, annotation.gr, TE_file, paste0(trnt_name, "_vs_", ctrl_name))

    # manhattan plots
    mclapply(
        list(list(cg_hd, "CG", fdr), list(chg_hd, "CHG", fdr), list(chh_hd, "CHH", fdr)),
        function(x) {
            plots_test(x[[1]], x[[2]], x[[3]])
        },
        mc.cores = 3
    )

    # NAs # # write bigWig file
    # NAs # write_bw <- function(x, cntx) {
    # NAs #     x %>%
    # NAs #         select(seqnames, start, end, sum_surprisal) %>%
    # NAs #         na.omit() %>%
    # NAs #         makeGRangesFromDataFrame(., keep.extra.columns = TRUE) %>%
    # NAs #     rtracklayer::export.bw(. con = paste0("surprisal_", cntx, "_", trnt_name, "_vs_", ctrl_name, ".bw"))
    # NAs # }
    # NAs # write_bw(cg_filtered, "CG")
    # NAs # write_bw(chg_filtered, "CHG")
    # NAs # write_bw(chh_filtered, "CHH")
}

########################################################################

dH_bins <- function(joint, cntx, min_coverage = 6) {
    thresh <- c(CG = 0.4, CHG = 0.2, CHH = 0.1)[cntx]
    min_C_sites <- c(CG = 1, CHG = 2, CHH = 4)[cntx]
    min_coverage <- 6

    cat(paste0("[", cntx, "]\tmin coverage: ", min_coverage, "; min proportion value: ", thresh, "; min dignificant sites: ", min_C_sites, "\n"))
    ### per-site surprisal value [S-Value = −log P_Bin(c | n, π0)]
    eps <- 1e-6 # avoid log(0)
    joint <- joint %>%
        filter(context == cntx, n_ctrl >= min_coverage, n_trnt >= min_coverage) %>%
        mutate(
            delta = abs(m_trnt / n_trnt - m_ctrl / n_ctrl),
        ) %>%
        filter(delta > thresh) %>%
        mutate(
            pi0 = pmin(pmax(m_ctrl / n_ctrl, eps), 1 - eps), # ctrl expectation π0
            surprisal = -(lchoose(n_trnt, m_trnt) +
                m_trnt * log(pi0) +
                (n_trnt - m_trnt) * log(1 - pi0))
        )

    ### 100bp windowSize
    gr_sites <- GRanges(joint$seqnames, IRanges(joint$start, width = 1))
    seqlevels <- unique(joint$seqnames)
    seqlengths_per_chr <- joint %>%
        group_by(seqnames) %>%
        summarise(seqlength = max(start)) %>%
        arrange(match(seqnames, seqlevels))

    seqlengths(gr_sites) <- seqlengths_per_chr$seqlength
    names(seqlengths(gr_sites)) <- seqlevels
    tiles100 <- unlist(tileGenome(seqlengths(gr_sites), tilewidth = 100))

    ### sum surprisal within each window                  ###
    ov <- findOverlaps(gr_sites, tiles100, ignore.strand = TRUE)
    joint$tile_id <- subjectHits(ov) # map each site → window

    window_surprisal <- joint %>%
        group_by(tile_id) %>%
        summarise(
            sum_surprisal = sum(surprisal, na.rm = TRUE),
            pi_ctrl = sum(m_ctrl) / sum(n_ctrl),
            pi_trnt = sum(m_trnt) / sum(n_trnt),
            log2FC = log2(pi_trnt / pi_ctrl),
            n_sites = n()
        ) %>%
        ungroup()

    # final df
    out <- cbind(
        as.data.frame(tiles100)[window_surprisal$tile_id, 1:3],
        window_surprisal
    ) %>%
        filter(n_sites >= min_C_sites)
    out
}

########################################################################

dHRs_circular_plot <- function(CG_df, CHG_df, CHH_df, ann.file, TE_4_dens, comparison_name) {
    chr_amount <- length(seqnames(ann.file)@values)

    genes_type <- ann.file[which(ann.file$type == "gene")]

    ############# the plot #############
    svg(paste0("dHR_Density_", comparison_name, ".svg"), width = 3.25, height = 3.25, family = "serif")

    circos.par(start.degree = 90)
    circos.genomicInitialize(as.data.frame(ann.file)[, 1:3], sector.names = paste0("Chr ", seq(chr_amount)), axis.labels.cex = 0.325, labels.cex = 1.35)

    circos.genomicDensity(
        list(
            CG_df[CG_df$log2FC > 0, 1:3],
            CG_df[CG_df$log2FC < 0, 1:3]
        ),
        bg.col = "#fafcff", bg.border = NA, count_by = "number",
        col = c("#FF000080", "#304ed180"), border = T, track.height = 0.165, track.margin = c(0, 0)
    )

    circos.genomicDensity(
        list(
            CHG_df[CHG_df$log2FC > 0, 1:3],
            CHG_df[CHG_df$log2FC < 0, 1:3]
        ),
        bg.col = "#fafcff", bg.border = NA, count_by = "number",
        col = c("#FF000080", "#304ed180"), border = T, track.height = 0.165, track.margin = c(0, 0)
    )

    circos.genomicDensity(
        list(
            CHH_df[CHH_df$log2FC > 0, 1:3],
            CHH_df[CHH_df$log2FC < 0, 1:3]
        ),
        bg.col = "#fafcff", bg.border = NA, count_by = "number",
        col = c("#FF000080", "#304ed180"), border = T, track.height = 0.165, track.margin = c(0, 0)
    )

    circos.genomicDensity(
        list(
            as.data.frame(genes_type)[, 1:3],
            as.data.frame(TE_4_dens)[, 1:3]
        ),
        bg.col = "#fafcff", bg.border = NA, count_by = "number",
        col = c("gray80", "#fcba0320"), border = T, track.height = 0.165, track.margin = c(0, 0)
    )

    circos.clear()
    dev.off()
}

########################################################################

plots_test <- function(df, cntx, fdr = 0.95) {
    S_threshold <- as.numeric(quantile(df$sum_surprisal, fdr))

    png(paste0(cntx, "_test_plots_310725.png"), width = 3.5, height = 2.75, units = "in", res = 300, family = "serif")

    # add cumulative genome coordinate
    chr_lengths <- tapply(df$end, df$seqnames, max)
    chr_offset <- c(0, cumsum(chr_lengths[-length(chr_lengths)]))
    names(chr_offset) <- names(chr_lengths)

    df$coord <- df$start + chr_offset[df$seqnames]

    print(
        ggplot(df, aes(coord, sum_surprisal,
            colour = as.factor(as.integer(seqnames) %% 2)
        )) +
            geom_point(size = 0.025) +
            geom_hline(yintercept = S_threshold, color = "#bf6828", linetype = "dashed", size = 0.7) +
            scale_colour_manual(values = c("grey40", "#7ca182"), guide = "none") +
            scale_x_continuous(
                breaks = chr_offset,
                labels = paste0("Chr", seq_along(chr_offset))
            ) +
            labs(
                x = "Chromosome",
                y = "Window surprisal"
            ) +
            theme_bw()
    )
    dev.off()
}
