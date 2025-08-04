run_sum_deltaH_CX <- function(ctrl_name, trnt_name, ctrl_pool, trnt_pool, annotation.gr, TE.gr, num_cores, fdr = 0.95) {
    ctrl_pool <- as.data.frame(ctrl_pool) %>%
        select(seqnames, start, context, m_ctrl = readsM, n_ctrl = readsN)
    trnt_pool <- as.data.frame(trnt_pool) %>%
        select(seqnames, start, context, m_trnt = readsM, n_trnt = readsN)

    merged_pool <- inner_join(
        ctrl_pool, trnt_pool,
        by = c("seqnames", "start", "context")
    )

    seqlengths_merged <- merged_pool %>%
        group_by(seqnames) %>%
        summarise(seqlength = max(start)) %>%
        arrange(match(seqnames, unique(merged_pool$seqnames)))

    # main loop
    cat("sum dH in 100bp windowsize\n")
    hd_list <- mclapply(
        c("CG", "CHG", "CHH"), function(cntx) {
            dH_bins(merged_pool, cntx)
        },
        mc.cores = ifelse(num_cores >= 3, 3, num_cores)
    )
    cg_hd <- hd_list[[1]]
    chg_hd <- hd_list[[2]]
    chh_hd <- hd_list[[3]]

    # filter by FDR threshold
    q_cg <- quantile(cg_hd$sum_surprisal, fdr) %>% as.numeric()
    q_chg <- quantile(chg_hd$sum_surprisal, fdr) %>% as.numeric()
    q_chh <- quantile(chh_hd$sum_surprisal, fdr) %>% as.numeric()
    cg_filtered <- cg_hd %>% filter(sum_surprisal > q_cg)
    chg_filtered <- chg_hd %>% filter(sum_surprisal > q_chg)
    chh_filtered <- chh_hd %>% filter(sum_surprisal > q_chh)

    # print df (both to terminal and .log file)
    df_2_print <- data.frame(
        context = c("CG", "CHG", "CHH"),
        total = c(nrow(cg_hd), nrow(chg_hd), nrow(chh_hd)),
        top_5 = c(nrow(cg_filtered), nrow(chg_filtered), nrow(chh_filtered)),
        threshold = paste0("S>", c(round(q_cg, 1), round(q_chg, 1), round(q_chh, 1)))
    )
    print(knitr::kable(df_2_print))
    message("\n", paste(
        c(
            paste(names(df_2_print), collapse = "\t\t"),
            paste(sapply(names(df_2_print), function(nm) paste(rep("=", nchar(nm)), collapse = "")), collapse = "\t\t"),
            apply(df_2_print, 1, function(row) paste(row, collapse = "\t\t"))
        ),
        collapse = "\n"
    ), "\n")

    # genome density plot
    dHRs_circular_plot(cg_filtered, chg_filtered, chh_filtered, annotation.gr, TE.gr, paste0(trnt_name, "_vs_", ctrl_name))

    # manhattan plots
    manH_par <- mclapply(
        list(list(cg_hd, "CG", fdr), list(chg_hd, "CHG", fdr), list(chh_hd, "CHH", fdr)),
        function(x) {
            mannh_plots(x[[1]], x[[2]], x[[3]])
        },
        mc.cores = ifelse(num_cores >= 3, 3, num_cores)
    )
    png(paste0("sum_dH_manhattan_plot_", trnt_name, "_vs_", ctrl_name, ".png"), width = 10, height = 2.75, units = "in", res = 300, family = "serif")
    multiplot(
        print(manH_par[[1]]), print(manH_par[[2]]), print(manH_par[[3]]),
        cols = 3
    )
    dev.off()

    # write bigWig file function
    write_bw <- function(x, cntx) {
        x <- x %>%
            select(seqnames, start, end, sum_surprisal) %>%
            filter(!is.na(sum_surprisal)) %>%
            mutate(score = as.numeric(sum_surprisal)) %>%
            makeGRangesFromDataFrame(., keep.extra.columns = TRUE, )
        seqlengths(x) <- setNames(seqlengths_merged$seqlength, seqlengths_merged$seqnames)
        rtracklayer::export.bw(x, con = paste0("surprisal_", cntx, "_", trnt_name, "_vs_", ctrl_name, ".bw"))
    }

    # annotate regions and write files
    ann_list <- genome_ann(annotation.gr, TE.gr)
    write_sum_dH <- function(filtered_df, cntx) {
        write.csv(filtered_df, paste0("surprisal_", cntx, "_", trnt_name, "_vs_", ctrl_name, ".csv"), row.names = F)
        write_bw(filtered_df, cntx)
        setwd("genome_annotation")
        suppressMessages(DMRs_ann(ann_list, makeGRangesFromDataFrame(filtered_df, keep.extra.columns = T), cntx, description_df))
        setwd("../")
    }
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
            m_ctrl = sum(m_ctrl),
            n_ctrl = sum(n_ctrl),
            m_trnt = sum(m_trnt),
            n_trnt = sum(n_trnt),
            # adding 1 m_reads and 2 n_reads (total) if the count is 0
            pi_ctrl = (m_ctrl + (m_ctrl == 0)) / (n_ctrl + 2 * (m_ctrl == 0)),
            pi_trnt = (m_trnt + (m_trnt == 0)) / (n_trnt + 2 * (m_trnt == 0)),
            pi_log2FC = log2(pi_trnt / pi_ctrl),
            n_sites = n()
        ) %>%
        ungroup()

    # final df
    out <- cbind(as.data.frame(tiles100)[window_surprisal$tile_id, 1:3], window_surprisal) %>%
        filter(n_sites >= min_C_sites) %>%
        mutate(context = cntx) %>%
        dplyr::relocate(context, .before = sum_surprisal) %>%
        select(-tile_id)
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
            CG_df[CG_df$pi_log2FC > 0, 1:3],
            CG_df[CG_df$pi_log2FC < 0, 1:3]
        ),
        bg.col = "#fafcff", bg.border = NA, count_by = "number",
        col = c("#FF000080", "#304ed180"), border = T, track.height = 0.165, track.margin = c(0, 0)
    )

    circos.genomicDensity(
        list(
            CHG_df[CHG_df$pi_log2FC > 0, 1:3],
            CHG_df[CHG_df$pi_log2FC < 0, 1:3]
        ),
        bg.col = "#fafcff", bg.border = NA, count_by = "number",
        col = c("#FF000080", "#304ed180"), border = T, track.height = 0.165, track.margin = c(0, 0)
    )

    circos.genomicDensity(
        list(
            CHH_df[CHH_df$pi_log2FC > 0, 1:3],
            CHH_df[CHH_df$pi_log2FC < 0, 1:3]
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

mannh_plots <- function(df, cntx, fdr = 0.95) {
    S_threshold <- as.numeric(quantile(df$sum_surprisal, fdr))

    # add cumulative genome coordinate
    chr_lengths <- tapply(df$end, df$seqnames, max)
    chr_offset <- c(0, cumsum(chr_lengths[-length(chr_lengths)]))
    names(chr_offset) <- names(chr_lengths)

    df$coord <- df$start + chr_offset[df$seqnames]

    manH_p <- ggplot(df, aes(coord, sum_surprisal,
        colour = as.factor(as.integer(seqnames) %% 2)
    )) +
        geom_point(size = 0.025) +
        geom_hline(yintercept = S_threshold, color = "#bf6828", linetype = "dashed", size = 0.7) +
        scale_colour_manual(values = c("grey40", "#7ca182"), guide = "none") +
        scale_x_continuous(
            breaks = chr_offset + chr_lengths / 2,
            labels = paste0("Chr", seq_along(chr_offset))
        ) +
        expand_limits(y = 0) +
        labs(
            x = "", # "Chromosome",
            y = "Window surprisal"
        ) +
        theme_bw()

    manH_p
}