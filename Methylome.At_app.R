#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(shiny)
})

# ----------------------------
# Helpers
# ----------------------------
is_abs_path <- function(p) {
  grepl("^(/|[A-Za-z]:)", p)
}

split_paths <- function(x) {
  if (is.null(x) || !nzchar(x)) return(character())
  parts <- unlist(strsplit(x, "[;,]"))
  parts <- trimws(parts)
  parts[nzchar(parts)]
}

resolve_paths <- function(paths, root) {
  if (length(paths) == 0) return(character())
  vapply(paths, function(p) {
    if (is_abs_path(p)) {
      normalizePath(p, winslash = "/", mustWork = FALSE)
    } else {
      normalizePath(file.path(root, p), winslash = "/", mustWork = FALSE)
    }
  }, character(1))
}

rel_path <- function(p, root) {
  if (!nzchar(p)) return(p)
  rp <- normalizePath(p, winslash = "/", mustWork = FALSE)
  rr <- normalizePath(root, winslash = "/", mustWork = FALSE)
  sub(paste0("^", rr, "/?"), "", rp)
}

list_experiments <- function(root) {
  results_dir <- file.path(root, "results")
  if (!dir.exists(results_dir)) return(character())
  dirs <- list.dirs(results_dir, full.names = FALSE, recursive = FALSE)
  dirs <- dirs[nzchar(dirs)]
  sort(dirs)
}

parse_experiment <- function(name) {
  if (!nzchar(name)) return(NULL)
  parts <- strsplit(name, "_vs_", fixed = TRUE)[[1]]
  if (length(parts) != 2) return(NULL)
  list(var2 = parts[1], var1 = parts[2])
}

is_media <- function(x) {
  ext <- tolower(tools::file_ext(x))
  ext %in% c("png", "jpg", "jpeg", "svg", "gif", "pdf")
}

list_media_files <- function(dir_path, recursive = FALSE) {
  if (!dir.exists(dir_path)) return(character())
  files <- list.files(dir_path, full.names = TRUE, recursive = recursive)
  files <- files[file.info(files)$isdir %in% c(FALSE, NA)]
  files <- files[sapply(files, is_media)]
  files[order(basename(files))]
}

order_by_context <- function(files, contexts = c("CG", "CHG", "CHH", "all_contexts", "all")) {
  b <- basename(files)
  ctx <- rep(NA_integer_, length(files))
  for (i in seq_along(contexts)) {
    hit <- grepl(paste0("\\b", contexts[i], "\\b"), b, ignore.case = TRUE)
    ctx[is.na(ctx) & hit] <- i
  }
  ctx[is.na(ctx)] <- length(contexts) + 1
  files[order(ctx, b)]
}

build_paths <- function(root, comparison_dir) {
  exp_path <- file.path(root, "results", comparison_dir)

  total_meth_path <- file.path(exp_path, "total_methylation_analysis")
  list(
    exp_path = exp_path,
    comparison_dir = comparison_dir,
    total_meth_path = total_meth_path,
    PCA_plots_path = file.path(total_meth_path, "PCA_plots"),
    meth_levels_path = file.path(total_meth_path, "methylation_levels"),
    ChrPlot_CX_path = file.path(total_meth_path, "ChrPlot_CX"),
    ChrPlot_subCX_path = file.path(total_meth_path, "ChrPlot_CX", "subCX"),
    TEs_distance_n_size_path = file.path(total_meth_path, "TE_size_n_distance"),
    total_meth_annotation_path = file.path(total_meth_path, "total_methylation_annotations"),
    TF_motifs_path = file.path(total_meth_path, "TF_motifs"),
    DMRs_analysis_path = file.path(exp_path, "DMR_analysis"),
    gainORloss_path = file.path(exp_path, "DMR_analysis", "gain_OR_loss"),
    genome_ann_path = file.path(exp_path, "DMR_analysis", "genome_annotation"),
    ChrPlots_DMRs_path = file.path(exp_path, "DMR_analysis", "ChrPlot_DMRs"),
    DMRs_bigWig_path = file.path(exp_path, "DMR_analysis", "DMRs_bigWig"),
    strand_asymmetry_path = file.path(exp_path, "Strand_Asymmetry_DMRs"),
    DMV_analysis_path = file.path(exp_path, "DMV_analysis"),
    dH_CX_path = file.path(exp_path, "deltaH"),
    metaPlot_pathA = file.path(exp_path, "MetaPlots"),
    metaPlot_pathB = file.path(exp_path, "metaPlots")
  )
}

media_dirs_from_paths <- function(p) {
  # Named list with path + default recursion
  candidates <- list(
    "PCA plots" = list(path = p$PCA_plots_path, recursive = FALSE),
    "Methylation levels" = list(path = p$meth_levels_path, recursive = FALSE),
    "ChrPlot_CX" = list(path = p$ChrPlot_CX_path, recursive = TRUE),
    "TE size and distance" = list(path = p$TEs_distance_n_size_path, recursive = TRUE),
    "Total methylation annotations" = list(path = p$total_meth_annotation_path, recursive = TRUE),
    "TF motifs" = list(path = p$TF_motifs_path, recursive = TRUE),
    "DMR density / gain-loss" = list(path = p$gainORloss_path, recursive = TRUE),
    "DMR chromosome plots" = list(path = p$ChrPlots_DMRs_path, recursive = TRUE),
    "DMR genome annotation" = list(path = p$genome_ann_path, recursive = TRUE),
    "Strand asymmetry" = list(path = p$strand_asymmetry_path, recursive = TRUE),
    "DMV analysis" = list(path = p$DMV_analysis_path, recursive = TRUE),
    "dH analysis" = list(path = p$dH_CX_path, recursive = TRUE),
    "MetaPlots" = list(path = p$metaPlot_pathA, recursive = TRUE),
    "metaPlots" = list(path = p$metaPlot_pathB, recursive = TRUE)
  )

  ok <- vapply(candidates, function(x) dir.exists(x$path), logical(1))
  candidates[ok]
}

ui <- fluidPage(
  titlePanel("Methylome.At - Overview Browser"),
  sidebarLayout(
    sidebarPanel(
      textInput("root", "Methylome.At_path", value = "."),
      uiOutput("experiment_selector"),
      textInput("var1", "var1 (control) - optional if experiment selected", value = ""),
      textInput("var2", "var2 (treatment) - optional if experiment selected", value = ""),
      actionButton("load", "Load overview", class = "btn-primary"),
      tags$hr(),
      helpText("Input paths (optional, separate multiple with ; or ,)"),
      textInput("var1_path", "var1_path", value = ""),
      textInput("var2_path", "var2_path", value = ""),
      textInput("annotation_file", "annotation_file", value = "./annotation_files/Methylome.At_annotations.csv.gz"),
      textInput("description_file", "description_file", value = "./annotation_files/Methylome.At_annotations.csv.gz"),
      textInput("TEs_file", "TEs_file", value = "./annotation_files/Methylome.At_annotations.csv.gz")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel(
          "Overview",
          uiOutput("overview_text"),
          tags$h3("Input summary"),
          tableOutput("paths_table"),
          tags$h3("Conversion rate (ChrC)"),
          tableOutput("conversion_table")
        ),
        tabPanel(
          "Media",
          uiOutput("media_controls"),
          uiOutput("media_gallery")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  experiments <- reactiveVal(character())

  observeEvent(input$root, {
    root <- normalizePath(input$root, winslash = "/", mustWork = FALSE)
    experiments(list_experiments(root))
  }, ignoreNULL = FALSE)

  observeEvent(input$refresh_experiments, {
    root <- normalizePath(input$root, winslash = "/", mustWork = FALSE)
    experiments(list_experiments(root))
  })

  output$experiment_selector <- renderUI({
    choices <- experiments()
    if (length(choices) == 0) {
      return(tagList(
        tags$em("No experiments found under results/."),
        actionButton("refresh_experiments", "Refresh list")
      ))
    }

    tagList(
      selectInput("experiment", "Experiment (results/*)", choices = c("", choices), selected = input$experiment),
      actionButton("refresh_experiments", "Refresh list")
    )
  })

  observeEvent(input$experiment, {
    vars <- parse_experiment(trimws(input$experiment))
    if (!is.null(vars)) {
      updateTextInput(session, "var1", value = vars$var1)
      updateTextInput(session, "var2", value = vars$var2)
    }
  })

  params <- eventReactive(input$load, {
    root <- normalizePath(input$root, winslash = "/", mustWork = FALSE)
    exp_name <- trimws(input$experiment)
    vars <- parse_experiment(exp_name)

    var1 <- if (!is.null(vars)) vars$var1 else trimws(input$var1)
    var2 <- if (!is.null(vars)) vars$var2 else trimws(input$var2)
    comparison_dir <- if (nzchar(exp_name)) {
      exp_name
    } else if (nzchar(var1) && nzchar(var2)) {
      paste0(var2, "_vs_", var1)
    } else {
      ""
    }

    list(
      root = root,
      experiment = exp_name,
      comparison_dir = comparison_dir,
      var1 = var1,
      var2 = var2,
      var1_path = split_paths(input$var1_path),
      var2_path = split_paths(input$var2_path),
      annotation_file = split_paths(input$annotation_file),
      description_file = split_paths(input$description_file),
      TEs_file = split_paths(input$TEs_file)
    )
  }, ignoreNULL = TRUE)

  paths <- reactive({
    p <- params()
    req(p)
    if (!nzchar(p$comparison_dir)) return(NULL)
    build_paths(p$root, p$comparison_dir)
  })

  observeEvent(params(), {
    p <- params()
    req(p)
    if (dir.exists(p$root)) {
      addResourcePath("methylome", p$root)
    }
  })

  output$overview_text <- renderUI({
    p <- params()
    if (is.null(p)) {
      return(tags$em("Set inputs and click 'Load overview'."))
    }
    if (!nzchar(p$comparison_dir)) {
      return(tags$div(style = "color:#a00;", "Select an experiment or provide var1 and var2."))
    }
    if (!dir.exists(p$root)) {
      return(tags$div(style = "color:#a00;", paste("Methylome.At_path does not exist:", p$root)))
    }

    pp <- paths()
    if (is.null(pp) || !dir.exists(pp$exp_path)) {
      return(tags$div(
        style = "color:#a00;",
        paste("Results folder not found:", file.path(p$root, "results", p$comparison_dir))
      ))
    }

    tagList(
      tags$p(tags$strong("Experiment folder:"), " ", tags$code(p$comparison_dir)),
      if (nzchar(p$var1) && nzchar(p$var2)) {
        tags$p(tags$strong("Comparison:"), " ", tags$code(paste0(p$var2, " vs ", p$var1)))
      },
      tags$p(tags$strong("Results folder:"), " ", tags$code(pp$exp_path)),
      tags$p("This app uses the same argument names as Methylome.At_main() and loads only available outputs."),
      tags$p(
        tags$strong("Contexts:"), " CG / CHG / CHH refer to the sequence context around the cytosine. ",
        "Methylation is summarized per context because the underlying maintenance pathways differ."
      )
    )
  })

  output$paths_table <- renderTable({
    p <- params()
    req(p)
    items <- list(
      var1_path = p$var1_path,
      var2_path = p$var2_path,
      annotation_file = p$annotation_file,
      description_file = p$description_file,
      TEs_file = p$TEs_file
    )

    data.frame(
      Item = names(items),
      Path = vapply(items, function(x) {
        if (length(x) == 0) return(NA_character_)
        resolved <- resolve_paths(x, p$root)
        paste(rel_path(resolved, p$root), collapse = "; ")
      }, character(1)),
      Exists = vapply(items, function(x) {
        if (length(x) == 0) return(NA)
        all(file.exists(resolve_paths(x, p$root)))
      }, logical(1)),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  output$conversion_table <- renderTable({
    p <- params()
    req(p)
    pp <- paths()
    req(pp)
    conv_path <- file.path(pp$exp_path, "conversion_rate.csv")
    if (!file.exists(conv_path)) {
      return(data.frame(Message = paste("Missing file:", rel_path(conv_path, p$root)), check.names = FALSE))
    }
    read.csv(conv_path, check.names = FALSE)
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  output$media_controls <- renderUI({
    p <- params()
    if (is.null(p)) return(NULL)
    pp <- paths()
    if (is.null(pp) || !dir.exists(pp$exp_path)) {
      return(tags$em("Load a valid overview to browse media."))
    }

    dirs <- media_dirs_from_paths(pp)
    if (length(dirs) == 0) {
      return(tags$em("No media folders found for this comparison."))
    }

    tagList(
      selectInput("dir_choice", "Media folder", choices = names(dirs)),
      checkboxInput("media_recursive", "Search subfolders", value = TRUE)
    )
  })

  output$media_gallery <- renderUI({
    p <- params()
    if (is.null(p)) return(NULL)
    pp <- paths()
    if (is.null(pp)) return(NULL)

    dirs <- media_dirs_from_paths(pp)
    if (length(dirs) == 0) return(NULL)

    choice <- input$dir_choice
    if (is.null(choice) || !choice %in% names(dirs)) return(NULL)
    dir_info <- dirs[[choice]]

    files <- list_media_files(dir_info$path, recursive = isTRUE(input$media_recursive))
    if (length(files) == 0) return(tags$em("No media files found in this folder."))

    files <- order_by_context(files)
    tagList(lapply(files, function(f) {
      ext <- tolower(tools::file_ext(f))
      rel <- rel_path(f, p$root)
      url <- paste0("methylome/", URLencode(rel))

      tag <- if (ext == "pdf") {
        tags$iframe(src = url, style = "width:100%; height:650px; border:0;")
      } else {
        tags$img(src = url, style = "max-width:100%; height:auto; display:block;")
      }

      tags$div(
        style = "margin: 12px 0 24px 0;",
        tags$div(style = "font-weight: 600; margin-bottom: 6px;", rel),
        tag
      )
    }))
  })
}

shinyApp(ui, server)
