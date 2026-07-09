#' Plot a tissue/cell-type drug sensitivity matrix
#'
#' Builds one figure per GWAS by collecting compound-level CDR results across
#' tissue/cell-type model directories, recalculating Benjamini-Hochberg adjusted
#' compound p-values within each GWAS, and plotting the requested compounds as a
#' matrix. Rows are TWAS model IDs/tissues, columns are `pert_iname` values, and
#' point color and size encode either `Compound.pseudo.zscore` or within-model
#' rank percentile, depending on `plot.value`.
#'
#' @param twas.df TWAS summary statistics file or data frame. Used only to
#'   derive model labels and hierarchical clustering order. Defaults to
#'   `paste0("output/2.METAXCAN/", basename(getwd()),
#'   "_df.all.annotated.onlytranscripts.csv.gz")`.
#' @param drug.results.suffix Optional file name for compound-level CDR results
#'   inside each GWAS/model directory. Defaults to
#'   `"_cdr_all_compound_level.csv"`. Other options are e.g.
#'   `_cdr_P3_and_launched_compound_level.csv` for P3 and launched compounds.
#' @param results.dir Directory containing GWAS/model CDR result folders.
#'   Defaults to `"results/GTP_CDR/"`.
#' @param drug.list Character vector of `pert_iname` values to include.
#' @param plot.value Character. Either `"pseudoz"` or `"rank"`. `"pseudoz"`
#'   plots signed `Compound.pseudo.zscore`. `"rank"` converts `Rank` within each
#'   GWAS-model_ID combination to a percentile, where rank 1 is 100% and the
#'   worst rank is 0%.
#' @param tissue.order.method Character. Either `"cluster"` or `"alphabetical"`.
#'   Ignored when `custom.tissueorder` contains more than one non-NA value.
#' @param drug.order.method Character. Either `"cluster"` or `"alphabetical"`.
#'   Ignored when `custom.drugorder` contains more than one non-NA value.
#' @param add.ALL Logical. If `TRUE`, add an `ALL` row at the top when an ALL
#'   result file exists; otherwise summarize the loaded model rows by compound.
#'   Defaults to `FALSE`.
#' @param show.stars Logical. If `TRUE`, display stars marking BH-adjusted
#'   `Compound.MW.FDR` thresholds. If `FALSE`, suppress stars.
#' @param simplify.model.labels Logical. If `TRUE`, simplify y-axis model labels
#'   by keeping only the text before the first `" ::"`.
#' @param custom.tissueorder Character vector with a custom model_ID/tissue
#'   order. If left as `NA`, rows are ordered by hierarchical clustering.
#' @param custom.drugorder Character vector with a custom `pert_iname` order. If
#'   left as `NA`, columns are ordered by hierarchical clustering of
#'   `Compound.pseudo.zscore` across models.
#' @param custom.height Numeric plot height in inches. If `NA`, height is chosen
#'   from the number of model rows.
#' @param custom.width Numeric plot width in inches. If `NA`, width is chosen
#'   from the number of drugs.
#' @param gwas.list Optional character vector of GWAS IDs to plot. If `NA`, all
#'   GWAS directories under `results.dir` are considered.
#' @param twas.trait.column Column in `twas.df` containing GWAS IDs.
#' @param twas.tissue.column Column in `twas.df` containing model IDs/tissues.
#' @param twas.feature.column Column in `twas.df` containing features used for
#'   tissue clustering. If absent, `feature_name` and then `gene` are tried.
#' @param twas.zscore.column Column in `twas.df` containing TWAS z-scores.
#' @param output.prefix Prefix for output files.
#' @param plot.date Date string to append to output files. Defaults to today's
#'   date.
#' @param dpi PNG resolution. Defaults to 600.
#' @return Invisibly returns a named list with each GWAS plot, plotting data,
#'   and output file paths.
#' @export
plot_tissue_drug_sensitivity_matrix <- function(
    twas.df = paste0(
      "output/2.METAXCAN/", basename(getwd()),
      "_df.all.annotated.onlytranscripts.csv.gz"),
    drug.results.suffix = "_cdr_all_compound_level.csv",
    results.dir = "results/GTP_CDR/",
    drug.list,
    plot.value = c("rank", "pseudoz"),
    tissue.order.method = c("cluster", "alphabetical"),
    drug.order.method = c("cluster", "alphabetical"),
    add.ALL = FALSE,
    show.stars = FALSE,
    simplify.model.labels = TRUE,
    custom.tissueorder = NA,
    custom.drugorder = NA,
    custom.height = NA,
    custom.width = NA,
    gwas.list = NA,
    twas.trait.column = "gwas",
    twas.tissue.column = "model_ID",
    twas.feature.column = "feature",
    twas.zscore.column = "zscore",
    output.prefix = "Model.Source.Drugs",
    plot.date = format(Sys.Date(), "%Y-%m-%d"),
    dpi = 600
) {
  if (missing(drug.list) || length(drug.list) < 1) {
    stop("Please provide at least one pert_iname in drug.list.")
  }
  if (!dir.exists(results.dir)) {
    stop("results.dir does not exist: ", results.dir)
  }
  plot.value <- match.arg(plot.value)
  tissue.order.method <- match.arg(tissue.order.method)
  drug.order.method <- match.arg(drug.order.method)

  if (!is.logical(show.stars) || length(show.stars) != 1 || is.na(show.stars)) {
    stop("show.stars must be TRUE or FALSE.")
  }

  if (
    !is.logical(simplify.model.labels) ||
    length(simplify.model.labels) != 1 ||
    is.na(simplify.model.labels)
  ) {
    stop("simplify.model.labels must be TRUE or FALSE.")
  }

  drug.list <- unique(as.character(drug.list))
  results.dir <- sub("/+$", "", results.dir)

  twas <- .prepare_twas_for_sensitivity_matrix(
    twas.df = twas.df,
    twas.trait.column = twas.trait.column,
    twas.tissue.column = twas.tissue.column,
    twas.feature.column = twas.feature.column,
    twas.zscore.column = twas.zscore.column
  )
  on.exit(rm(twas), add = TRUE)

  gwas.dirs <- list.dirs(results.dir, full.names = TRUE, recursive = FALSE)
  gwas.ids <- basename(gwas.dirs)
  gwas.ids <- gwas.ids[gwas.ids != "intermediate.files"] # Exclude internal working directory
  if (!is.na(gwas.list[1])) {
    gwas.ids <- gwas.ids[gwas.ids %in% gwas.list]
  }
  if (length(gwas.ids) < 1) {
    stop("No GWAS directories found to plot under results.dir.")
  }

  output <- list()
  for (this.gwas in gwas.ids) {
    message("Generating tissue-drug sensitivity matrix for GWAS: ", this.gwas)
    # this loads all, no filtering.
    model.results <- .load_gwas_cdr_compound_results(
      results.dir = results.dir,
      this.gwas = this.gwas,
      drug.results.suffix = drug.results.suffix
    )
    if (nrow(model.results) < 1) {
      warning("No requested compounds found for GWAS: ", this.gwas)
      next
    }


    # Map Java-safe model_IDs back to original TWAS model_IDs
    id.map <- unique(
      twas[gwas == this.gwas, .(
        model_ID_java = make_java_safe(model_ID),
        model_ID_twas = model_ID
      )]
    )
    if (nrow(id.map) == 0) {
      stop("No TWAS model_IDs found for this.gwas.")
    }
    bad.map <- id.map[, .N, by = model_ID_java][N > 1]
    if (nrow(bad.map) > 0) {
      stop("Some Java-safe model_ID values map to multiple TWAS model_ID values.")
    }
    allowed.unmatched <- "ALL" # this for the aggregation of all tissues
    missing <- setdiff(model.results$model_ID, c(id.map$model_ID_java, allowed.unmatched))
    if (length(missing) > 0) {
      stop("Some model.results$model_ID values were not found in TWAS after make_java_safe().")
    }
    model.results[
      id.map,
      model_ID := i.model_ID_twas,
      on = .(model_ID = model_ID_java)
    ]

    # Remove ALL if we are done with it
    if (!add.ALL) {
      model.results <- model.results[model_ID != "ALL"]
    }

    # Estimate Compound.MW.FDR now that we have our final interest list
    model.results[, Compound.MW.FDR := p.adjust(Compound.MW.p, method = "BH")]


    # Add plotting value: pseudo-z or rank percentile
    model.results <- .add_sensitivity_plot_values(
      model.results = model.results,
      plot.value = plot.value
    )

    # Limit to compounds for plotting
    plot.dt <- model.results[pert_iname %in% drug.list]

    if (nrow(plot.dt) < 1) {
      warning("No requested compounds survived filtering for GWAS: ", this.gwas)
      next
    }

    ##### TISSUE ORDER BLOCK #####


    has.custom.tissueorder <- length(custom.tissueorder[!is.na(custom.tissueorder)]) > 1
    has.custom.drugorder <- length(custom.drugorder[!is.na(custom.drugorder)]) > 1

    if (!has.custom.tissueorder && tissue.order.method == "alphabetical") {
      tissue.order <- sort(unique(as.character(plot.dt$model_ID)))
      tissue.order <- setdiff(tissue.order, "ALL")

      if (isTRUE(add.ALL) && "ALL" %in% plot.dt$model_ID) {
        tissue.order <- c("ALL", tissue.order)
      }
    } else {
      tissue.order <- .get_sensitivity_tissue_order(
        twas = twas[gwas == this.gwas],
        plot.dt = plot.dt,
        custom.tissueorder = custom.tissueorder,
        add.ALL = add.ALL
      )
    }

    ##### DRUG ORDER BLOCK #####

    if (!has.custom.drugorder && drug.order.method == "alphabetical") {
      drug.order <- sort(unique(as.character(plot.dt$pert_iname)))
    } else {
      drug.order <- .get_sensitivity_drug_order(
        plot.dt = plot.dt,
        custom.drugorder = custom.drugorder,
        drug.list = drug.list
      )
    }

    plot.dt <- plot.dt[
      model_ID %in% tissue.order & pert_iname %in% drug.order
    ]

    if (nrow(plot.dt) < 1) {
      warning("No plot rows remain after applying tissue and drug order for GWAS: ", this.gwas)
      next
    }


    plot.dt <- .add_sensitivity_model_labels(
      plot.dt = plot.dt,
      tissue.order = tissue.order,
      simplify.model.labels = simplify.model.labels
    )

    plot.dt[, pert_iname := factor(pert_iname, levels = drug.order)]

    if (isTRUE(show.stars)) {
      plot.dt[, significance := .compound_fdr_to_stars(Compound.MW.FDR)]
    } else {
      plot.dt[, significance := ""]
    }

    plot.width <- ifelse(
      is.na(custom.width),
      max(7, min(36, 3.5 + 0.35 * length(drug.order))),
      custom.width
    )
    plot.height <- ifelse(
      is.na(custom.height),
      max(4, min(36, 1.8 + 0.28 * length(tissue.order))),
      custom.height
    )

    plot.limit <- max(plot.dt$Plot.abs.value, na.rm = TRUE)
    if (!is.finite(plot.limit) || plot.limit == 0) plot.limit <- 1

    if (plot.value == "pseudoz") {
      plot.limit <- ceiling(plot.limit)
    } else if (plot.value == "rank") {
      plot.limit <- 100
    }

##### PLOTTING STARTS HERE #####

    matrix.plot <- ggplot2::ggplot(
      plot.dt,
      ggplot2::aes(x = pert_iname, y = model_ID_plot)
    ) +
      ggplot2::geom_tile(
        fill = "white",
        color = "grey85",
        linewidth = 0.25
      )

    if (plot.value == "pseudoz") {
      matrix.plot <- matrix.plot +
        ggplot2::geom_point(
          ggplot2::aes(
            size = Plot.abs.value,
            fill = Plot.value
          ),
          shape = 22,
          color = "grey80",
          stroke = 0.15
        ) +
        ggplot2::scale_fill_gradient2(
          name = "Compound\npseudo z-score",
          low = "#053061",
          mid = "#F7F7F7",
          high = "#67001F",
          midpoint = 0,
          limits = c(-plot.limit, plot.limit),
          oob = scales::squish
        ) +
        ggplot2::scale_size_continuous(
          name = "|Compound\npseudo z-score|",
          range = c(0.4, 5.5),
          limits = c(0, plot.limit)
        )
    }

    if (plot.value == "rank") {
      matrix.plot <- matrix.plot +
        ggplot2::geom_tile(
          ggplot2::aes(fill = Plot.value),
          color = "grey85",
          linewidth = 0.25
        ) +
        ggplot2::scale_fill_gradient2(
          name = "Rank\npercentile",
          low = "#053061",
          mid = "#F7F7F7",
          high = "#67001F",
          midpoint = 50,
          limits = c(0, 100),
          labels = function(x) paste0(x, "%"),
          oob = scales::squish
        )
    }

    if (isTRUE(show.stars)) {
      matrix.plot <- matrix.plot +
        ggplot2::geom_text(
          ggplot2::aes(label = significance),
          color = "black",
          size = 2.8,
          vjust = 0.5,
          hjust = 0.5
        )
    }

    matrix.plot <- matrix.plot +
      ggplot2::labs(
        x = "Perturbagen",
        y = "GFI Model",
        title = paste0(this.gwas, " CDR tissue/cell-type sensitivity matrix"),
        subtitle = paste0(
          ifelse(
            plot.value == "pseudoz",
            "Fill/size show compound pseudo z-score",
            "Fill shows within-model rank percentile"
          ),
          ifelse(
            isTRUE(show.stars),
            "; stars indicate BH-adjusted Compound.MW.FDR thresholds",
            ""
          )
        )
      ) +
      ggplot2::coord_fixed(ratio = 1) +
      ggplot2::theme_classic(base_size = 10) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(
          angle = 45,
          hjust = 1,
          vjust = 1,
          color = "black"
        ),
        axis.text.y = ggplot2::element_text(color = "black"),
        axis.line = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(face = "bold"),
        legend.position = "right"
      )

    fig.dir <- file.path(results.dir, this.gwas, "Figures")
    gv_dir.create(fig.dir)
    pdf.file <- file.path(
      fig.dir,
      paste0(output.prefix, "_", plot.date, ".pdf")
    )
    png.file <- file.path(
      fig.dir,
      paste0(output.prefix, "_", plot.date, ".png")
    )

    ggplot2::ggsave(
      filename = pdf.file,
      plot = matrix.plot,
      width = plot.width,
      height = plot.height,
      units = "in"
    )
    ggplot2::ggsave(
      filename = png.file,
      plot = matrix.plot,
      width = plot.width,
      height = plot.height,
      units = "in",
      dpi = dpi
    )

    output[[this.gwas]] <- list(
      plot = matrix.plot,
      data = plot.dt,
      pdf = pdf.file,
      png = png.file,
      width = plot.width,
      height = plot.height
    )
  }

  invisible(output)
}

.prepare_twas_for_sensitivity_matrix <- function(
    twas.df,
    twas.trait.column,
    twas.tissue.column,
    twas.feature.column,
    twas.zscore.column
) {
  twas <- return_df(twas.df)
  twas <- data.table::as.data.table(twas)
  if (!(twas.trait.column %in% names(twas))) {
    stop("Could not find TWAS trait/GWAS column: ", twas.trait.column)
  }
  if (!(twas.tissue.column %in% names(twas))) {
    stop("Could not find TWAS tissue/model column: ", twas.tissue.column)
  }
  if (!(twas.zscore.column %in% names(twas))) {
    stop("Could not find TWAS z-score column: ", twas.zscore.column)
  }
  if (!(twas.feature.column %in% names(twas))) {
    fallback.feature <- c("feature_name", "gene")
    fallback.feature <- fallback.feature[fallback.feature %in% names(twas)]
    if (length(fallback.feature) < 1) {
      stop("Could not find TWAS feature column: ", twas.feature.column)
    }
    twas.feature.column <- fallback.feature[1]
  }
  data.table::setnames(
    twas,
    c(
      twas.trait.column,
      twas.tissue.column,
      twas.feature.column,
      twas.zscore.column
    ),
    c("gwas", "model_ID", "feature", "zscore"),
    skip_absent = TRUE
  )
  twas <- twas[
    !is.na(gwas) & !is.na(model_ID) & !is.na(feature) & !is.na(zscore),
    .(gwas, model_ID, feature, zscore)
  ]
  twas[, zscore := as.numeric(zscore)]
  twas
}

.load_gwas_cdr_compound_results <- function(
    results.dir,
    this.gwas,
    drug.results.suffix = "_cdr_all_compound_level.csv"
) {
  gwas.dir <- file.path(results.dir, this.gwas)
  model.dirs <- list.dirs(gwas.dir, full.names = TRUE, recursive = FALSE)
  model.dirs <- model.dirs[basename(model.dirs) != "Figures"]
  model.dirs <- model.dirs[basename(model.dirs) != "ALL"]

  if (length(model.dirs) < 1) {
    return(data.table::data.table())
  }

  drug.results.suffix <- as.character(drug.results.suffix[1])

  if (!startsWith(drug.results.suffix, "_")) {
    drug.results.suffix <- paste0("_", drug.results.suffix)
  }

  data.table::rbindlist(
    lapply(model.dirs, function(model.dir) {
      model.id <- basename(model.dir)
      this.file <- file.path(
        model.dir,
        paste0(this.gwas, drug.results.suffix)
      )

      if (!file.exists(this.file)) {
        warning("Skipping missing CDR file: ", this.file)
        return(data.table::data.table())
      }

      x <- data.table::fread(this.file)

      needed <- c("pert_iname", "Compound.MW.p", "Compound.pseudo.zscore")
      if (any(!(needed %in% names(x)))) {
        warning("Skipping file with missing required columns: ", this.file)
        return(data.table::data.table())
      }

      if (nrow(x) < 1) {
        return(data.table::data.table())
      }

      x[, gwas := this.gwas]
      x[, model_ID := model.id]

      x
    }),
    fill = TRUE
  )
}

.load_or_summarize_all_cdr_results <- function(
    results.dir,
    this.gwas,
    model.results
) {
  all.dir <- file.path(results.dir, this.gwas, "ALL")
  this.file <- file.path(all.dir, paste0(this.gwas, "_cdr_all_compound_level.csv"))
  if (!is.na(this.file)) {
    x <- data.table::fread(this.file)
    needed <- c("pert_iname", "Compound.MW.p", "Compound.pseudo.zscore")
    if (all(needed %in% names(x))) {
      x[, gwas := this.gwas]
      x[, model_ID := "ALL"]
      x[, Compound.MW.p := as.numeric(Compound.MW.p)]
      x[, Compound.pseudo.zscore := as.numeric(Compound.pseudo.zscore)]
      if (!("Rank" %in% names(x))) {
        x[
          ,
          Rank := data.table::frank(
            -Compound.pseudo.zscore,
            ties.method = "average",
            na.last = "keep"
          )
        ]
      }
      return(x)
    }
    warning("ALL file is missing required columns and will be ignored: ", this.file)
  }

  summary.dt <- model.results[
    ,
    .(
      Compound.MW.p = min(Compound.MW.p, na.rm = TRUE),
      Compound.pseudo.zscore = mean(Compound.pseudo.zscore, na.rm = TRUE),
      N_models = data.table::uniqueN(model_ID)
    ),
    by = .(gwas, pert_iname)
  ]
  summary.dt[!is.finite(Compound.MW.p), Compound.MW.p := NA_real_]
  summary.dt[!is.finite(Compound.pseudo.zscore), Compound.pseudo.zscore := NA_real_]
  summary.dt[, model_ID := "ALL"]

  summary.dt[
    ,
    Rank := data.table::frank(
      -Compound.pseudo.zscore,
      ties.method = "average",
      na.last = "keep"
    ),
    by = gwas
  ]

  summary.dt
}

##### .get_sensitivity_tissue_order #####

.get_sensitivity_tissue_order <- function(
    twas,
    plot.dt,
    custom.tissueorder,
    add.ALL
) {
  model.ids <- unique(as.character(plot.dt$model_ID))
  non.all.model.ids <- setdiff(model.ids, "ALL")

  if (length(custom.tissueorder) > 1) {

    ordered.models <- c(
      as.character(custom.tissueorder[custom.tissueorder %in% model.ids]),
      setdiff(model.ids, custom.tissueorder)
    )

  } else if (length(non.all.model.ids) > 1 && nrow(twas) > 0) {

    hc.dt <- twas[model_ID %in% non.all.model.ids]

    hc.dt <- hc.dt[
      ,
      .(zscore = mean(zscore, na.rm = TRUE)),
      by = .(model_ID, feature)
    ]

    complete.features <- hc.dt[
      !is.na(zscore),
      .(n_models = uniqueN(model_ID)),
      by = feature
    ][
      n_models == length(non.all.model.ids),
      feature
    ]

    hc.dt <- hc.dt[feature %in% complete.features]

    if (length(complete.features) > 0) {

      wide <- data.table::dcast(
        hc.dt,
        model_ID ~ feature,
        value.var = "zscore"
      )

      row.names <- wide$model_ID
      mat <- as.matrix(wide[, setdiff(names(wide), "model_ID"), with = FALSE])
      rownames(mat) <- row.names

      if (anyNA(mat)) {
        stop("Unexpected NA values remain after restricting to complete features.")
      }

      if (nrow(mat) > 1) {
        ordered.models <- rownames(mat)[
          stats::hclust(stats::dist(mat), method = "ward.D2")$order
        ]
      } else {
        ordered.models <- row.names
      }

    } else {

      wide <- data.table::dcast(
        plot.dt[model_ID != "ALL"],
        model_ID ~ pert_iname,
        value.var = "Compound.pseudo.zscore",
        fill = 0
      )

      row.names <- wide$model_ID
      mat <- as.matrix(wide[, setdiff(names(wide), "model_ID"), with = FALSE])
      rownames(mat) <- row.names
      mat[is.na(mat)] <- 0

      ordered.models <- rownames(mat)[
        stats::hclust(stats::dist(mat), method = "ward.D2")$order
      ]
    }

  } else if (length(non.all.model.ids) > 1) {

    wide <- data.table::dcast(
      plot.dt[model_ID != "ALL"],
      model_ID ~ pert_iname,
      value.var = "Compound.pseudo.zscore",
      fill = 0
    )

    row.names <- wide$model_ID
    mat <- as.matrix(wide[, setdiff(names(wide), "model_ID"), with = FALSE])
    rownames(mat) <- row.names
    mat[is.na(mat)] <- 0

    ordered.models <- rownames(mat)[
      stats::hclust(stats::dist(mat), method = "ward.D2")$order
    ]

  } else {

    ordered.models <- non.all.model.ids
  }

  if (isTRUE(add.ALL) && "ALL" %in% model.ids) {
    ordered.models <- c("ALL", setdiff(ordered.models, "ALL"))
  }

  ordered.models
}

.get_sensitivity_drug_order <- function(plot.dt, custom.drugorder, drug.list) {
  drug.ids <- unique(as.character(plot.dt$pert_iname))

  .scale_safe <- function(x) {
    sx <- stats::sd(x, na.rm = TRUE)
    if (is.na(sx) || sx == 0) {
      return(rep(0, length(x)))
    }
    as.numeric((x - mean(x, na.rm = TRUE)) / sx)
  }

  if (length(custom.drugorder) > 1) {
    return(c(
      as.character(custom.drugorder[custom.drugorder %in% drug.ids]),
      setdiff(drug.ids, custom.drugorder)
    ))
  }

  if (length(drug.ids) <= 2) {
    return(c(
      drug.list[drug.list %in% drug.ids],
      setdiff(drug.ids, drug.list)
    ))
  }

  hc.dt <- plot.dt[
    ,
    .(
      Compound.pseudo.zscore = mean(Compound.pseudo.zscore, na.rm = TRUE)
    ),
    by = .(pert_iname, model_ID)
  ]

  complete.models <- hc.dt[
    !is.na(Compound.pseudo.zscore),
    .(n_drugs = uniqueN(pert_iname)),
    by = model_ID
  ][
    n_drugs == length(drug.ids),
    model_ID
  ]

  hc.dt <- hc.dt[model_ID %in% complete.models]

  if (length(complete.models) > 0) {
    hc.dt[
      ,
      Compound.pseudo.zscore.scaled := .scale_safe(Compound.pseudo.zscore),
      by = model_ID
    ]

    wide <- data.table::dcast(
      hc.dt,
      pert_iname ~ model_ID,
      value.var = "Compound.pseudo.zscore.scaled"
    )

    row.names <- wide$pert_iname
    mat <- as.matrix(wide[, setdiff(names(wide), "pert_iname"), with = FALSE])
    rownames(mat) <- row.names

    if (anyNA(mat)) {
      stop("Unexpected NA values remain after restricting to complete model_IDs and scaling.")
    }

    return(rownames(mat)[
      stats::hclust(stats::dist(mat), method = "ward.D2")$order
    ])
  }

  c(
    drug.list[drug.list %in% drug.ids],
    setdiff(drug.ids, drug.list)
  )
}

.compound_fdr_to_stars <- function(fdr) {
  data.table::fifelse(
    !is.na(fdr) & fdr <= 0.0001, "****",
    data.table::fifelse(
      !is.na(fdr) & fdr <= 0.001, "***",
      data.table::fifelse(
        !is.na(fdr) & fdr <= 0.01, "**",
        data.table::fifelse(!is.na(fdr) & fdr <= 0.05, "*", "")
      )
    )
  )
}

.add_sensitivity_plot_values <- function(model.results, plot.value) {
  if (plot.value == "pseudoz") {
    model.results[, Plot.value := as.numeric(Compound.pseudo.zscore)]
    model.results[, Plot.abs.value := abs(Plot.value)]
    model.results[, Plot.fill.label := "Compound\npseudo z-score"]
    model.results[, Plot.size.label := "|Compound\npseudo z-score|"]
    return(model.results)
  }

  if (plot.value == "rank") {
    if (!("Rank" %in% names(model.results))) {
      stop("plot.value = 'rank' requires a Rank column in model.results.")
    }

    model.results[, Rank := as.numeric(Rank)]

    model.results[
      ,
      Rank.percentile := {
        max.rank <- suppressWarnings(max(Rank, na.rm = TRUE))

        if (!is.finite(max.rank)) {
          rep(NA_real_, .N)
        } else if (max.rank <= 1) {
          rep(100, .N)
        } else {
          100 * (1 - ((Rank - 1) / (max.rank - 1)))
        }
      },
      by = .(gwas, model_ID)
    ]

    model.results[, Rank.percentile := pmin(100, pmax(0, Rank.percentile))]
    model.results[, Plot.value := Rank.percentile]
    model.results[, Plot.abs.value := Rank.percentile]
    model.results[, Plot.fill.label := "Rank\npercentile"]
    model.results[, Plot.size.label := "Rank\npercentile"]
    return(model.results)
  }

  stop("Unsupported plot.value: ", plot.value)
}

.add_sensitivity_model_labels <- function(
    plot.dt,
    tissue.order,
    simplify.model.labels
) {
  label.map <- data.table::data.table(
    model_ID_raw = as.character(tissue.order),
    model_ID_plot = as.character(tissue.order)
  )

  if (isTRUE(simplify.model.labels)) {
    label.map[, model_ID_plot := sub("\\s*::.*$", "", model_ID_plot)]

    if (anyDuplicated(label.map$model_ID_plot)) {
      warning(
        "Simplified model_ID labels are not unique; adding suffixes to preserve plotting rows."
      )
      label.map[, model_ID_plot := make.unique(model_ID_plot, sep = " #")]
    }
  }

  plot.dt[, model_ID_raw := as.character(model_ID)]

  plot.dt[
    label.map,
    model_ID_plot := i.model_ID_plot,
    on = .(model_ID_raw)
  ]

  if (anyNA(plot.dt$model_ID_plot)) {
    stop("Some model_ID values could not be mapped to plotting labels.")
  }

  plot.dt[
    ,
    model_ID_plot := factor(
      model_ID_plot,
      levels = rev(label.map$model_ID_plot)
    )
  ]

  plot.dt
}
