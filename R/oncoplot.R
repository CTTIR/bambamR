#' Create Publication-Ready Oncoplot
#'
#' Generates a waterfall-style oncoplot showing the mutation landscape
#' across samples. Accepts either a simple `data.frame` with `sample`,
#' `gene`, `mutation_type` columns or a MAF-format data.frame.
#'
#' @param data A data.frame with columns `sample`, `gene`, `mutation_type`.
#'   Alternatively, a MAF-like data.frame with `Hugo_Symbol`,
#'   `Tumor_Sample_Barcode`, and `Variant_Classification` columns.
#' @param genes Character vector. Specific genes to display. If `NULL`,
#'   the top `n_genes` most frequently mutated genes are shown. Duplicated
#'   entries are shown once.
#' @param n_genes Integer. Number of top mutated genes if `genes` is `NULL`.
#'   Must be a positive whole number, or `Inf` for all genes. Default `20`.
#' @param sort_by Character. How to sort samples: `"frequency"` (by number of
#'   mutated genes shown) or `"cluster"` (by co-occurrence clustering).
#'   Default `"frequency"`.
#' @param annotation_df A data.frame for sample annotations (bottom tracks).
#'   Rownames must be sample identifiers. Each column becomes an annotation
#'   track. Annotation tracks require `patchwork`. A warning is raised when
#'   no row name matches a sample, because the tracks are then omitted.
#' @param mutation_colors Named character vector mapping mutation types to
#'   colours. If `NULL`, a default nature-style palette is used. Mutation
#'   types without an entry take the `"Other"` colour of this vector or, if it
#'   has none, a neutral grey.
#' @param show_pct Logical. Show mutation percentage per gene on the right
#'   margin. Default `TRUE`.
#' @param show_barplot Logical. Show top (mutated genes per sample) and side
#'   (mutated samples per gene) barplots. Default `TRUE`.
#' @param title Character or NULL. Plot title.
#' @param border_color Character or NULL. Colour of tile borders. Default
#'   `"white"`.
#'
#' @return A [ggplot2::ggplot] object (composed with patchwork if barplots
#'   are shown and patchwork is available).
#'
#' @details
#' The oncoplot encodes the mutation landscape of a cohort as a tiled grid:
#' columns = samples, rows = genes. Tiles are coloured by mutation type.
#' Multi-hit genes (multiple mutation types in one sample) are labelled
#' `"Multi_Hit"`.
#'
#' When `show_barplot = TRUE` and the `patchwork` package is available,
#' the plot is composed of three panels: top barplot (mutated genes per
#' sample), main tile plot, and right barplot (mutated samples per gene).
#' Without patchwork, the main tile plot is returned unless annotations were
#' requested. Percentages use all samples represented in `data` before gene
#' selection; samples with no recorded mutations cannot be inferred.
#'
#' The legend is placed below the plot in three columns so that it fits
#' common device widths. When moving it, reset the layout as well. A plot
#' with marginal barplots or annotation tracks is a patchwork composition, in
#' which `&` applies a change to every panel:
#' `p & ggplot2::theme(legend.position = "right") &
#' ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1))`. For a single
#' ggplot (`show_barplot = FALSE` and no `annotation_df`), use `+` instead
#' of `&`.
#'
#' @examples
#' set.seed(42)
#' mut_data <- data.frame(
#'   sample = sample(paste0("TCGA-", 1:30), 150, replace = TRUE),
#'   gene = sample(c("TP53","KRAS","PIK3CA","PTEN","APC",
#'                    "BRAF","EGFR","NRAS","CDKN2A","RB1"), 150, replace = TRUE),
#'   mutation_type = sample(c("Missense_Mutation","Nonsense_Mutation",
#'                            "Frame_Shift_Del","Splice_Site"), 150, replace = TRUE)
#' )
#' bb_oncoplot(mut_data, n_genes = 10)
#'
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_manual labs theme
#' @importFrom ggplot2 theme_minimal element_blank element_text geom_bar geom_text
#' @export
bb_oncoplot <- function(data,
                        genes = NULL,
                        n_genes = 20L,
                        sort_by = c("frequency", "cluster"),
                        annotation_df = NULL,
                        mutation_colors = NULL,
                        show_pct = TRUE,
                        show_barplot = TRUE,
                        title = NULL,
                        border_color = "white") {

  sort_by <- match.arg(sort_by)

  # Default mutation color palette (nature-style)
  default_colors <- c(
    "Missense_Mutation"      = "#2166AC",
    "Nonsense_Mutation"      = "#B2182B",
    "Frame_Shift_Del"        = "#1B7837",
    "Frame_Shift_Ins"        = "#762A83",
    "Splice_Site"            = "#E08214",
    "In_Frame_Del"           = "#F4A582",
    "In_Frame_Ins"           = "#92C5DE",
    "Translation_Start_Site" = "#D6604D",
    "Multi_Hit"              = "#4393C3",
    "Other"                  = "#999999"
  )
  if (is.null(mutation_colors)) {
    mutation_colors <- default_colors
  } else if (!is.character(mutation_colors) || is.null(names(mutation_colors)) ||
               anyNA(mutation_colors) || any(!nzchar(names(mutation_colors)))) {
    stop("'mutation_colors' must be a named character vector of colors.",
         call. = FALSE)
  }
  if (!is.null(annotation_df) && !is.data.frame(annotation_df)) {
    stop("'annotation_df' must be a data.frame with samples as rownames.",
         call. = FALSE)
  }

  # 1. Standardize input -----------------------------------------------
  df <- .standardize_onco_input(data)
  all_samples <- unique(df$sample)
  n_samples <- length(all_samples)

  # 2. Handle multi-hit ------------------------------------------------
  df <- .handle_multi_hit(df)

  # 3. Select genes ----------------------------------------------------
  gene_freq <- sort(table(unique(df[, c("sample", "gene")])$gene),
                    decreasing = TRUE)
  if (!is.null(genes)) {
    if (is.factor(genes)) genes <- as.character(genes)
    if (!is.character(genes) || !length(genes) || anyNA(genes)) {
      stop("'genes' must be a character vector of gene names without missing values.",
           call. = FALSE)
    }
    gene_list <- unique(genes)
  } else {
    .validate_n(n_genes, "n_genes")
    gene_list <- names(gene_freq)[seq_len(min(n_genes, length(gene_freq)))]
  }
  df <- df[df$gene %in% gene_list, , drop = FALSE]

  if (nrow(df) == 0L) {
    stop("No mutations found for the specified genes.", call. = FALSE)
  }

  # All samples and genes
  # 4. Sort samples ----------------------------------------------------
  sample_order <- .sort_onco_samples(df, gene_list, sort_by)
  sample_order <- c(sample_order, setdiff(all_samples, sample_order))

  # 5. Gene order by frequency -----------------------------------------
  gene_counts <- table(factor(df$gene, levels = gene_list))
  gene_order <- names(sort(gene_counts, decreasing = FALSE))

  # 6. Factor levels for ordering --------------------------------------
  df$sample <- factor(df$sample, levels = sample_order)
  df$gene <- factor(df$gene, levels = gene_order)

  # Map unknown mutation types to "Other"
  unknown <- setdiff(unique(df$mutation_type), names(mutation_colors))
  if (length(unknown) > 0L) {
    other_color <- if ("Other" %in% names(mutation_colors)) {
      mutation_colors[["Other"]]
    } else {
      default_colors[["Other"]]
    }
    mutation_colors[unknown] <- other_color
  }

  # 7. Build main tile plot --------------------------------------------
  p_main <- ggplot2::ggplot(df, ggplot2::aes(x = sample, y = gene,
                                             fill = mutation_type)) +
    ggplot2::geom_tile(color = border_color, linewidth = 0.5) +
    ggplot2::scale_x_discrete(limits = sample_order, drop = FALSE) +
    ggplot2::scale_fill_manual(
      values = mutation_colors,
      name = "Mutation Type",
      drop = FALSE
    ) +
    # Three legend columns with the title on top stay inside a 7 in device
    # even for the longest default label ("Translation_Start_Site").
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = 3)) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.key.size = ggplot2::unit(0.4, "cm"),
      legend.text = ggplot2::element_text(size = 8),
      plot.margin = ggplot2::margin(5, 5, 5, 5)
    ) +
    .legend_title_top()

  # 8. Percentage annotation -------------------------------------------
  if (show_pct) {
    gene_pct <- as.data.frame(gene_counts, stringsAsFactors = FALSE)
    colnames(gene_pct) <- c("gene", "gene_count")
    gene_pct$pct <- round(gene_pct$gene_count / n_samples * 100, 0)
    gene_pct$pct_label <- paste0(gene_pct$pct, "%")
    gene_pct$gene <- factor(gene_pct$gene, levels = gene_order)

    p_pct <- ggplot2::ggplot(gene_pct, ggplot2::aes(x = 1, y = gene,
                                                    label = pct_label)) +
      ggplot2::geom_text(size = 3, hjust = 0.5) +
      ggplot2::theme_void() +
      ggplot2::theme(plot.margin = ggplot2::margin(5, 2, 5, 2))
  }

  # 9. Top barplot (mutated genes shown per sample) --------------------
  if (show_barplot) {
    sample_burden <- as.data.frame(
      table(factor(df$sample, levels = sample_order)),
      stringsAsFactors = FALSE
    )
    colnames(sample_burden) <- c("sample", "n_mutations")
    sample_burden$sample <- factor(sample_burden$sample, levels = sample_order)

    p_top <- ggplot2::ggplot(sample_burden,
                             ggplot2::aes(x = sample, y = n_mutations)) +
      ggplot2::geom_bar(stat = "identity", fill = "#636363", width = 0.8) +
      # The bars count genes, so the axis breaks are whole numbers.
      ggplot2::scale_y_continuous(breaks = function(limits) {
        breaks <- unique(round(pretty(limits)))
        breaks[breaks >= limits[1] & breaks <= limits[2]]
      }) +
      # Two lines: the panel is short, a one-line title would run past it.
      ggplot2::labs(y = "Mutated\ngenes") +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_text(size = 9),
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(5, 5, 0, 5)
      )

    # Right barplot (gene mutation count) --------------------------------
    gene_bar <- as.data.frame(gene_counts, stringsAsFactors = FALSE)
    colnames(gene_bar) <- c("gene", "gene_count")
    gene_bar$gene <- factor(gene_bar$gene, levels = gene_order)

    p_right <- ggplot2::ggplot(gene_bar,
                               ggplot2::aes(x = gene_count, y = gene)) +
      ggplot2::geom_bar(stat = "identity", fill = "#636363", width = 0.8) +
      ggplot2::labs(x = "Samples") +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(5, 5, 5, 0)
      )
  }

  # 10. Annotation tracks (bottom) ------------------------------------
  if (!is.null(annotation_df)) {
    anno_samples <- intersect(sample_order, rownames(annotation_df))
    if (length(anno_samples) > 0L) {
      anno_long <- .melt_annotation(annotation_df, anno_samples, sample_order)

      p_anno <- ggplot2::ggplot(anno_long,
                                ggplot2::aes(x = sample, y = variable,
                                             fill = annotation_value)) +
        ggplot2::geom_tile(color = border_color, linewidth = 0.5) +
        ggplot2::labs(fill = "Annotation") +
        ggplot2::theme_minimal(base_size = 10) +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 90, hjust = 1,
                                              vjust = 0.5, size = 6),
          axis.title = ggplot2::element_blank(),
          panel.grid = ggplot2::element_blank(),
          legend.position = "bottom",
          legend.key.size = ggplot2::unit(0.4, "cm"),
          legend.text = ggplot2::element_text(size = 8),
          plot.margin = ggplot2::margin(0, 5, 5, 5)
        ) +
        .legend_title_top()
    } else {
      warning("No row name of 'annotation_df' matches a sample; ",
              "annotation tracks are omitted.", call. = FALSE)
      annotation_df <- NULL
    }
  }

  # 11. Compose all panels on one grid so sample columns remain aligned.
  if (show_barplot && requireNamespace("patchwork", quietly = TRUE)) {
    p_top <- p_top + ggplot2::scale_x_discrete(limits = sample_order, drop = FALSE)
    panels <- list(T = p_top, M = p_main, R = p_right)
    design <- if (show_pct) "T##\nMPR" else "T#\nMR"
    widths <- if (show_pct) c(1, 0.12, 0.25) else c(1, 0.25)
    heights <- c(0.25, 1)
    if (show_pct) panels$P <- p_pct
    if (!is.null(annotation_df)) {
      panels$A <- p_anno + ggplot2::scale_x_discrete(limits = sample_order, drop = FALSE)
      design <- paste0(design, if (show_pct) "\nA##" else "\nA#")
      heights <- c(heights, 0.25)
    }
    # Collect the legends below the whole figure: a legend centred under the
    # (narrower) tile panel alone is clipped at the left device edge.
    layout <- patchwork::wrap_plots(panels, design = design,
                                    widths = widths, heights = heights,
                                    guides = "collect") &
      ggplot2::theme(legend.position = "bottom", legend.box = "vertical")
    if (!is.null(title)) layout <- layout + patchwork::plot_annotation(title = title)
    return(layout)
  }

  # Fallback: return main tile plot only
  if (!is.null(title)) {
    p_main <- p_main + ggplot2::labs(title = title)
  }
  if (show_pct) {
    labels <- stats::setNames(paste0(gene_pct$gene, " (", gene_pct$pct_label, ")"),
                              gene_pct$gene)
    p_main <- p_main +
      ggplot2::scale_y_discrete(labels = labels) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))
  }
  if (!is.null(annotation_df)) {
    check_pkg("patchwork", reason = "for oncoplot annotation tracks")
    p_anno <- p_anno + ggplot2::scale_x_discrete(limits = sample_order, drop = FALSE)
    return(
      p_main / p_anno +
        patchwork::plot_layout(heights = c(4, 1), guides = "collect") &
        ggplot2::theme(legend.position = "bottom", legend.box = "vertical")
    )
  }

  p_main
}


# Internal helpers for oncoplot
# ==================================================================

#' Standardize oncoplot input to sample/gene/mutation_type columns
#' @noRd
.standardize_onco_input <- function(data) {
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame.", call. = FALSE)
  }

  # Check for MAF format
  maf_cols <- c("Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification")
  simple_cols <- c("sample", "gene", "mutation_type")

  if (all(maf_cols %in% colnames(data))) {
    df <- data.frame(
      sample = as.character(data$Tumor_Sample_Barcode),
      gene = as.character(data$Hugo_Symbol),
      mutation_type = as.character(data$Variant_Classification),
      stringsAsFactors = FALSE
    )
  } else if (all(simple_cols %in% colnames(data))) {
    df <- data.frame(
      sample = as.character(data$sample),
      gene = as.character(data$gene),
      mutation_type = as.character(data$mutation_type),
      stringsAsFactors = FALSE
    )
  } else {
    stop(
      "'data' must have columns {sample, gene, mutation_type} or ",
      "{Hugo_Symbol, Tumor_Sample_Barcode, Variant_Classification}.",
      call. = FALSE
    )
  }

  # Remove empty rows
  df <- df[!is.na(df$gene) & !is.na(df$sample) &
             nchar(df$gene) > 0L & nchar(df$sample) > 0L, , drop = FALSE]
  if (anyNA(df$mutation_type)) {
    stop("Mutation types must not be missing.", call. = FALSE)
  }
  df
}


#' Handle multi-hit mutations (multiple types per sample-gene pair)
#' @noRd
.handle_multi_hit <- function(df) {
  # Find duplicated sample-gene pairs
  pairs <- df[, c("sample", "gene"), drop = FALSE]
  is_dup <- duplicated(pairs) | duplicated(pairs, fromLast = TRUE)

  if (any(is_dup)) {
    # Keep one row per duplicate, mark as Multi_Hit
    df_unique <- df[!is_dup, , drop = FALSE]
    df_multi <- unique(df[is_dup, c("sample", "gene"), drop = FALSE])
    df_multi$mutation_type <- "Multi_Hit"
    df <- rbind(df_unique, df_multi)
  }

  df
}


#' Sort oncoplot samples
#' @noRd
.sort_onco_samples <- function(df, gene_list, sort_by) {
  all_samples <- unique(df$sample)

  if (sort_by == "frequency") {
    # Sort by number of mutated genes shown (descending)
    burden <- table(df$sample)
    sample_order <- names(sort(burden, decreasing = TRUE))
  } else {
    # Cluster by co-occurrence
    # Build binary matrix (genes x samples)
    mat <- matrix(0L, nrow = length(gene_list), ncol = length(all_samples),
                  dimnames = list(gene_list, all_samples))
    keep <- df$gene %in% gene_list
    mat[cbind(as.character(df$gene[keep]), as.character(df$sample[keep]))] <- 1L
    if (ncol(mat) > 2L) {
      hc <- stats::hclust(stats::dist(t(mat), method = "binary"))
      sample_order <- colnames(mat)[hc$order]
    } else {
      sample_order <- all_samples
    }
  }

  sample_order
}


#' Melt annotation data.frame to long format
#' @noRd
.melt_annotation <- function(annotation_df, anno_samples, sample_order) {
  anno_sub <- annotation_df[anno_samples, , drop = FALSE]
  anno_list <- lapply(colnames(anno_sub), function(col) {
    data.frame(
      sample = rownames(anno_sub),
      variable = col,
      annotation_value = as.character(anno_sub[[col]]),
      stringsAsFactors = FALSE
    )
  })
  anno_long <- do.call(rbind, anno_list)
  anno_long$sample <- factor(anno_long$sample, levels = sample_order)
  anno_long
}
