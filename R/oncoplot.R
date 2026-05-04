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
#'   the top `n_genes` most frequently mutated genes are shown.
#' @param n_genes Integer. Number of top mutated genes if `genes` is `NULL`.
#'   Default `20`.
#' @param sort_by Character. How to sort samples: `"frequency"` (by mutation
#'   burden) or `"cluster"` (by co-occurrence clustering). Default `"frequency"`.
#' @param annotation_df A data.frame for sample annotations (bottom tracks).
#'   Rownames must be sample identifiers. Each column becomes an annotation
#'   track.
#' @param mutation_colors Named character vector mapping mutation types to
#'   colors. If `NULL`, a default nature-style palette is used.
#' @param show_pct Logical. Show mutation percentage per gene on the right
#'   margin. Default `TRUE`.
#' @param show_barplot Logical. Show top (sample burden) and side (gene count)
#'   barplots. Default `TRUE`.
#' @param title Character or NULL. Plot title.
#' @param border_color Character or NULL. Color of tile borders. Default
#'   `"white"`.
#'
#' @return A [ggplot2::ggplot] object (composed with patchwork if barplots
#'   are shown and patchwork is available).
#'
#' @details
#' The oncoplot encodes the mutation landscape of a cohort as a tiled grid:
#' columns = samples, rows = genes. Tiles are colored by mutation type.
#' Multi-hit genes (multiple mutation types in one sample) are labeled
#' `"Multi_Hit"`.
#'
#' When `show_barplot = TRUE` and the `patchwork` package is available,
#' the plot is composed of three panels: top barplot (mutation burden per
#' sample), main tile plot, and right barplot (mutations per gene).
#' Without patchwork, only the main tile plot is returned.
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
#'   theme_minimal element_blank element_text geom_bar geom_text
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
  }

  # 1. Standardize input -----------------------------------------------
  df <- .standardize_onco_input(data)

  # 2. Handle multi-hit ------------------------------------------------
  df <- .handle_multi_hit(df)

  # 3. Select genes ----------------------------------------------------
  gene_freq <- sort(table(unique(df[, c("sample", "gene")])$gene),
                    decreasing = TRUE)
  if (!is.null(genes)) {
    gene_list <- genes
  } else {
    gene_list <- names(gene_freq)[seq_len(min(n_genes, length(gene_freq)))]
  }
  df <- df[df$gene %in% gene_list, , drop = FALSE]

  if (nrow(df) == 0L) {
    stop("No mutations found for the specified genes.", call. = FALSE)
  }

  # All samples and genes
  all_samples <- unique(df$sample)
  n_samples <- length(all_samples)

  # 4. Sort samples ----------------------------------------------------
  sample_order <- .sort_onco_samples(df, gene_list, sort_by)

  # 5. Gene order by frequency -----------------------------------------
  gene_counts <- table(factor(df$gene, levels = gene_list))
  gene_order <- names(sort(gene_counts, decreasing = FALSE))

  # 6. Factor levels for ordering --------------------------------------
  df$sample <- factor(df$sample, levels = sample_order)
  df$gene <- factor(df$gene, levels = gene_order)

  # Map unknown mutation types to "Other"
  unknown <- setdiff(unique(df$mutation_type), names(mutation_colors))
  if (length(unknown) > 0L) {
    for (u in unknown) {
      mutation_colors[u] <- mutation_colors[["Other"]]
    }
  }

  # 7. Build main tile plot --------------------------------------------
  p_main <- ggplot2::ggplot(df, ggplot2::aes(x = sample, y = gene,
                                              fill = mutation_type)) +
    ggplot2::geom_tile(color = border_color, linewidth = 0.5) +
    ggplot2::scale_fill_manual(
      values = mutation_colors,
      name = "Mutation Type",
      drop = FALSE
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.key.size = ggplot2::unit(0.4, "cm"),
      plot.margin = ggplot2::margin(5, 5, 5, 5)
    )

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

  # 9. Top barplot (sample mutation burden) ----------------------------
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
      ggplot2::labs(y = "# Mutations") +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
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
      ggplot2::labs(x = "# Mutations") +
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
          plot.margin = ggplot2::margin(0, 5, 5, 5)
        )
    } else {
      annotation_df <- NULL
    }
  }

  # 11. Compose with patchwork -----------------------------------------
  if (show_barplot && requireNamespace("patchwork", quietly = TRUE)) {
    # Build layout
    empty <- ggplot2::ggplot() + ggplot2::theme_void()

    if (show_pct) {
      # Top row: empty + top barplot + empty + empty
      # Mid row: empty + main tile + pct + right barplot
      layout <- (
        (empty | p_top | empty | empty) /
        (empty | p_main | p_pct | p_right)
      ) +
        patchwork::plot_layout(
          widths = c(0.02, 1, 0.08, 0.15),
          heights = c(0.2, 1)
        )
    } else {
      layout <- (
        (empty | p_top | empty) /
        (empty | p_main | p_right)
      ) +
        patchwork::plot_layout(
          widths = c(0.02, 1, 0.15),
          heights = c(0.2, 1)
        )
    }

    if (!is.null(annotation_df)) {
      if (show_pct) {
        layout <- layout / (empty | p_anno | empty | empty) +
          patchwork::plot_layout(heights = c(0.2, 1, 0.15))
      } else {
        layout <- layout / (empty | p_anno | empty) +
          patchwork::plot_layout(heights = c(0.2, 1, 0.15))
      }
    }

    if (!is.null(title)) {
      layout <- layout + patchwork::plot_annotation(title = title)
    }

    return(layout)
  }

  # Fallback: return main tile plot only
  if (!is.null(title)) {
    p_main <- p_main + ggplot2::labs(title = title)
  }
  if (show_pct) {
    p_main <- p_main +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))
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
  df <- df[nchar(df$gene) > 0L & nchar(df$sample) > 0L, , drop = FALSE]
  df
}


#' Handle multi-hit mutations (multiple types per sample-gene pair)
#' @noRd
.handle_multi_hit <- function(df) {
  # Find duplicated sample-gene pairs
  sg <- paste(df$sample, df$gene, sep = "___")
  dup_sg <- unique(sg[duplicated(sg)])

  if (length(dup_sg) > 0L) {
    is_dup <- sg %in% dup_sg
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
    # Sort by total mutation burden (descending)
    burden <- table(df$sample)
    sample_order <- names(sort(burden, decreasing = TRUE))
  } else {
    # Cluster by co-occurrence
    # Build binary matrix (genes x samples)
    mat <- matrix(0L, nrow = length(gene_list), ncol = length(all_samples),
                  dimnames = list(gene_list, all_samples))
    for (i in seq_len(nrow(df))) {
      g <- as.character(df$gene[i])
      s <- as.character(df$sample[i])
      if (g %in% gene_list) mat[g, s] <- 1L
    }
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
