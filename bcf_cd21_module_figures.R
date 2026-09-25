# BCF anti-CD21 figures and the Nature source-data workbook.
#
# Authors: GW McElfresh, Daniel Newhouse, Leanne Whitmore
#
# Reads the deposited BCF voom log2 matrix, recomputes the row z-scores that
# color the heatmaps, and writes figures/ plus data/Source_data.xlsx.
# Figure 3A is the BCF heatmap. Figure 3B is the volcano, with anti-CD21 in
# the positive direction. Figure 3D is the immediate-early heatmap and score.
# Figure 3E is the Tfh heatmap and score. Supplemental Figure 11B is NF-kB
# and 11C is IFN/ISG.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
repo_root <- if (length(file_arg) == 1L) {
  dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  normalizePath(getwd())
}

data_dir <- file.path(repo_root, "data")
fig_dir <- file.path(repo_root, "figures")
if (!file.exists(file.path(data_dir, "voom_log2_bcf.tsv"))) {
  stop("Run this script from a checkout that contains data/voom_log2_bcf.tsv")
}
if (!dir.exists(fig_dir)) {
  dir.create(fig_dir)
}

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(ggplot2)
  library(grid)
  library(openxlsx)
  library(patchwork)
  library(pheatmap)
})

treatment_order <- c("No treatment", "anti-CD21")
treatment_colors <- c("No treatment" = "skyblue", "anti-CD21" = "orange")
animal_colors <- c(
  "37172" = "#4A4A4A",
  "35684" = "#C5C5C5",
  "35661" = "#6B5B95",
  "37170" = "#8C510A"
)
gene_colors <- c("#2166AC", "#F7F7F7", "#B2182B")
heatmap_tx_colors <- c("anti-CD21" = "black", "No treatment" = "grey")

module_specs <- data.frame(
  key = c("immediate_early", "tfh", "nfkb", "ifn"),
  title = c("Immediate-early", "Tfh", "NF-kB", "IFN / ISG"),
  panel = c(
    "Figure 3D", "Figure 3E",
    "Supplemental Figure 11B", "Supplemental Figure 11C"
  ),
  ylabel = c(
    "Immediate-early score", "Tfh score", "NF-kB score", "IFN / ISG score"
  ),
  heatmap_sheet = c(
    "Fig 3D heatmap", "Fig 3E heatmap",
    "Fig S11B heatmap", "Fig S11C heatmap"
  ),
  score_sheet = c(
    "Fig 3D scores", "Fig 3E scores",
    "Fig S11B scores", "Fig S11C scores"
  ),
  stem = c(
    "figure3D_immediate_early", "figure3E_tfh",
    "figureS11B_nfkb", "figureS11C_ifn_isg"
  ),
  stats_module = c("Immediate-early", "Tfh module", "NF-kB", "IFN / ISG"),
  stringsAsFactors = FALSE
)

#' Require column names and return the data frame unchanged.
#'
#' @param frame Data frame to check.
#' @param columns Character vector of required column names.
#' @param name Label used in the error.
#' @return The input data frame.
RequireColumns <- function(frame, columns, name) {
  missing <- setdiff(columns, colnames(frame))
  if (length(missing) > 0L) {
    stop(name, " is missing columns: ", paste(missing, collapse = ", "))
  }
  frame
}

#' Row z-scores, clustered heatmap matrix, and per-ROI mean scores.
#'
#' Z-scores use \code{scale}: subtract the gene mean across ROIs and divide by
#' the sample standard deviation. Non-finite z-scores are omitted from the
#' mean and stored as 0 in the heatmap matrix.
#'
#' @param expression Numeric genes-by-samples matrix.
#' @return List with \code{heat} (clustered genes by samples) and \code{scores}.
ModuleMatrices <- function(expression) {
  if (is.null(expression) || length(dim(expression)) != 2L) {
    stop("expression must be a genes-by-samples matrix")
  }
  if (ncol(expression) < 2L || nrow(expression) < 2L) {
    stop("expression needs at least 2 genes and 2 samples")
  }
  values <- as.matrix(expression)
  storage.mode(values) <- "double"
  z <- suppressWarnings(t(scale(t(values))))
  dimnames(z) <- dimnames(values)
  z[!is.finite(z)] <- NA_real_
  if (all(!is.finite(z))) {
    stop("Row z-score matrix has no finite values")
  }
  scores <- colMeans(z, na.rm = TRUE)
  names(scores) <- colnames(values)
  heat <- z
  heat[!is.finite(heat)] <- 0
  if (nrow(heat) >= 2L) {
    row_order <- stats::hclust(stats::dist(heat), method = "complete")$order
    heat <- heat[row_order, , drop = FALSE]
  }
  list(heat = heat, scores = scores)
}

#' ComplexHeatmap of a pre-clustered z-score matrix, split by treatment.
#'
#' @param heat Genes-by-samples numeric matrix.
#' @param metadata Sample metadata aligned to the full ROI table.
#' @param color_limit Optional symmetric color limit. Defaults to max abs value.
#' @param show_gene_labels Draw gene symbols.
#' @param legend_at Optional legend breaks. Defaults to the color limits.
#' @return A ComplexHeatmap Heatmap object.
BuildHeatmap <- function(
    heat,
    metadata,
    color_limit = NULL,
    show_gene_labels = TRUE,
    legend_at = NULL
) {
  sample_info <- metadata[match(colnames(heat), metadata$sample), , drop = FALSE]
  if (anyNA(sample_info$treatment)) {
    stop("Heatmap columns are missing treatment metadata")
  }
  group <- factor(sample_info$treatment, levels = treatment_order)
  if (is.null(color_limit)) {
    color_limit <- max(abs(heat), na.rm = TRUE)
  }
  if (!is.finite(color_limit) || color_limit == 0) {
    color_limit <- 1
  }
  row_fontsize <- if (nrow(heat) > 20L) 8 else 10
  if (is.null(legend_at)) {
    legend_at <- c(-color_limit, 0, color_limit)
  }
  legend_labels <- formatC(legend_at, format = "fg", digits = 2)
  ComplexHeatmap::Heatmap(
    heat,
    name = "gene z-score",
    col = circlize::colorRamp2(
      c(-color_limit, 0, color_limit),
      gene_colors
    ),
    column_split = group,
    cluster_columns = FALSE,
    cluster_column_slices = FALSE,
    cluster_rows = FALSE,
    show_column_names = FALSE,
    show_row_names = isTRUE(show_gene_labels),
    row_names_gp = grid::gpar(fontsize = row_fontsize),
    column_title = treatment_order,
    column_title_gp = grid::gpar(fontsize = 10),
    column_gap = grid::unit(2, "mm"),
    border = TRUE,
    use_raster = FALSE,
    heatmap_legend_param = list(
      title = "gene z-score",
      at = legend_at,
      labels = legend_labels
    )
  )
}

#' Box and point plot of per-ROI module scores.
#'
#' @param scores Named numeric vector of ROI scores.
#' @param metadata Sample metadata.
#' @param title Plot title.
#' @param ylabel Y-axis label.
#' @param holm_p Holm-adjusted p-value drawn on the bracket.
#' @return A ggplot object.
BuildScorePlot <- function(scores, metadata, title, ylabel, holm_p) {
  if (is.null(names(scores)) || any(!nzchar(names(scores)))) {
    stop("scores must be named by ROI")
  }
  plot_data <- metadata[match(names(scores), metadata$sample), , drop = FALSE]
  plot_data$score <- as.numeric(scores)
  plot_data$treatment <- factor(plot_data$treatment, levels = treatment_order)
  plot_data$animal <- factor(plot_data$animal, levels = names(animal_colors))
  if (anyNA(plot_data$treatment) || anyNA(plot_data$score)) {
    stop("Score plot is missing treatment or score values")
  }
  score_min <- min(plot_data$score)
  score_max <- max(plot_data$score)
  span <- max(score_max - score_min, 1e-6)
  bracket_y <- score_max + 0.14 * span
  bracket <- data.frame(
    x = c(1, 1, 2, 2),
    y = c(
      bracket_y - 0.03 * span, bracket_y,
      bracket_y, bracket_y - 0.03 * span
    )
  )
  ggplot2::ggplot(plot_data, ggplot2::aes(x = treatment, y = score)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey55", linewidth = 0.3) +
    ggplot2::geom_boxplot(
      ggplot2::aes(fill = treatment),
      width = 0.55,
      outlier.shape = NA,
      color = "grey20",
      linewidth = 0.4,
      alpha = 0.85
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = animal),
      position = ggplot2::position_jitter(width = 0.12, height = 0, seed = 1),
      size = 2
    ) +
    ggplot2::geom_path(
      data = bracket,
      ggplot2::aes(x = x, y = y),
      inherit.aes = FALSE,
      linewidth = 0.4
    ) +
    ggplot2::annotate(
      "text",
      x = 1.5,
      y = bracket_y + 0.04 * span,
      label = paste0("Holm P = ", formatC(holm_p, format = "e", digits = 1)),
      size = 3.4
    ) +
    ggplot2::scale_fill_manual(values = treatment_colors, drop = FALSE) +
    ggplot2::scale_color_manual(values = animal_colors, drop = FALSE) +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::coord_cartesian(
      ylim = c(score_min - 0.1 * span, bracket_y + 0.22 * span),
      clip = "off"
    ) +
    ggplot2::guides(fill = "none") +
    ggplot2::labs(title = title, x = NULL, y = ylabel, color = "Animal") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(size = 12, face = "plain"),
      plot.margin = ggplot2::margin(8, 8, 8, 8)
    )
}

#' Record a ComplexHeatmap as a grob.
#'
#' @param heat_plot Undrawn Heatmap.
#' @param width Device width in inches.
#' @param height Device height in inches.
#' @param title Optional title drawn above the treatment slices.
#' @return A grob.
HeatmapGrob <- function(heat_plot, width, height, title = NULL) {
  grid::grid.grabExpr(
    {
      if (is.null(title)) {
        ComplexHeatmap::draw(heat_plot, newpage = TRUE)
      } else {
        ComplexHeatmap::draw(heat_plot, newpage = TRUE, column_title = title)
      }
    },
    width = width,
    height = height
  )
}

#' Write one heatmap sheet with annotation rows above the genes.
#'
#' @param workbook openxlsx workbook.
#' @param sheet Sheet name, at most 31 characters.
#' @param heat Genes-by-samples z-score matrix in figure order.
#' @param metadata Sample metadata.
#' @param annotation_columns Metadata columns written above the gene rows.
#' @return The sheet name, invisibly.
WriteHeatmapSheet <- function(
    workbook,
    sheet,
    heat,
    metadata,
    annotation_columns = c("treatment", "animal")
) {
  if (nchar(sheet) > 31L) {
    stop("Sheet name exceeds 31 characters: ", sheet)
  }
  annotation_labels <- c(treatment = "Treatment", animal = "Animal")
  unknown <- setdiff(annotation_columns, names(annotation_labels))
  if (length(unknown) > 0L) {
    stop("Unknown heatmap annotation columns: ", paste(unknown, collapse = ", "))
  }
  sample_info <- metadata[match(colnames(heat), metadata$sample), , drop = FALSE]
  if (anyNA(sample_info$sample)) {
    stop(sheet, " heatmap columns are missing from metadata")
  }
  openxlsx::addWorksheet(workbook, sheet)
  openxlsx::writeData(
    workbook, sheet,
    data.frame(t(c("gene", colnames(heat))), stringsAsFactors = FALSE),
    colNames = FALSE
  )
  for (annotation_index in seq_along(annotation_columns)) {
    column <- annotation_columns[[annotation_index]]
    openxlsx::writeData(
      workbook, sheet,
      data.frame(
        t(c(annotation_labels[[column]], sample_info[[column]])),
        stringsAsFactors = FALSE
      ),
      startRow = annotation_index + 1L,
      colNames = FALSE
    )
  }
  openxlsx::writeData(
    workbook, sheet,
    data.frame(gene = rownames(heat), heat, check.names = FALSE),
    startRow = length(annotation_columns) + 2L,
    colNames = FALSE
  )
  invisible(sheet)
}

#' Save a ggplot or patchwork object to png and pdf.
#'
#' @param plot Plot object accepted by \code{ggsave}.
#' @param stem File stem under figures/.
#' @param width Width in inches.
#' @param height Height in inches.
#' @return The png and pdf paths, invisibly.
SavePlot <- function(plot, stem, width, height) {
  png_path <- file.path(fig_dir, paste0(stem, ".png"))
  pdf_path <- file.path(fig_dir, paste0(stem, ".pdf"))
  ggplot2::ggsave(
    png_path, plot = plot, width = width, height = height, dpi = 300, bg = "white"
  )
  ggplot2::ggsave(
    pdf_path, plot = plot, width = width, height = height,
    device = grDevices::cairo_pdf, bg = "white"
  )
  message("Wrote ", basename(png_path), " and ", basename(pdf_path))
  invisible(c(png_path, pdf_path))
}

#' Row z-scores of a genes-by-samples matrix.
#'
#' Uses \code{scale}: subtract the gene mean and divide by the sample standard
#' deviation. This is the matrix \code{pheatmap} clusters when \code{scale = "row"}.
#'
#' @param expression Numeric genes-by-samples matrix.
#' @param name Label used in the error.
#' @return Numeric matrix of the same dimensions.
RowZScores <- function(expression, name) {
  values <- as.matrix(expression)
  storage.mode(values) <- "double"
  if (nrow(values) < 2L || ncol(values) < 2L) {
    stop(name, " needs at least 2 genes and 2 samples")
  }
  z <- suppressWarnings(t(scale(t(values))))
  dimnames(z) <- dimnames(values)
  if (any(!is.finite(z))) {
    stop(name, " row z-scores are not all finite")
  }
  z
}

#' Dark-blue to dark-red palette for Figure 3A.
#'
#' \code{pheatmap} requires one fewer color than the number of breaks.
#'
#' @param breaks Numeric color breaks.
#' @return Character vector of colors.
HeatmapPalette <- function(breaks) {
  grDevices::colorRampPalette(c(
    "darkblue", "mediumblue", "dodgerblue", "white",
    "orange", "red", "darkred"
  ))(length(breaks) - 1L)
}

#' Treatment annotation for the Figure 3A column bar.
#'
#' @param metadata Sample metadata with \code{sample} and \code{treatment}.
#' @return A one-column data frame aligned to \code{metadata$sample}.
TreatmentAnnotation <- function(metadata) {
  if (!all(metadata$treatment %in% c("No treatment", "anti-CD21"))) {
    stop("Treatment labels must be 'No treatment' or 'anti-CD21'")
  }
  data.frame(
    Tx = factor(metadata$treatment, levels = c("No treatment", "anti-CD21")),
    row.names = metadata$sample
  )
}

#' Cluster a row z-score matrix and save the Figure 3A heatmap.
#'
#' \code{pheatmap} opens a null PDF and closes the current device, so the
#' gtable is drawn here onto the png and pdf devices. Color is clipped at
#' +/- 3. Gene and ROI order match the dendrogram. ROI ids are omitted
#' because they do not fit; the source-data sheet lists them in this order.
#'
#' @param expression Genes-by-samples voom log2 matrix.
#' @param annotation Column annotation aligned to \code{expression}.
#' @param annotation_colors Named list of annotation colors.
#' @param main Heatmap title.
#' @param stem File stem under figures/.
#' @return List with \code{heat} (clustered z-scores).
DrawBcfHeatmap <- function(
    expression,
    annotation,
    annotation_colors,
    main,
    stem
) {
  if (!is.data.frame(annotation)) {
    stop("annotation must be a data frame")
  }
  if (!identical(rownames(annotation), colnames(expression))) {
    stop(stem, " annotation rows do not match expression columns")
  }
  missing_colors <- setdiff(colnames(annotation), names(annotation_colors))
  if (length(missing_colors) > 0L) {
    stop(
      stem, " annotation is missing colors for: ",
      paste(missing_colors, collapse = ", ")
    )
  }
  breaks <- seq(-3, 3, by = 0.1)
  colors <- HeatmapPalette(breaks)
  z <- RowZScores(expression, stem)
  built <- pheatmap::pheatmap(
    z,
    scale = "none",
    color = colors,
    breaks = breaks,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = FALSE,
    show_colnames = FALSE,
    border_color = "black",
    annotation_col = annotation,
    annotation_colors = annotation_colors,
    main = main,
    silent = TRUE
  )
  if (is.null(built$tree_row$order) || is.null(built$tree_col$order)) {
    stop(stem, " pheatmap did not return a cluster order")
  }
  ordered <- z[built$tree_row$order, built$tree_col$order, drop = FALSE]
  for (extension in c("png", "pdf")) {
    path <- file.path(fig_dir, paste0(stem, ".", extension))
    if (extension == "png") {
      grDevices::png(
        path, width = 8, height = 10, units = "in", res = 900, type = "cairo"
      )
    } else {
      grDevices::cairo_pdf(path, width = 8, height = 10)
    }
    grid::grid.newpage()
    grid::grid.draw(built$gtable)
    grDevices::dev.off()
    written <- file.info(path)$size
    if (is.na(written) || written < 1000) {
      stop("Failed to write ", path)
    }
    message("Wrote ", basename(path))
  }
  list(heat = ordered)
}

#' Figure 3B volcano: grey points, chocolate1 higher in anti-CD21, steelblue lower.
#'
#' A point is colored when Benjamini-Hochberg adjusted P < 0.05 and absolute
#' log2 fold change is greater than 0.58. Positive log2 fold change is higher
#' in anti-CD21. The y axis is -log10 of the adjusted P.
#'
#' @param de_table Differential expression table.
#' @param stem File stem under figures/.
#' @param main Plot title.
#' @return The plotted table with a \code{point} column (\code{up}, \code{down}, or \code{ns}).
DrawBcfVolcano <- function(de_table, stem, main) {
  de_table <- RequireColumns(
    de_table,
    c("gene", "logFC", "P.Value", "adj.P.Val"),
    stem
  )
  if (
    any(!is.finite(de_table$logFC)) ||
      any(!is.finite(de_table$P.Value)) ||
      any(!is.finite(de_table$adj.P.Val))
  ) {
    stop(stem, " has non-finite logFC, P.Value, or adj.P.Val")
  }
  de_table$point <- ifelse(
    de_table$logFC > 0.58 & de_table$adj.P.Val < 0.05,
    "up",
    ifelse(de_table$logFC < -0.58 & de_table$adj.P.Val < 0.05, "down", "ns")
  )
  for (extension in c("png", "pdf")) {
    path <- file.path(fig_dir, paste0(stem, ".", extension))
    if (extension == "png") {
      grDevices::png(path, width = 6, height = 4, units = "in", res = 300, type = "cairo")
    } else {
      grDevices::cairo_pdf(path, width = 6, height = 4)
    }
    graphics::plot(
      de_table$logFC,
      -log10(de_table$adj.P.Val),
      pch = 21,
      col = "gray90",
      bg = "grey",
      xlab = "Log2 fold change (positive = anti-CD21)",
      ylab = "-log10(adj. P)",
      main = main
    )
    up <- de_table$point == "up"
    down <- de_table$point == "down"
    graphics::points(
      de_table$logFC[up],
      -log10(de_table$adj.P.Val[up]),
      pch = 21,
      col = "black",
      bg = "chocolate1"
    )
    graphics::points(
      de_table$logFC[down],
      -log10(de_table$adj.P.Val[down]),
      pch = 21,
      col = "black",
      bg = "steelblue"
    )
    graphics::abline(h = 1.3, v = c(-0.58, 0.58), lty = 2)
    grDevices::dev.off()
    message("Wrote ", basename(path))
  }
  de_table
}

#' Write the volcano coordinates and the color class of each gene.
#'
#' @param workbook openxlsx workbook.
#' @param sheet Sheet name, at most 31 characters.
#' @param volcano Table returned by \code{DrawBcfVolcano}.
#' @return The sheet name, invisibly.
WriteVolcanoSheet <- function(workbook, sheet, volcano) {
  if (nchar(sheet) > 31L) {
    stop("Sheet name exceeds 31 characters: ", sheet)
  }
  plotted <- data.frame(
    gene = volcano$gene,
    logFC = volcano$logFC,
    neg_log10_adj.P.Val = -log10(volcano$adj.P.Val),
    P.Value = volcano$P.Value,
    adj.P.Val = volcano$adj.P.Val,
    point = volcano$point,
    stringsAsFactors = FALSE
  )
  openxlsx::addWorksheet(workbook, sheet)
  openxlsx::writeData(workbook, sheet, plotted)
  invisible(sheet)
}

voom <- read.delim(
  file.path(data_dir, "voom_log2_bcf.tsv"),
  check.names = FALSE,
  row.names = 1
)
metadata <- RequireColumns(
  read.delim(file.path(data_dir, "sample_metadata.tsv"), check.names = FALSE),
  c("sample", "animal", "treatment"),
  "sample_metadata.tsv"
)
module_genes <- RequireColumns(
  read.delim(file.path(data_dir, "module_genes.tsv"), check.names = FALSE),
  c("module", "gene"),
  "module_genes.tsv"
)
limma <- RequireColumns(
  read.delim(file.path(data_dir, "bcf_limma_topTable.tsv"), check.names = FALSE),
  c("gene", "logFC", "P.Value", "adj.P.Val"),
  "bcf_limma_topTable.tsv"
)
statistics <- RequireColumns(
  read.delim(file.path(data_dir, "module_statistics.tsv"), check.names = FALSE),
  c("figure", "module", "n_genes", "t", "P.Value", "holm_p"),
  "module_statistics.tsv"
)

if (anyDuplicated(metadata$sample) > 0L) {
  stop("sample_metadata.tsv has duplicate ROI ids")
}
missing_samples <- setdiff(metadata$sample, colnames(voom))
if (length(missing_samples) > 0L) {
  stop(
    "Metadata ROIs absent from voom: ",
    paste(missing_samples, collapse = ", ")
  )
}
voom <- voom[, metadata$sample, drop = FALSE]
n_control <- sum(metadata$treatment == "No treatment")
if (
  n_control < 1L ||
    n_control >= nrow(metadata) ||
    any(metadata$treatment[seq_len(n_control)] != "No treatment") ||
    any(metadata$treatment[(n_control + 1L):nrow(metadata)] != "anti-CD21")
) {
  stop("ROIs must be ordered no treatment, then anti-CD21")
}

figure3a_genes <- limma$gene[
  limma$adj.P.Val < 0.05 & abs(limma$logFC) >= 0.58
]
if (length(figure3a_genes) < 2L) {
  stop("Figure 3A gene filter returned fewer than 2 genes")
}
if (anyDuplicated(figure3a_genes) > 0L) {
  stop("Figure 3A gene list contains duplicates")
}
missing_figure3a <- setdiff(figure3a_genes, rownames(voom))
if (length(missing_figure3a) > 0L) {
  stop(
    "Figure 3A genes absent from voom: ",
    paste(missing_figure3a, collapse = ", ")
  )
}

message("BCF ROIs: ", ncol(voom))
message("Figure 3A genes: ", length(figure3a_genes))

modules <- vector("list", nrow(module_specs))
names(modules) <- module_specs$key
for (row_index in seq_len(nrow(module_specs))) {
  spec <- module_specs[row_index, ]
  genes <- module_genes$gene[module_genes$module == spec$key]
  if (length(genes) < 2L) {
    stop(spec$key, " needs at least 2 genes")
  }
  missing_genes <- setdiff(genes, rownames(voom))
  if (length(missing_genes) > 0L) {
    stop(spec$key, " genes absent from voom: ", paste(missing_genes, collapse = ", "))
  }
  stat_row <- statistics[statistics$module == spec$stats_module, , drop = FALSE]
  if (nrow(stat_row) != 1L) {
    stop("Expected one statistics row for ", spec$stats_module)
  }
  matrices <- ModuleMatrices(voom[genes, , drop = FALSE])
  if (nrow(matrices$heat) != stat_row$n_genes[[1]]) {
    stop(
      spec$key, " gene count ", nrow(matrices$heat),
      " != ", stat_row$n_genes[[1]]
    )
  }
  modules[[spec$key]] <- list(
    heat = matrices$heat,
    scores = matrices$scores,
    stats = stat_row
  )
  message(
    spec$title, ": ", nrow(matrices$heat), " genes, score range ",
    sprintf("%.3f to %.3f", min(matrices$scores), max(matrices$scores))
  )
}

figure3a <- DrawBcfHeatmap(
  expression = voom[figure3a_genes, , drop = FALSE],
  annotation = TreatmentAnnotation(metadata),
  annotation_colors = list(Tx = heatmap_tx_colors),
  main = "Figure 3A",
  stem = "figure3A_bcf_heatmap"
)
bcf_volcano <- DrawBcfVolcano(
  limma,
  stem = "figure3B_bcf_volcano",
  main = "Figure 3B"
)

for (row_index in seq_len(nrow(module_specs))) {
  spec <- module_specs[row_index, ]
  module <- modules[[spec$key]]
  figure_height <- max(4.6, 0.28 * nrow(module$heat) + 2)
  figure_width <- 11
  heat_width <- figure_width * 1.7 / 2.7
  heat_plot <- BuildHeatmap(
    heat = module$heat,
    metadata = metadata,
    show_gene_labels = TRUE
  )
  score_plot <- BuildScorePlot(
    scores = module$scores,
    metadata = metadata,
    title = paste(spec$title, "score"),
    ylabel = spec$ylabel,
    holm_p = module$stats$holm_p[[1]]
  )
  combined <- patchwork::wrap_elements(
    HeatmapGrob(heat_plot, width = heat_width, height = figure_height)
  ) +
    score_plot +
    patchwork::plot_layout(widths = c(1.7, 1)) +
    patchwork::plot_annotation(title = paste(spec$panel, spec$title))
  SavePlot(combined, spec$stem, width = figure_width, height = figure_height)
}

statistics$figure <- module_specs$panel[
  match(statistics$module, module_specs$stats_module)
]
if (anyNA(statistics$figure)) {
  stop("module_statistics.tsv has a module that is not in the figure list")
}

notes <- data.frame(
  sheet = c(
    "Notes",
    "Figure 3A",
    "Figure 3B",
    as.vector(rbind(module_specs$heatmap_sheet, module_specs$score_sheet)),
    "Module statistics"
  ),
  contents = c(
    paste(
      "BCF ROIs only. Heatmap values are row z-scores of voom log2 expression",
      "(mean subtracted, sample standard deviation, as in R scale).",
      "Score values are the mean of those z-scores across module genes.",
      "Figure 3D, 3E, 11B, and 11C columns are ordered no treatment,",
      "then anti-CD21. Figure 3A clusters genes and ROIs.",
      "Figure 3B point is up, down, or ns when Benjamini-Hochberg adjusted",
      "P < 0.05 and absolute log2 fold change > 0.58. Positive logFC is",
      "higher in anti-CD21."
    ),
    paste(
      "Figure 3A. Unclipped row z-scores for BCF genes with adjusted P < 0.05",
      "and absolute log2 fold change >= 0.58. Rows and columns are in the",
      "clustered order drawn in figures/figure3A_bcf_heatmap. Colors are",
      "clipped at +/- 3."
    ),
    paste(
      "Figure 3B. All genes in the BCF anti-CD21 versus no-treatment contrast.",
      "Positive logFC is higher in anti-CD21. neg_log10_adj.P.Val is the",
      "plotted y coordinate. point up is higher in anti-CD21. The y axis is",
      "-log10 of the Benjamini-Hochberg adjusted P."
    ),
    "Figure 3D. Row z-scores for the immediate-early gene heatmap.",
    "Figure 3D. Per-ROI immediate-early scores drawn in the score plot.",
    "Figure 3E. Row z-scores for the Tfh gene heatmap.",
    "Figure 3E. Per-ROI Tfh scores drawn in the score plot.",
    "Supplemental Figure 11B. Row z-scores for the NF-kB gene heatmap.",
    "Supplemental Figure 11B. Per-ROI NF-kB scores drawn in the score plot.",
    "Supplemental Figure 11C. Row z-scores for the IFN/ISG gene heatmap.",
    "Supplemental Figure 11C. Per-ROI IFN/ISG scores drawn in the score plot.",
    paste(
      "Gene counts, limma t statistics, raw P values, and Holm-adjusted P",
      "values annotated on the four score plots. Holm adjustment is across",
      "the immediate-early, Tfh, NF-kB, and IFN/ISG modules."
    )
  ),
  stringsAsFactors = FALSE
)
if (any(nchar(notes$sheet) > 31L)) {
  stop("A Notes sheet label exceeds 31 characters")
}

workbook_path <- file.path(data_dir, "Source_data.xlsx")
workbook <- openxlsx::createWorkbook()
openxlsx::addWorksheet(workbook, "Notes")
openxlsx::writeData(workbook, "Notes", notes)
openxlsx::setColWidths(workbook, "Notes", cols = 1:2, widths = c(28, 88))
openxlsx::addStyle(
  workbook, "Notes",
  openxlsx::createStyle(wrapText = TRUE, valign = "top"),
  rows = seq_len(nrow(notes) + 1L),
  cols = 1:2,
  gridExpand = TRUE,
  stack = TRUE
)
WriteHeatmapSheet(workbook, "Figure 3A", figure3a$heat, metadata)
WriteVolcanoSheet(workbook, "Figure 3B", bcf_volcano)
for (row_index in seq_len(nrow(module_specs))) {
  spec <- module_specs[row_index, ]
  module <- modules[[spec$key]]
  WriteHeatmapSheet(workbook, spec$heatmap_sheet, module$heat, metadata)
  score_table <- metadata[match(names(module$scores), metadata$sample), , drop = FALSE]
  score_table$score <- as.numeric(module$scores)
  openxlsx::addWorksheet(workbook, spec$score_sheet)
  openxlsx::writeData(
    workbook,
    spec$score_sheet,
    score_table[, c("sample", "animal", "treatment", "score")]
  )
}
openxlsx::addWorksheet(workbook, "Module statistics")
openxlsx::writeData(workbook, "Module statistics", statistics)
openxlsx::saveWorkbook(workbook, workbook_path, overwrite = TRUE)
message("Wrote ", workbook_path)
stray_plot <- file.path(repo_root, "Rplots.pdf")
if (file.exists(stray_plot)) {
  unlink(stray_plot)
}
message("Sheets: ", paste(openxlsx::getSheetNames(workbook_path), collapse = ", "))
