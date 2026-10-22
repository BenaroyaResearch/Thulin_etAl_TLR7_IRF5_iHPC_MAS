# This code was generated using Claude Code, Anthropic (Opus 5.5) and carefully by a bioinformatician.

# Overview: This code generates the transcription factor motif similarity matrix in Fig. S2 B. 

# Run from the Motif-Sim-Matrix_FigS2B folder:
#   Rscript Code/FigS2B_motif_similarity_heatmap.R

suppressPackageStartupMessages({
  library(readr); library(readxl)
  library(ComplexHeatmap); library(circlize); library(grid)
})

data_dir <- "Data"
out_dir  <- "Output"
dir.create(out_dir, showWarnings = FALSE)

sim_file <- "CellTypeConsensus_TFmotif_SimilarityMatrix.csv"
cat_file <- "TF-Motif-Categories.xlsx"

# ── Reload ─────────────────────────────────────────────────────────────────
mat_raw  <- read_csv(file.path(data_dir, sim_file), show_col_types = FALSE)
mat_data <- as.matrix(mat_raw[, -1])
rownames(mat_data) <- mat_raw[[1]]
colnames(mat_data) <- mat_raw[[1]]

cats       <- read_excel(file.path(data_dir, cat_file))
group_col  <- cats[[1]]
motif_cols <- cats[, 2:11]

normalize <- function(x) gsub(" ", "", x)
motif_group <- list()
for (i in seq_len(nrow(cats))) {
  grp    <- group_col[i]
  motifs <- as.character(unlist(motif_cols[i, ]))
  motifs <- motifs[!is.na(motifs)]
  for (m in motifs) motif_group[[normalize(m)]] <- grp
}

all_motifs <- rownames(mat_data)
missing    <- all_motifs[!normalize(all_motifs) %in% names(motif_group)]
if (length(missing) > 0)
  stop("Motifs missing from ", cat_file, ": ", paste(missing, collapse = ", "))
group_vec  <- sapply(all_motifs, function(m) motif_group[[normalize(m)]])

# ── Ordering ───────────────────────────────────────────────────────────────
unique_groups <- unique(group_vec)
group_profiles <- sapply(unique_groups, function(grp) {
  members <- all_motifs[group_vec == grp]
  colMeans(mat_data[members, , drop = FALSE])
})
group_hc    <- hclust(dist(t(group_profiles), method = "euclidean"), method = "complete")
group_order <- unique_groups[group_hc$order]

final_order <- c()
for (grp in group_order) {
  members <- all_motifs[group_vec == grp]
  if (length(members) == 1) {
    final_order <- c(final_order, members)
  } else {
    hc <- hclust(as.dist(1 - mat_data[members, members]), method = "complete")
    final_order <- c(final_order, members[hc$order])
  }
}

group_blocks <- data.frame(group = group_order, start = NA_integer_, end = NA_integer_,
                            stringsAsFactors = FALSE)
for (i in seq_along(group_order)) {
  idx <- which(final_order %in% all_motifs[group_vec == group_order[i]])
  group_blocks$start[i] <- min(idx); group_blocks$end[i] <- max(idx)
}

mat_ordered <- mat_data[final_order, final_order]
mat_tri     <- mat_ordered
mat_tri[lower.tri(mat_tri, diag = TRUE)] <- NA

# ── Colours ────────────────────────────────────────────────────────────────
col_fun <- colorRamp2(c(-0.1, 0, 0.5, 1),
                      c("#2166AC", "#F7F7F7", "#F4A582", "#B2182B"))

grp_colours <- c(
  "#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00",
  "#A65628","#F781BF","#999999","#66C2A5","#FC8D62",
  "#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494",
  "#B3B3B3","#1B9E77","#D95F02","#7570B3","#E7298A",
  "#66A61E","#E6AB02","#A6761D","#666666","#8DD3C7",
  "#FFFFB3","#BEBADA","#FB8072","#80B1D3"
)
names(grp_colours) <- group_order

n <- nrow(mat_tri)
group_label_vec <- group_vec[final_order]

row_ann <- rowAnnotation(
  Group = group_label_vec,
  col   = list(Group = grp_colours),
  show_legend        = FALSE,
  annotation_name_gp = gpar(fontsize = 0),
  simple_anno_size   = unit(4, "mm"),
  annotation_label   = ""
)

ht <- Heatmap(
  mat_tri,
  name              = "Similarity",
  col               = col_fun,
  na_col            = "white",
  cluster_rows      = FALSE,
  cluster_columns   = FALSE,
  show_row_names    = TRUE,
  show_column_names = FALSE,
  row_names_side    = "left",
  row_names_gp      = gpar(fontsize = 5.5),
  row_names_max_width = unit(3.5, "cm"),
  rect_gp           = gpar(col = NA),
  heatmap_legend_param = list(
    title         = "Motif similarity",
    title_gp      = gpar(fontsize = 9, fontface = "bold"),
    labels_gp     = gpar(fontsize = 8),
    legend_width  = unit(6, "cm"),
    direction     = "horizontal",
    at            = c(-0.1, 0, 0.25, 0.5, 0.75, 1),
    labels        = c("-0.1", "0", "0.25", "0.5", "0.75", "1")
  ),
  right_annotation = row_ann,
  width  = unit(20, "cm"),
  height = unit(20, "cm")
)

# ── Label push-apart ───────────────────────────────────────────────────────
ideal_y <- sapply(seq_len(nrow(group_blocks)), function(i) {
  s <- group_blocks$start[i]; e <- group_blocks$end[i]
  1 - (s + e - 1) / (2 * n)
})
line_h <- 7.5 / (20 * 28.35)   # 7.5pt font on 20cm body
adj_y  <- ideal_y
for (i in 2:length(adj_y)) {
  if (adj_y[i-1] - adj_y[i] < line_h)
    adj_y[i] <- adj_y[i-1] - line_h
}
adj_y <- pmax(0.005, pmin(0.995, adj_y))

# ── Group label placement ──────────────────────────────────────────────────
# Labels sit to the right of the colour bar: annotation padding (1 mm) +
# bar (4 mm) + gap. Black text so pale group colours stay legible; the bar
# carries the colour.
lab_gp  <- gpar(fontsize = 7.5, fontface = "bold", col = "black")
lab_gap <- 1.5
lab_x   <- unit(1, "npc") + unit(1 + 4 + lab_gap, "mm")

pdf(NULL)
lab_w <- max(sapply(group_blocks$group, function(g)
  convertWidth(grobWidth(textGrob(g, gp = lab_gp)), "mm", valueOnly = TRUE)))
invisible(dev.off())

# draw() padding order is bottom, left, top, right; the right side holds labels
plot_padding <- unit(c(5, 5, 5, 5 + 1 + 4 + lab_gap + lab_w), "mm")

# Size the page to the drawn heatmap so there is no surrounding white space
pdf(NULL)
ht_drawn <- draw(ht, heatmap_legend_side = "bottom", padding = plot_padding)
page_w <- convertWidth(ComplexHeatmap:::width(ht_drawn), "inches", valueOnly = TRUE)
page_h <- convertHeight(ComplexHeatmap:::height(ht_drawn), "inches", valueOnly = TRUE)
invisible(dev.off())

# ── Render (PDF with editable text, plus SVG) ──────────────────────────────
render <- function() {
  draw(ht,
       heatmap_legend_side = "bottom",
       padding = plot_padding)

  for (i in seq_len(nrow(group_blocks))) {
    s <- group_blocks$start[i]; e <- group_blocks$end[i]
    decorate_heatmap_body("Similarity", {
      x0 <- (s-1)/n; x1 <- e/n; y0 <- 1-e/n; y1 <- 1-(s-1)/n
      grid.rect(x=unit(x0,"npc"), y=unit(y0,"npc"),
                width=unit(x1-x0,"npc"), height=unit(y1-y0,"npc"),
                just=c("left","bottom"),
                gp=gpar(col="grey20", lwd=1.2, fill=NA))
    })
  }

  for (i in seq_len(nrow(group_blocks))) {
    grp <- group_blocks$group[i]
    decorate_heatmap_body("Similarity", {
      grid.text(
        grp,
        x    = lab_x,
        y    = unit(adj_y[i], "npc"),
        just = c("left", "centre"),
        gp   = lab_gp
      )
    })
  }
}


pdf(file.path(out_dir, "TFmotif_similarity_heatmap.pdf"),
    width = page_w, height = page_h, family = "Helvetica", useDingbats = FALSE)
render()
dev.off()

svglite::svglite(file.path(out_dir, "TFmotif_similarity_heatmap.svg"),
                 width = page_w, height = page_h)
render()
dev.off()

write.csv(data.frame(motif = final_order, group = group_label_vec, row.names = NULL),
          file.path(out_dir, "TFmotif_order_and_groups.csv"), row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
cat("Wrote PDF, SVG, motif order CSV and sessionInfo to", out_dir, "\n")
