# Purpose: Generate relative-expression heatmaps for six selected species.
# Input: Excel table of expression summaries and a directory of gene tree files.
# Output: Tree-ordered heatmaps saved in PNG format.

library(readxl)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(ape)
library(tidyr)
library(dendextend)
library(tibble)

# File paths
file_path <- "/Users/hyunohlee/Analysis/SNF_genes_phylogeny/Chamaecrista_RNA-seq/Redo_06-18-25/Counts/snf_families_counts_2025-07-01.xlsx"
tree_dir <- "/Users/hyunohlee/Analysis/SNF_genes_phylogeny/Chamaecrista_RNA-seq/Redo_06-18-25/Trees"
output_dir <- "Relative_heatmaps_PNG"
dir.create(output_dir, showWarnings = FALSE)

# Load input data
df <- read_excel(file_path, sheet = "genes in tree order")

# Relative-expression columns for each species
expr_cols <- list(
  medtr = c("Medicago nodule; relative mean", "Medicago root; relative mean", "Medicago shoot tip relative percent"),
  glyma = c("Glycine nodule relative percent", "Glycine root relative percent", "Glycine shoot tip relative percent"),
  phavu = c("Phaseolus nodule relative mean", "Phaseolus root relative mean", "Phaseolus shoot relative mean"),
  chafa = c("Chamaecrista nodule sum of rates", "Chamaecrista root; sum of rates", "Chamaecrista shoot; sum of rates")
)

# Filter the target species and build three-tissue expression tables
# for each gene.
df_filtered <- df %>%
  filter(prefix %in% c("medtr", "glyma", "phavu", "chafa", "vitvi", "arath")) %>%
  rowwise() %>%
  mutate(
    gene_id_clean = `gene ID`,
    Nodule = ifelse(prefix %in% names(expr_cols), as.numeric(get(expr_cols[[prefix]][1])), NA_real_),
    Root = ifelse(prefix %in% names(expr_cols), as.numeric(get(expr_cols[[prefix]][2])), NA_real_),
    Shoot = ifelse(prefix %in% names(expr_cols), as.numeric(get(expr_cols[[prefix]][3])), NA_real_)
  ) %>%
  ungroup()

# Load all tree files
tree_files <- list.files(tree_dir, pattern = "\.rt$", full.names = TRUE)
gene_trees <- lapply(tree_files, read.tree)

# Generate one heatmap per tree
for (i in seq_along(gene_trees)) {
  tree <- gene_trees[[i]]
  tree_file <- tree_files[[i]]
  clean_name <- gsub("\.rt$", "", basename(tree_file))

  common_genes <- intersect(tree$tip.label, df_filtered$gene_id_clean)

  if (length(common_genes) < 2) next

  tree <- keep.tip(tree, common_genes)

  gene_order <- tree$tip.label[tree$edge[tree$edge[, 2] <= length(tree$tip.label), 2]]

  expr_df <- df_filtered %>%
    filter(gene_id_clean %in% common_genes) %>%
    select(gene_id_clean, Nodule, Root, Shoot) %>%
    column_to_rownames("gene_id_clean")

  expr_df <- expr_df[gene_order, , drop = FALSE]
  expr_df <- expr_df %>% select(where(~ !all(is.na(.))))
  expr_mat <- as.matrix(expr_df)

  # Build a dendrogram based on tree cophenetic distances and rotate it
  # to match the tree tip order.
  dist_matrix <- cophenetic(tree)
  hc <- hclust(as.dist(dist_matrix))
  dend <- as.dendrogram(hc)
  dend <- dendextend::rotate(dend, order = gene_order)

  col_fun <- colorRamp2(seq(0, 100, length.out = 5), c("white", "lightyellow", "yellow", "orange", "red"))

  prefix_labels_ordered <- ifelse(grepl("^PAP\.", gene_order), "PAP",
                                  ifelse(grepl("^CAE\.", gene_order), "CAE", "OUT"))

  fontsize_row <- ifelse(nrow(expr_mat) > 60, 5, ifelse(nrow(expr_mat) > 40, 6, ifelse(nrow(expr_mat) > 20, 8, 10)))

  out_file <- file.path(output_dir, paste0("heatmap_", clean_name, ".png"))

  tryCatch({
    png(out_file, width = 12, height = 15, units = "in", res = 300)

    main_heatmap <- Heatmap(expr_mat,
                            name = "Expression",
                            col = col_fun,
                            cluster_rows = dend,
                            cluster_columns = FALSE,
                            show_row_names = FALSE,
                            na_col = "white",
                            column_title = paste("Relative Expression Heatmap -", clean_name),
                            rect_gp = gpar(col = "black", lwd = 0.5),
                            
                            row_dend_width = unit(12, "cm"))

    prefix_annotation <- Heatmap(matrix(prefix_labels_ordered, ncol = 1),
                                 name = "Subfam",
                                 col = c("PAP" = "blue", "CAE" = "red", "OUT" = "gray"),
                                 width = unit(2.5, "mm"),
                                 show_row_names = FALSE,
                                 cluster_rows = dend,
                                 cluster_columns = FALSE,
                                 rect_gp = gpar(col = "black", lwd = 0.5))

    gene_annotation <- rowAnnotation(Gene = anno_text(gene_order, gp = gpar(fontsize = fontsize_row)))

    ht_list <- main_heatmap + prefix_annotation + gene_annotation

    draw(ht_list, heatmap_legend_side = "bottom", padding = unit(c(2, 2, 2, 2), "mm"), auto_adjust = TRUE)

    dev.off()

    cat("Saved:", out_file, "\n\n")
  }, error = function(e) {
    dev.off()
    cat("Error while processing:", clean_name, "\n", conditionMessage(e), "\n")
  })
}

cat("All PNG heatmaps were generated successfully.\n")
