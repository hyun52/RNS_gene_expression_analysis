# Purpose: Generate PCA and UMAP plots from superfamily-level expression summaries.
# Input: Excel file containing expression columns and bias classifications.
# Output: PCA and UMAP plots in PNG, PDF, and SVG formats.

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(umap)
library(patchwork)
library(scales)
# If needed for stable SVG output, load svglite explicitly.
# library(svglite)

# Load input data
df <- read_excel(
  "/Users/hyunohlee/Analysis/SNF_genes_phylogeny/Chamaecrista_RNA-seq/Redo_06-18-25/Counts/snf_families_counts_2025-07-02.xlsx",
  sheet = "genes in tree order"
)

# Expression columns used for PCA and UMAP
expr_cols <- c(
  "Medicago nodule; relative mean", "Medicago root; relative mean", "Medicago shoot tip relative percent",
  "Glycine nodule relative percent", "Glycine root relative percent", "Glycine shoot tip relative percent",
  "Phaseolus nodule relative mean", "Phaseolus root relative mean", "Phaseolus shoot relative mean",
  "Chamaecrista nodule sum of rates", "Chamaecrista root; sum of rates", "Chamaecrista shoot; sum of rates"
)

# Summarize at the SFam ID level and remove entries lacking expression values
df_sf <- df %>%
  group_by(`SFam ID`) %>%
  summarise(
    across(all_of(expr_cols), ~ mean(as.numeric(.), na.rm = TRUE)),
    medtr_bias   = first(na.omit(medtr_bias)),
    glyma_bias   = first(na.omit(glyma_bias)),
    phavu_bias   = first(na.omit(phavu_bias)),
    chafa_bias   = first(na.omit(chafa_bias)),
    overall_bias = first(na.omit(overall_bias))
  ) %>%
  ungroup() %>%
  drop_na(all_of(expr_cols))

# Expression matrix used for dimensionality reduction
df_expr <- df_sf %>% select(all_of(expr_cols))

# PCA
pca_res <- prcomp(df_expr, center = TRUE, scale. = TRUE)
pca_df  <- cbind(as.data.frame(pca_res$x), df_sf)

# UMAP
set.seed(123)
umap_res <- umap(df_expr)
umap_df  <- cbind(as.data.frame(umap_res$layout), df_sf)
colnames(umap_df)[1:2] <- c("UMAP1", "UMAP2")

# Color and shape settings
bias_palette <- c(
  "nodule"  = "#D55E00",
  "root"    = "#009E73",
  "shoot"   = "#0072B2",
  "neutral" = "#999999"
)

shape_mapping <- c(
  "medtr_bias" = 16,
  "glyma_bias" = 17,
  "phavu_bias" = 15,
  "chafa_bias" = 18
)

# Species-specific bias panels
pca_species_df <- pca_df %>%
  select(PC1, PC2, medtr_bias, glyma_bias, phavu_bias, chafa_bias) %>%
  pivot_longer(cols = c(medtr_bias, glyma_bias, phavu_bias, chafa_bias),
               names_to = "species", values_to = "bias")

umap_species_df <- umap_df %>%
  select(UMAP1, UMAP2, medtr_bias, glyma_bias, phavu_bias, chafa_bias) %>%
  pivot_longer(cols = c(medtr_bias, glyma_bias, phavu_bias, chafa_bias),
               names_to = "species", values_to = "bias")

pca_species_plot <- ggplot(pca_species_df, aes(PC1, PC2, color = bias)) +
  geom_point(size = 2) +
  facet_wrap(~ species) +
  theme_minimal() +
  scale_color_manual(values = bias_palette) +
  labs(title = "PCA: Species-specific Bias (SFam)", color = "Bias")

umap_species_plot <- ggplot(umap_species_df, aes(UMAP1, UMAP2, color = bias)) +
  geom_point(size = 2) +
  facet_wrap(~ species) +
  theme_minimal() +
  scale_color_manual(values = bias_palette) +
  labs(title = "UMAP: Species-specific Bias (SFam)", color = "Bias")

# Combined species-bias panel using shape distinctions
pca_species_single <- ggplot(pca_species_df, aes(PC1, PC2, color = bias, shape = species)) +
  geom_point(size = 2.5, alpha = 0.8) +
  theme_minimal() +
  scale_color_manual(values = bias_palette) +
  scale_shape_manual(values = shape_mapping) +
  labs(title = "PCA: Combined Species Bias (SFam)", color = "Bias", shape = "Species")

umap_species_single <- ggplot(umap_species_df, aes(UMAP1, UMAP2, color = bias, shape = species)) +
  geom_point(size = 2.5, alpha = 0.8) +
  theme_minimal() +
  scale_color_manual(values = bias_palette) +
  scale_shape_manual(values = shape_mapping) +
  labs(title = "UMAP: Combined Species Bias (SFam)", color = "Bias", shape = "Species")

# Overall bias plots
pca_overall_plot <- ggplot(pca_df, aes(PC1, PC2, color = overall_bias)) +
  geom_point(size = 2.5, alpha = 0.8) +
  theme_minimal() +
  scale_color_manual(values = bias_palette) +
  labs(title = "PCA: SFam Overall Bias", color = "Bias")

umap_overall_plot <- ggplot(umap_df, aes(UMAP1, UMAP2, color = overall_bias)) +
  geom_point(size = 2.5, alpha = 0.8) +
  theme_minimal() +
  scale_color_manual(values = bias_palette) +
  labs(title = "UMAP: SFam Overall Bias", color = "Bias")

# Save plots in multiple formats
plots <- list(
  PCA_species_4panels   = pca_species_plot,
  PCA_species_combined  = pca_species_single,
  PCA_overall_bias      = pca_overall_plot,
  UMAP_species_4panels  = umap_species_plot,
  UMAP_species_combined = umap_species_single,
  UMAP_overall_bias     = umap_overall_plot
)

outdir <- "Figures"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

for (name in names(plots)) {
  p <- plots[[name]]
  # High-resolution PNG
  ggsave(file.path(outdir, paste0(name, ".png")), p, width = 10, height = 8, dpi = 300)
  # Vector PDF
  ggsave(file.path(outdir, paste0(name, ".pdf")), p, width = 10, height = 8)
  # Vector SVG for editing
  ggsave(file.path(outdir, paste0(name, ".svg")), p, width = 10, height = 8)
}

message("Saved output to: ", normalizePath(outdir))
