#################### Medicago (medtr) Mac version ####################

# Load required packages and install them if necessary
packages <- c("tximport", "edgeR", "statmod", "dplyr", "tibble", "tidyr", "readr", "ggplot2", "pheatmap")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "http://cran.us.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

# Install Bioconductor packages if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("tximport", ask = FALSE)

# Set input paths
base_dir <- "/Users/hyunohlee/Analysis/SNF_genes_phylogeny/Chamaecrista_RNA-seq/Redo_06-18-25/medtr/Plant physiology/salmon_output"
tx2gene_path <- "/Users/hyunohlee/Analysis/SNF_genes_phylogeny/Chamaecrista_RNA-seq/Redo_06-18-25/medtr/Plant physiology/edgeR/tx2gene.csv"

# Define sample metadata
sample_info <- data.frame(
  sample = c("K-FN-Root_1", "K-FN-Root_2", "K-FN-Root_3", "K-FN-Shoot_1", "K-FN-Shoot_2", "K-FN-Shoot_3",
             "Mt4wkNod_1", "Mt4wkNod_2", "Mt4wkNod_3", "MtNod0dpi_1", "MtNod0dpi_2", "MtNod0dpi_3",
             "MtNod10dpi_1", "MtNod10dpi_2", "MtNod10dpi_3", "MtNod14dpi_1", "MtNod14dpi_2", "MtNod14dpi_3",
             "MtNod14dpi_12h_1", "MtNod14dpi_12h_2", "MtNod14dpi_12h_3", "MtNod14dpi_48h_1", "MtNod14dpi_48h_2", "MtNod14dpi_48h_3",
             "MtNod4dpi_1", "MtNod4dpi_2", "MtNod4dpi_3", "N-FN-Root_1", "N-FN-Root_2", "N-FN-Root_3",
             "N-FN-Shoot_1", "N-FN-Shoot_2", "N-FN-Shoot_3", "P-FN-Root_1", "P-FN-Root_2", "P-FN-Root_3",
             "P-FN-Shoot_1", "P-FN-Shoot_2", "P-FN-Shoot_3", "S-FN-Root_1", "S-FN-Root_2", "S-FN-Root_3",
             "S-FN-Shoot_1", "S-FN-Shoot_2", "S-FN-Shoot_3"),
  group = c(rep("K_FN_Root", 3), rep("K_FN_Shoot", 3), rep("4wkNodule", 3), rep("Nodule0dpi", 3), rep("Nodule10dpi", 3), rep("Nodule14dpi", 3),
            rep("Nodule14dpi_12h", 3), rep("Nodule14dpi_48h", 3), rep("Nodule4dpi", 3), rep("N_FN_Root", 3), rep("N_FN_Shoot", 3),
            rep("P_FN_Root", 3), rep("P_FN_Shoot", 3), rep("S_FN_Root", 3), rep("S_FN_Shoot", 3)),
  stringsAsFactors = FALSE
)

# Locate quant.sf files
files <- file.path(base_dir, sample_info$sample, "quant.sf")
names(files) <- sample_info$sample

# Check file existence
missing_files <- files[!file.exists(files)]
if (length(missing_files) > 0) {
  stop("Missing quant.sf files:
", paste(names(missing_files), collapse = ", "))
}

# Load tx2gene mapping file
if (!file.exists(tx2gene_path)) stop("tx2gene.csv was not found. Check the path.")
tx2gene <- read_csv(tx2gene_path)

# Import transcript counts using tximport
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "no")

# Run edgeR
group <- factor(sample_info$group)
dge <- DGEList(counts = txi$counts, group = group)
dge <- calcNormFactors(dge, method = "TMM")
design <- model.matrix(~group)
dge <- estimateDisp(dge, design, robust = TRUE)
fit <- glmFit(dge, design)
lrt <- glmLRT(fit)

# Save differential expression results
write_csv(topTags(lrt, n = Inf)$table, "1.medtr_DEG_results_edgeR.csv")
cat("
Analysis completed. Results were saved to '1.medtr_DEG_results_edgeR.csv'.
")

# Calculate group means and standard errors
count_data <- as.data.frame(txi$counts)
group_mapping <- setNames(sample_info$group, sample_info$sample)

group_means <- count_data %>%
  rownames_to_column("gene_id") %>%
  pivot_longer(-gene_id, names_to = "sample", values_to = "count") %>%
  mutate(group = group_mapping[sample]) %>%
  group_by(gene_id, group) %>%
  summarise(
    mean_expr = mean(count, na.rm = TRUE),
    se_expr = sd(count, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = group, values_from = c(mean_expr, se_expr))

# Calculate relative expression percentages
relative_expr <- group_means %>%
  mutate(total_expr = rowSums(select(., starts_with("mean_expr_")), na.rm = TRUE)) %>%
  mutate(across(starts_with("mean_expr_"), ~ (. / total_expr) * 100, .names = "RELATIVE_{.col}")) %>%
  mutate(across(starts_with("se_expr_"), ~ (. / total_expr) * 100, .names = "RELATIVE_SE_{.col}")) %>%
  select(-total_expr)

# Prepare edgeR result table
edgeR_results <- as.data.frame(topTags(lrt, n = Inf)$table)
edgeR_results$gene_id <- rownames(edgeR_results)
edgeR_results <- edgeR_results %>%
  select(gene_id, logFC, FDR, PValue)

# Merge expression summaries with edgeR results and annotation
final_table <- left_join(relative_expr, edgeR_results, by = "gene_id")
final_table <- left_join(final_table, tx2gene, by = "gene_id")

# Clean duplicated columns after merging
final_table <- final_table %>%
  rename_with(~ gsub("_x$", "", .), ends_with("_x")) %>%
  select(-ends_with("_y"))

# Save the final table
write_csv(final_table, "2.GeneExpression_GroupStats.csv")
cat("
Group-level expression summary with edgeR statistics was saved to '2.GeneExpression_GroupStats.csv'.
")

# Retain only transcript IDs ending in .1
final_table_isoform1 <- final_table %>%
  filter(grepl("\.1$", transcript_id))
write_csv(final_table_isoform1, "2.GeneExpression_GroupStats_isoform1.csv")
cat("
Isoform-1-only table was saved to '2.GeneExpression_GroupStats_isoform1.csv'.
")

# Filter significant DE genes (FDR < 0.05 and absolute logFC > 1)
DE_genes <- final_table %>% filter(FDR < 0.05 & abs(logFC) > 1)
write_csv(DE_genes, "3.Significant_DEGs.csv")
cat("
Significant DE genes were saved to '3.Significant_DEGs.csv'.
")

# Save the isoform-1 subset of significant DE genes
DE_genes_isoform1 <- DE_genes %>%
  filter(grepl("\.1$", transcript_id))
write_csv(DE_genes_isoform1, "3.Significant_DEGs_isoform1.csv")
cat("
Isoform-1-only significant DE genes were saved to '3.Significant_DEGs_isoform1.csv'.
")
