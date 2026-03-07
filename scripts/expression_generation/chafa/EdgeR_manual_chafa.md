#################### Chamaecrista (chafa) Mac version ####################

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
base_dir <- "/Users/hyunohlee/Analysis/SNF_genes_phylogeny/Chamaecrista_RNA-seq/Redo_06-18-25/chafa/salmon_output"
tx2gene_path <- "/Users/hyunohlee/Analysis/SNF_genes_phylogeny/Chamaecrista_RNA-seq/Redo_06-18-25/chafa/edgeR/tx2gene.csv"

# Define sample metadata
sample_info <- data.frame(
  sample = c(
    "nodule_and_root_1", "nodule_and_root_2", "nodule_and_root_3",
    "nodule_whole_1", "nodule_whole_2", "nodule_whole_3", "nodule_whole_4",
    "nodule_nonscen_1", "nodule_nonscen_2",
    "nodule_scen_1", "nodule_scen_2", "nodule_scen_3", "nodule_scen_4",
    "root_tip_1", "root_tip_2", "root_tip_3", "root_tip_4",
    "root_nobranch_1", "root_nobranch_2",
    "root_branched_3", "root_branched_4", "root_branched_5", "root_branched_6",
    "shoot_tip_4d_1", "shoot_tip_4d_2", "shoot_tip_4d_3",
    "shoot_tip_age2_1", "shoot_tip_age2_2", "shoot_tip_age2_3", "shoot_tip_age2_4", "shoot_tip_age2_5",
    "shoot_tip_age4_1", "shoot_tip_age4_2",
    "shoot_tip_age6_1", "shoot_tip_age6_2", "shoot_tip_age6_3", "shoot_tip_age6_4",
    "shoot_tip_age8_1", "shoot_tip_age8_2", "shoot_tip_age8_3", "shoot_tip_age8_4",
    "shoot_tip_age12_16_1", "shoot_tip_age12_16_2", "shoot_tip_age12_16_3", "shoot_tip_age12_16_4", "shoot_tip_age12_16_5",
    "shoot_tip_age18_21_1", "shoot_tip_age18_21_2", "shoot_tip_age18_21_3",
    "shoot_tip_age28_36_KS_1", "shoot_tip_age28_36_KS_2", "shoot_tip_age28_36_KS_3",
    "shoot_tip_age24_31_OK_1", "shoot_tip_age24_31_OK_2", "shoot_tip_age24_31_OK_3"
  ),
  group = c(
    rep("nodule_and_root", 3),
    rep("nodule_whole", 4),
    rep("nodule_nonscen", 2),
    rep("nodule_scen", 4),
    rep("root_tip", 4),
    rep("root_nobranch", 2),
    rep("root_branched", 4),
    rep("shoot_tip_4d", 3),
    rep("shoot_tip_age2", 5),
    rep("shoot_tip_age4", 2),
    rep("shoot_tip_age6", 4),
    rep("shoot_tip_age8", 4),
    rep("shoot_tip_age12_16", 5),
    rep("shoot_tip_age18_21", 3),
    rep("shoot_tip_age28_36_KS", 3),
    rep("shoot_tip_age24_31_OK", 3)
  ),
  stringsAsFactors = FALSE
)

# Locate quant.sf files
files <- file.path(base_dir, paste0(sample_info$sample, "_quant"), "quant.sf")
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

# Remove samples with zero library size
lib_sizes <- colSums(txi$counts)
valid_samples <- lib_sizes > 0

if (any(!valid_samples)) {
  cat("Removing samples with zero library size:
", paste(names(lib_sizes[!valid_samples]), collapse = ", "), "
")
  flush.console()
  txi$counts <- txi$counts[, valid_samples, drop = FALSE]
  txi$length <- txi$length[, valid_samples, drop = FALSE]
  sample_info <- sample_info[valid_samples, , drop = FALSE]
}

# Inspect the filtered dataset before running edgeR
print(dim(txi$counts))
print(table(sample_info$group))

if (ncol(txi$counts) == 0) {
  stop("All samples were removed. Check the dataset.")
}

# Run edgeR
group <- factor(sample_info$group)
dge <- DGEList(counts = txi$counts, group = group)
dge <- calcNormFactors(dge, method = "TMM")
design <- model.matrix(~group)
dge <- estimateDisp(dge, design, robust = TRUE)
fit <- glmFit(dge, design)
lrt <- glmLRT(fit)

# Save differential expression results
write_csv(topTags(lrt, n = Inf)$table, "1.DEG_results_edgeR.csv")
cat("
Analysis completed. Results were saved to '1.DEG_results_edgeR.csv'.
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

####################################################
## Generate the Chamaecrista expression table
python3 Apply_expression_data_for_chafa_RNA-seq.py
