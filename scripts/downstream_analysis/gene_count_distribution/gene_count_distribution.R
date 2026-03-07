# Purpose: Generate gene-count distribution plots across species and perform
# nonparametric statistical tests.
# Input: Excel file containing per-species gene counts.
# Output: Boxplots (PNG/PDF), Kruskal-Wallis summary, and pairwise Wilcoxon p-values.

library(readxl)
library(tidyr)
library(ggplot2)
library(dplyr)

# Load input data
file_path <- "count_only_07-10-25.xlsx"  # Update with the actual file path if needed.
df <- read_excel(file_path, sheet = "Sheet1")

# Reshape the gene-count columns (columns 5 onward) and standardize species labels
# from abbreviated prefixes to genus names.
df_long <- df %>%
  pivot_longer(cols = 5:ncol(df), names_to = "Species", values_to = "Gene_Count") %>%
  mutate(Species = sub("\..*$", "", Species))

# Map abbreviated prefixes to genus names
species_names <- c(
  "acacr" = "Acacia",
  "aesev" = "Aeschynomene",
  "arath" = "Arabidopsis",
  "bauva" = "Bauhinia",
  "cerca" = "Cercis",
  "chafa" = "Chamaecrista",
  "glyma" = "Glycine",
  "lotja" = "Lotus",
  "medtr" = "Medicago",
  "paran" = "Parasponia",
  "phach" = "Phanera",
  "phavu" = "Phaseolus",
  "prupe" = "Prunus",
  "quisa" = "Quillaja",
  "sento" = "Senna",
  "singl" = "Sindora",
  "treor" = "Trema",
  "vitvi" = "Vitis"
)
df_long$Species <- species_names[df_long$Species]

# Apply the requested species order
desired_order <- c("Lotus", "Medicago", "Glycine", "Phaseolus", "Aeschynomene", "Chamaecrista",
                   "Senna", "Acacia", "Bauhinia", "Phanera", "Cercis", "Sindora",
                   "Quillaja", "Parasponia", "Trema", "Prunus", "Arabidopsis", "Vitis")

df_long$Species <- factor(df_long$Species, levels = desired_order)

# Statistical testing: Kruskal-Wallis followed by pairwise Wilcoxon tests (BH correction)
df_kw <- df_long %>%
  filter(!is.na(Species), !is.na(Gene_Count))

# Sample counts per species
species_n <- df_kw %>% count(Species) %>% arrange(match(Species, desired_order))
print(species_n)

# Kruskal-Wallis test across species
kw_res <- kruskal.test(Gene_Count ~ Species, data = df_kw)
print(kw_res)

# Effect size (epsilon^2) for the Kruskal-Wallis test
k <- nlevels(droplevels(df_kw$Species))
N <- nrow(df_kw)
H <- as.numeric(kw_res$statistic)
epsilon2 <- (H - (k - 1)) / (N - 1)
cat(sprintf("Kruskal-Wallis epsilon^2 (effect size): %.4f\n", epsilon2))

# Pairwise Wilcoxon tests with BH correction
pw_res <- pairwise.wilcox.test(df_kw$Gene_Count, df_kw$Species, p.adjust.method = "BH")
print(pw_res)

# Save statistical outputs
sink("KW_test_summary.txt")
cat("=== Kruskal-Wallis Test Summary ===\n")
print(kw_res)
cat(sprintf("\nEpsilon^2 (effect size): %.4f\n", epsilon2))
cat("\n=== Sample sizes per species ===\n")
print(species_n)
sink()

# Save the pairwise Wilcoxon p-value matrix
pw_mat <- pw_res$p.value
write.csv(pw_mat, file = "Pairwise_Wilcoxon_BH_pvalues.csv", row.names = TRUE)

# Visualization: boxplot on a log scale
plot_log <- ggplot(df_long, aes(x = Species, y = Gene_Count + 1, fill = Species)) +
  geom_boxplot(outlier.size = 1, outlier.alpha = 0.7) +
  scale_y_log10() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        legend.position = "none") +
  labs(x = "Species", y = "Gene Counts (log scale)",
       title = "Distribution of Gene Counts per Species (log scale)")

ggsave("gene_counts_distribution_log_08-05-25.png", plot = plot_log, width = 10, height = 7, dpi = 300)
ggsave("gene_counts_distribution_log_08-05-25.pdf", plot = plot_log, width = 10, height = 7)

# Visualization: boxplot with counts capped at 50
cap_value <- 50
df_long_capped <- df_long %>%
  mutate(Gene_Count = ifelse(Gene_Count > cap_value, cap_value, Gene_Count))

plot_original_capped <- ggplot(df_long_capped, aes(x = Species, y = Gene_Count, fill = Species)) +
  geom_boxplot(outlier.size = 1, outlier.alpha = 0.7) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        legend.position = "none") +
  labs(x = "Species", y = "Gene Counts (Capped at 50)",
       title = "Distribution of Gene Counts per Species (Capped at 50)")

ggsave("gene_counts_distribution_original_capped_08-05-25.png", plot = plot_original_capped, width = 10, height = 7, dpi = 300)
ggsave("gene_counts_distribution_original_capped_08-05-25.pdf", plot = plot_original_capped, width = 10, height = 7)
