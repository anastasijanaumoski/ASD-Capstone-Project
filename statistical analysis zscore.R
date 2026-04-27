# Libraries
library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(ggplot2)

# Load datasets
brainspan_matrix <- read_csv(
  "/N/u/annaum/Quartz/BrainSpan_data_preprocess/BrainSpan_RNAseq_v10/expression_matrix.csv",
  show_col_types = FALSE
)
brainspan_columns <- read_csv(
  "/N/u/annaum/Quartz/BrainSpan_data_preprocess/BrainSpan_RNAseq_v10/columns_metadata.csv",
  show_col_types = FALSE
)
brainspan_rows <- read_csv(
  "/N/u/annaum/Quartz/BrainSpan_data_preprocess/BrainSpan_RNAseq_v10/rows_metadata.csv",
  show_col_types = FALSE
)

# Fix expression matrix
# First column in expression_matrix.csv is gene_id
colnames(brainspan_matrix)[1] <- "gene_id"

# All other columns are samples; name them using column_num from columns_metadata
sample_ids <- as.character(brainspan_columns$column_num)
stopifnot(length(sample_ids) == ncol(brainspan_matrix) - 1)
colnames(brainspan_matrix)[-1] <- sample_ids

# Attach gene metadata 
brainspan_rows_small <- brainspan_rows %>%
  dplyr::select(gene_id, ensembl_gene_id, gene_symbol)

expr_annot <- brainspan_rows_small %>%
  right_join(brainspan_matrix, by = "gene_id")

# Keep only rows with valid Ensembl IDs
expr_annot <- expr_annot %>% filter(!is.na(ensembl_gene_id))

# Log-transform expression
expr_cols <- sample_ids  # vector of all sample column names

# Make sure these columns are numeric
expr_annot[expr_cols] <- lapply(expr_annot[expr_cols], as.numeric)

# Log1p transform
expr_annot[expr_cols] <- log1p(as.matrix(expr_annot[expr_cols]))

# Drop genes with all-zero expression after log transform (safety)
gene_sums <- rowSums(expr_annot[expr_cols], na.rm = TRUE)
expr_annot <- expr_annot[gene_sums > 0, ]

# Add region info to samples
brainspan_columns <- brainspan_columns %>%
  mutate(
    sample_id = as.character(column_num),
    region = structure_name,
    donor_id = as.character(donor_id)
  )

# Long format for SD calc
long_df <- expr_annot %>%
  dplyr::select(ensembl_gene_id, all_of(expr_cols)) %>%
  pivot_longer(
    cols = all_of(expr_cols),
    names_to = "sample_id",
    values_to = "expression"
  ) %>%
  left_join(
    brainspan_columns %>% dplyr::select(sample_id, region, donor_id),
    by = "sample_id"
  )

# Identify regions with at least 2 donors
region_donor_counts <- long_df %>%
  group_by(region) %>%
  summarise(n_donors = n_distinct(donor_id), .groups = "drop")

valid_regions <- region_donor_counts %>%
  filter(n_donors >= 2) %>%
  pull(region)

# Restrict long_df to valid regions only
long_df <- long_df %>%
  filter(region %in% valid_regions)

# SD per gene per region
gene_sd <- long_df %>%
  group_by(region, ensembl_gene_id) %>%
  summarise(
    sd_expression = sd(expression, na.rm = TRUE),
    .groups = "drop"
  )

# Z-score variability within each region
gene_sd_z <- gene_sd %>%
  group_by(region) %>%
  mutate(
    sd_mean = mean(sd_expression, na.rm = TRUE),
    sd_sd   = sd(sd_expression, na.rm = TRUE),
    sd_z    = ifelse(sd_sd == 0, NA_real_, (sd_expression - sd_mean) / sd_sd)
  ) %>%
  ungroup()

# Z-score cutoff
z_cutoff <- 2

# "highly variable" genes using a statistical cutoff
hv_genes <- gene_sd_z %>%
  filter(!is.na(sd_z), sd_z >= z_cutoff)

# Counts per region (answers "how many genes?")
hv_counts <- hv_genes %>%
  count(region, name = "n_hv_genes") %>%
  arrange(desc(n_hv_genes))

# Save outputs (spring2026/brainspan)
out_dir <- "/N/u/annaum/Quartz/spring2026/brainspan"

write_csv(region_donor_counts,
          file.path(out_dir, "region_donor_counts_z2.csv"))

write_csv(gene_sd_z,
          file.path(out_dir, "gene_sd_with_zscore2.csv"))

write_csv(hv_genes,
          file.path(out_dir, "highlyvariable_genes_zscore_sd_ge2.csv"))

write_csv(hv_counts,
          file.path(out_dir, "highlyvariable_gene_counts_by_region_z_ge2.csv"))
