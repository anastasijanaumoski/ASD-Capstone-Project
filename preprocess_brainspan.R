# Libraries
library(dplyr)
library(tidyr)
library(readr)
library(tibble)

# Load datasets
brainspan_matrix <- read_csv(
  "data/BrainSpan_RNAseq_v10/expression_matrix.csv",
  col_names = FALSE,
  show_col_types = FALSE
)

brainspan_columns <- read_csv(
  "data/BrainSpan_RNAseq_v10/columns_metadata.csv",
  show_col_types = FALSE
)

brainspan_rows <- read_csv(
  "data/BrainSpan_RNAseq_v10/rows_metadata.csv",
  show_col_types = FALSE
)

# Fix expression matrix
# First column in expression_matrix.csv is gene_id
colnames(brainspan_matrix)[1] <- "gene_id"

# All other columns are samples
sample_ids <- as.character(brainspan_columns$column_num)
colnames(brainspan_matrix)[-1] <- sample_ids

# Attach gene metadata 
brainspan_rows_small <- brainspan_rows %>%
  select(gene_id, ensembl_gene_id, gene_symbol)

expr_annot <- brainspan_rows_small %>%
  right_join(brainspan_matrix, by = "gene_id")

# Keep only valid Ensembl IDs
expr_annot <- expr_annot %>% filter(!is.na(ensembl_gene_id))

# Make expression numeric
expr_cols <- sample_ids
expr_annot[expr_cols] <- lapply(expr_annot[expr_cols], as.numeric)

# Log transform
expr_annot[expr_cols] <- log1p(as.matrix(expr_annot[expr_cols]))

# Remove genes with zero expression
gene_sums <- rowSums(expr_annot[expr_cols], na.rm = TRUE)
expr_annot <- expr_annot[gene_sums > 0, ]

# Clean metadata
brainspan_columns_processed <- brainspan_columns %>%
  mutate(
    sample_id = as.character(column_num),
    region = structure_name,
    donor_id = as.character(donor_id)
  )

# Long format
long_df <- expr_annot %>%
  select(ensembl_gene_id, all_of(expr_cols)) %>%
  pivot_longer(
    cols = all_of(expr_cols),
    names_to = "sample_id",
    values_to = "expression"
  ) %>%
  left_join(
    brainspan_columns_processed %>%
      select(sample_id, region, donor_id),
    by = "sample_id"
  )

# Filter regions with at least 2 donors
region_donor_counts <- long_df %>%
  group_by(region) %>%
  summarise(n_donors = n_distinct(donor_id), .groups = "drop")

valid_regions <- region_donor_counts %>%
  filter(n_donors >= 2) %>%
  pull(region)

long_df <- long_df %>%
  filter(region %in% valid_regions)

# SD per gene per region
gene_sd <- long_df %>%
  group_by(region, ensembl_gene_id) %>%
  summarise(
    sd_expression = sd(expression, na.rm = TRUE),
    .groups = "drop"
  )

# Z-score variability
gene_sd_z <- gene_sd %>%
  group_by(region) %>%
  mutate(
    sd_mean = mean(sd_expression, na.rm = TRUE),
    sd_sd   = sd(sd_expression, na.rm = TRUE),
    sd_z    = ifelse(sd_sd == 0, NA_real_,
                     (sd_expression - sd_mean) / sd_sd)
  ) %>%
  ungroup()

# Z cutoff
z_cutoff <- 2

hv_genes <- gene_sd_z %>%
  filter(!is.na(sd_z), sd_z >= z_cutoff)

# Output directory
out_dir <- "data/processed"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Save outputs (IMPORTANT for downstream scripts)
write_csv(hv_genes,
          file.path(out_dir, "highlyvariable_genes_zscore_sd_ge2.csv"))

saveRDS(expr_annot,
        file.path(out_dir, "expr_annot.rds"))

saveRDS(brainspan_columns_processed,
        file.path(out_dir, "brainspan_columns_processed.rds"))