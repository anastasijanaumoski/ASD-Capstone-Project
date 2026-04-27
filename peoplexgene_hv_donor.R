## PEOPLE × GENE (Highly Variable Genes Per Region; DONOR × GENE)
library(dplyr)
library(readr)

# Load processed datasets
expr_annot <- readRDS("data/processed/expr_annot.rds")
brainspan_columns_processed <- readRDS("data/processed/brainspan_columns_processed.rds")

# Load HV genes
hv_genes <- read_csv(
  "data/processed/highlyvariable_genes_zscore_sd_ge2.csv",
  show_col_types = FALSE
)

# Output directory
outdir <- "data/peoplexgene_hv_donor"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

regions <- unique(hv_genes$region)

for (reg in regions) {
  cat("Processing region:", reg, "\n")
  
  clean_name <- gsub("[^A-Za-z0-9_]", "_", reg)
  
  # Metadata for region
  region_meta <- brainspan_columns_processed %>%
    filter(region == reg) %>%
    select(sample_id, donor_id)
  
  # Skip small regions
  if (n_distinct(region_meta$donor_id) < 2) {
    cat("Skipping:", reg, "\n")
    next
  }
  
  region_sample_ids <- region_meta$sample_id
  
  # HV genes
  region_genes <- hv_genes %>%
    filter(region == reg) %>%
    pull(ensembl_gene_id) %>%
    as.character()
  
  # Subset expression
  submat <- expr_annot %>%
    filter(ensembl_gene_id %in% region_genes) %>%
    select(ensembl_gene_id, all_of(region_sample_ids))
  
  # Matrix
  mat <- as.matrix(submat[, -1])
  rownames(mat) <- submat$ensembl_gene_id
  
  # Samples × genes
  sxg <- t(mat)
  
  sxg_df <- as.data.frame(sxg)
  sxg_df$sample_id <- rownames(sxg_df)
  
  # Ensure numeric
  sxg_df <- sxg_df %>%
    mutate(across(-sample_id, as.numeric))
  
  # Aggregate to donor × gene
  donor_x_gene <- sxg_df %>%
    left_join(region_meta, by = "sample_id") %>%
    select(-sample_id) %>%
    group_by(donor_id) %>%
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE)),
              .groups = "drop")
  
  pxg <- as.matrix(donor_x_gene[, -1])
  rownames(pxg) <- donor_x_gene$donor_id
  
  # Save
  saveRDS(pxg,
          file.path(outdir,
                    paste0("peoplexgene_", clean_name, "_hv_donor.rds")))
}