library(dplyr)
library(readr)
library(tidyr)

# Load data
sfari <- read_csv(
  "data/SFARI/SFARI-Gene_genes_01-14-2026release_02-17-2026export.csv",
  show_col_types = FALSE
)

brainspan_rows <- read_csv(
  "data/BrainSpan_RNAseq_v10/rows_metadata.csv",
  show_col_types = FALSE
) %>%
  mutate(
    gene_symbol = as.character(gene_symbol),
    ensembl_gene_id = as.character(ensembl_gene_id)
  )

gene_across_regions_z <- read_csv(
  "results/sensitivity_analysis/gene_across_regions_thr085_zge2_sensitivity_topPcts.csv",
  show_col_types = FALSE
) %>%
  filter(top_pct == 0.05) %>%
  rename(
    n_regions_high_degree_z = regions_high_degree,
    n_regions_high_bet_z = regions_high_bet,
    n_regions_high_both_z = regions_high_both
  )

# Map SFARI gene symbols to BrainSpan Ensembl IDs
sfari_map <- sfari %>%
  transmute(
    sfari_symbol = as.character(`gene-symbol`),
    gene_score = suppressWarnings(as.numeric(`gene-score`)),
    status = as.character(status)
  ) %>%
  distinct(sfari_symbol, .keep_all = TRUE) %>%
  left_join(
    brainspan_rows %>% select(gene_symbol, ensembl_gene_id),
    by = c("sfari_symbol" = "gene_symbol")
  )

# High-confidence SFARI genes
sfari_ids_hc <- sfari_map %>%
  filter(!is.na(ensembl_gene_id), gene_score %in% c(1, 2)) %>%
  pull(ensembl_gene_id) %>%
  unique()

# All mapped SFARI genes
sfari_ids_all <- sfari_map %>%
  filter(!is.na(ensembl_gene_id)) %>%
  pull(ensembl_gene_id) %>%
  unique()

# Annotate z-score hub recurrence table
gene_across_regions_z2 <- gene_across_regions_z %>%
  mutate(
    is_sfari_hc = gene %in% sfari_ids_hc,
    is_sfari_any = gene %in% sfari_ids_all
  )

# Output directory
out_dir <- "results/sfari_overlap"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

write_csv(
  gene_across_regions_z2,
  file.path(out_dir, "gene_across_regions_thr085_zscore_ge2_sfari.csv")
)

# Enrichment analysis
run_sfari_enrichment <- function(df, min_regions, sfari_col = "is_sfari_any") {
  
  conserved <- df %>%
    mutate(conserved = n_regions_high_degree_z >= min_regions)
  
  a <- sum(conserved$conserved & conserved[[sfari_col]], na.rm = TRUE)
  b <- sum(conserved$conserved & !conserved[[sfari_col]], na.rm = TRUE)
  c <- sum(!conserved$conserved & conserved[[sfari_col]], na.rm = TRUE)
  d <- sum(!conserved$conserved & !conserved[[sfari_col]], na.rm = TRUE)
  
  mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  rownames(mat) <- c("Conserved", "Not_conserved")
  colnames(mat) <- c("SFARI", "Not_SFARI")
  
  ft <- fisher.test(mat)
  
  tibble(
    min_regions = min_regions,
    sfari_set = sfari_col,
    conserved_genes = a + b,
    sfari_in_conserved = a,
    odds_ratio = unname(ft$estimate),
    p_value = ft$p.value,
    conf_low = ft$conf.int[1],
    conf_high = ft$conf.int[2]
  )
}

# Enrichment for cutoffs
cutoffs <- c(2, 3, 4, 5)

enrichment_results <- bind_rows(
  lapply(cutoffs, function(k) run_sfari_enrichment(gene_across_regions_z2, k, "is_sfari_any")),
  lapply(cutoffs, function(k) run_sfari_enrichment(gene_across_regions_z2, k, "is_sfari_hc"))
)

print(enrichment_results)

write_csv(
  enrichment_results,
  file.path(out_dir, "sfari_enrichment_conserved_degree_hubs_zge2.csv")
)

# Extract conserved hub lists for reporting
conserved_hubs_ge3 <- gene_across_regions_z2 %>%
  filter(n_regions_high_degree_z >= 3) %>%
  arrange(desc(n_regions_high_degree_z), desc(is_sfari_hc), desc(is_sfari_any))

conserved_hubs_ge4 <- gene_across_regions_z2 %>%
  filter(n_regions_high_degree_z >= 4) %>%
  arrange(desc(n_regions_high_degree_z), desc(is_sfari_hc), desc(is_sfari_any))

print(conserved_hubs_ge3)
print(conserved_hubs_ge4)

write_csv(
  conserved_hubs_ge3,
  file.path(out_dir, "conserved_degree_hubs_zge2_ge3regions.csv")
)

write_csv(
  conserved_hubs_ge4,
  file.path(out_dir, "conserved_degree_hubs_zge2_ge4regions.csv")
)