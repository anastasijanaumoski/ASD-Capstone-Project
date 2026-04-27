# Libraries
library(dplyr)
library(tidyr)
library(readr)
library(tibble)

# Load dataset(s)
brainspan_matrix <- read_csv("/N/slate/annaum/data/transcriptomics/BrainSpan_RNAseq_v10/expression_matrix.csv",
                             col_names = FALSE,
                             show_col_types = FALSE)

brainspan_columns <- read_csv("/N/slate/annaum/data/transcriptomics/BrainSpan_RNAseq_v10/columns_metadata.csv",
                              show_col_types = FALSE)

brainspan_rows <- read_csv("/N/slate/annaum/data/transcriptomics/BrainSpan_RNAseq_v10/rows_metadata.csv",
                           show_col_types = FALSE)

# Figure out number of subjects
unique_subjects <- unique(brainspan_columns$donor_id)
num_unique_subjects <- length(unique_subjects)

print(num_unique_subjects) # gives us 42 subjects in total

## Separate the matrix to be person, time, and region specific

# need to remove a column from the matrix dataset (525 columns total)
# because brainspan_columns has 524 samples
brainspan_matrix <- brainspan_matrix[, -1] # removes first column (index column)

# convert matrix to dataframe, tibble causes a warning message
brainspan_matrix <- as.data.frame(brainspan_matrix)

# assign names to matrix & columns
colnames(brainspan_matrix) <- brainspan_columns$column_num
rownames(brainspan_matrix) <- brainspan_rows$ensembl_gene_id

# check distribution of raw expression values
hist(as.numeric(unlist(brainspan_matrix)),
     breaks = 100,
     main = "Distribution of Raw Expression Values",
     xlab = "Expression Level")

# log transformation (only done once)
brainspan_matrix <- log1p(brainspan_matrix)

# check distribution after log transformation
hist(as.numeric(unlist(brainspan_matrix)),
     breaks = 100,
     main = "Distribution of Log-Transformed Expression Values",
     xlab = "Expression Level")

# clean metadata + create useful columns
brainspan_columns <- brainspan_columns %>%
  mutate(sample_id = as.character(column_num),
         person = donor_id,
         time = age,
         region = structure_name,
         period = case_when(
           grepl("pcw", time) ~ "Prenatal",
           grepl("mos|yrs", time) & as.numeric(gsub("[^0-9]", "", time)) < 24 ~ "Infant",
           TRUE ~ "Child/Adult"
         ))

## For each region and time period, calculate the SD (variance) of each gene across the 42 people

# 1. Convert matrix to long format and join metadata
long_df <- brainspan_matrix %>%
  rownames_to_column("ensembl_gene_id") %>%
  pivot_longer(
    cols = -ensembl_gene_id,
    names_to = "sample_id",
    values_to = "expression"
  ) %>%
  left_join(brainspan_columns %>%
              select(sample_id, person, region, time, period),
            by = "sample_id")

# 2. Compute SD across people for each gene and region
group_by_sd_region <- long_df %>%
  group_by(ensembl_gene_id, region) %>%
  summarize(
    mean_expression = mean(expression, na.rm = TRUE),
    sd_expression   = sd(expression, na.rm = TRUE),
    var_expression  = var(expression, na.rm = TRUE),
    n               = n(), # number of samples per region
    n_donors        = n_distinct(person)
  ) %>%
  ungroup()

# 3. Extract the top 100 genes based on SD
top_genes_sd_region <- group_by_sd_region %>%
  group_by(region) %>%
  arrange(desc(sd_expression)) %>%
  slice_head(n = 100) %>%
  ungroup()

write.csv(top_genes_sd_region,
          "/N/slate/annaum/data/transcriptomics/logtrans_top_genes_sd_region.csv",
          row.names = FALSE)

### sum expression across all samples for each gene
gene_sums <- rowSums(brainspan_matrix, na.rm = TRUE)
summary(gene_sums)

hist(gene_sums,
     breaks = 100,
     main = "Distribution of Gene Expression Sums Across All Samples",
     xlab = "Summed Expression per Gene")

sum(gene_sums == 0)
mean(gene_sums == 0)