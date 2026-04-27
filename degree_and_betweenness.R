## Thresholded BrainSpan co-expression networks (HV genes, z >= 2)
## Builds graphs at multiple correlation thresholds and computes degree + betweenness
## Outputs per-threshold combined tables:
##   - node_metrics_ALL_thrX.csv
##   - graph_summaries_ALL_thrX.csv

# Libraries
library(dplyr)
library(readr)
library(tibble)
library(igraph)

BASE <- "data"

genexgene_dir <- file.path(BASE, "genexgene_hv")
out_base <- file.path(BASE, "network_threshold")

# Output: thresholded networks + centrality metrics
out_base <- file.path(BASE, "network_threshold")
dir.create(out_base, showWarnings = FALSE, recursive = TRUE)

thresholds <- c(0.7, 0.85, 0.9)

# Helper: region name from filename
# Assumes files named like: genexgene_<REGION>_hv.rds  OR genexgene_<REGION>.rds
get_region <- function(path) {
  bn <- basename(path)
  bn <- sub("^genexgene_", "", bn)
  bn <- sub("_hv\\.rds$", "", bn)
  bn <- sub("\\.rds$", "", bn)
  bn
}

files <- list.files(genexgene_dir, pattern = "\\.rds$", full.names = TRUE)
cat("Input dir:", genexgene_dir, "\n")
cat("Found", length(files), "RDS gxg files\n")
if (length(files) == 0) stop("No .rds files found in genexgene_dir. Check path and filenames.")

cat("Example files:\n")
print(head(files))

for (thr in thresholds) {
  
  out_dir <- file.path(out_base, paste0("thr_", thr))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  all_node_metrics <- list()
  all_graph_summaries <- list()
  
  cat("\n===== Running threshold:", thr, "=====\n")
  
  for (f in files) {
    
    region_clean <- get_region(f)
    cat("Processing region:", region_clean, "\n")
    
    gxg <- readRDS(f)
    
# adjacency: keep strong absolute correlations
    adj <- abs(gxg) >= thr
    diag(adj) <- FALSE
    adj[is.na(adj)] <- FALSE
    
    # Ensure plain matrix (prevents class-related issues)
    adj <- as.matrix(adj)
    
# build thresholded graph (undirected, unweighted)
    g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
    
# node metrics
    deg <- igraph::degree(g)
    btw <- igraph::betweenness(g, directed = FALSE, normalized = TRUE)
    
    node_metrics <- tibble::tibble(
      region = region_clean,
      gene = names(deg),
      degree = as.numeric(deg),
      betweenness = as.numeric(btw)
    )
    
# force igraph
    comps <- igraph::components(g)
    
    summ <- tibble::tibble(
      region = region_clean,
      threshold = thr,
      n_genes = igraph::vcount(g),
      n_edges = igraph::ecount(g),
      density = igraph::edge_density(g, loops = FALSE),
      n_components = comps$no,
      largest_component_size = max(comps$csize)
    )
    
    all_node_metrics[[region_clean]] <- node_metrics
    all_graph_summaries[[region_clean]] <- summ
  }
  
  # --- combine and save per-threshold outputs ---
  node_metrics_all <- dplyr::bind_rows(all_node_metrics)
  graph_summaries_all <- dplyr::bind_rows(all_graph_summaries)
  
  write_csv(node_metrics_all, file.path(out_dir, paste0("node_metrics_ALL_thr", thr, ".csv")))
  write_csv(graph_summaries_all, file.path(out_dir, paste0("graph_summaries_ALL_thr", thr, ".csv")))
  
  cat("Saved outputs to:", out_dir, "\n")
  cat("  -", paste0("node_metrics_ALL_thr", thr, ".csv"), "\n")
  cat("  -", paste0("graph_summaries_ALL_thr", thr, ".csv"), "\n")
}