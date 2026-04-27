library(igraph)
library(dplyr)
library(purrr)
library(tibble)
library(readr)

# setting
network_dir <- "data/network_threshold/thr_0.85"
out_dir <- "results/random_network_validation"
n_random <- 1000

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# functions
safe_clustering <- function(g) {
  if (ecount(g) == 0) return(NA_real_)
  transitivity(g, type = "average")
}

safe_mean_path <- function(g) {
  if (ecount(g) == 0 || vcount(g) < 2) return(NA_real_)
  
  comps <- components(g)
  giant_nodes <- V(g)[comps$membership == which.max(comps$csize)]
  g_giant <- induced_subgraph(g, giant_nodes)
  
  if (vcount(g_giant) < 2 || ecount(g_giant) == 0) return(NA_real_)
  
  suppressWarnings(mean_distance(g_giant, directed = FALSE))
}

safe_assortativity <- function(g) {
  if (ecount(g) == 0) return(NA_real_)
  suppressWarnings(assortativity_degree(g, directed = FALSE))
}

compute_metrics <- function(g) {
  tibble(
    n_nodes = vcount(g),
    n_edges = ecount(g),
    clustering = safe_clustering(g),
    mean_path = safe_mean_path(g),
    assortativity = safe_assortativity(g)
  )
}

randomize_graph <- function(g) {
  rewire(g, with = keeping_degseq(niter = ecount(g) * 10))
}

z_score <- function(obs, rand_vals) {
  mu <- mean(rand_vals, na.rm = TRUE)
  sdv <- sd(rand_vals, na.rm = TRUE)
  if (is.na(obs) || is.na(mu) || is.na(sdv) || sdv == 0) return(NA_real_)
  (obs - mu) / sdv
}

empirical_p_upper <- function(obs, rand_vals) {
  if (is.na(obs)) return(NA_real_)
  (sum(rand_vals >= obs, na.rm = TRUE) + 1) / (sum(!is.na(rand_vals)) + 1)
}

empirical_p_lower <- function(obs, rand_vals) {
  if (is.na(obs)) return(NA_real_)
  (sum(rand_vals <= obs, na.rm = TRUE) + 1) / (sum(!is.na(rand_vals)) + 1)
}

# main
files <- list.files(
  network_dir,
  pattern = "^graph_.*thr0.85\\.rds$",
  full.names = TRUE
)

if (length(files) == 0) {
  stop("No graph RDS files found. Check that graph objects were saved in data/network_threshold/thr_0.85.")
}

results <- map_dfr(files, function(file) {
  region <- basename(file)
  region <- sub("^graph_", "", region)
  region <- sub("_thr0.85\\.rds$", "", region)
  
  message("Processing: ", region)
  
  g <- readRDS(file)
  
  if (vcount(g) < 10 || ecount(g) < 10) {
    return(tibble(
      region = region,
      n_nodes = vcount(g),
      n_edges = ecount(g),
      obs_clustering = NA_real_,
      rand_clustering_mean = NA_real_,
      clustering_z = NA_real_,
      clustering_p = NA_real_,
      obs_mean_path = NA_real_,
      rand_mean_path_mean = NA_real_,
      mean_path_z = NA_real_,
      mean_path_p_upper = NA_real_,
      mean_path_p_lower = NA_real_,
      obs_assortativity = NA_real_,
      rand_assortativity_mean = NA_real_,
      assortativity_z = NA_real_,
      assortativity_p_upper = NA_real_,
      assortativity_p_lower = NA_real_,
      note = "too small"
    ))
  }
  
  obs <- compute_metrics(g)
  
  rand <- replicate(n_random, {
    g_rand <- randomize_graph(g)
    compute_metrics(g_rand)
  }, simplify = FALSE) %>%
    bind_rows()
  
  tibble(
    region = region,
    n_nodes = obs$n_nodes,
    n_edges = obs$n_edges,
    
    obs_clustering = obs$clustering,
    rand_clustering_mean = mean(rand$clustering, na.rm = TRUE),
    clustering_z = z_score(obs$clustering, rand$clustering),
    clustering_p = empirical_p_upper(obs$clustering, rand$clustering),
    
    obs_mean_path = obs$mean_path,
    rand_mean_path_mean = mean(rand$mean_path, na.rm = TRUE),
    mean_path_z = z_score(obs$mean_path, rand$mean_path),
    mean_path_p_upper = empirical_p_upper(obs$mean_path, rand$mean_path),
    mean_path_p_lower = empirical_p_lower(obs$mean_path, rand$mean_path),
    
    obs_assortativity = obs$assortativity,
    rand_assortativity_mean = mean(rand$assortativity, na.rm = TRUE),
    assortativity_z = z_score(obs$assortativity, rand$assortativity),
    assortativity_p_upper = empirical_p_upper(obs$assortativity, rand$assortativity),
    assortativity_p_lower = empirical_p_lower(obs$assortativity, rand$assortativity),
    
    note = NA_character_
  )
})

print(results)

write_csv(results,
          file.path(out_dir, "random_network_validation_1000.csv"))