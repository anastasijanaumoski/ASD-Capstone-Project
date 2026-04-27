# Libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(poweRlaw)
})

# Load node metrics from final threshold
df <- read_csv(
  "data/network_threshold/thr_0.85/node_metrics_ALL_thr0.85.csv",
  show_col_types = FALSE
) %>%
  mutate(
    region = as.character(region),
    gene = as.character(gene),
    degree = as.numeric(degree),
    betweenness = as.numeric(betweenness)
  )

# Output folders
plot_dir <- "figures/scaling"
fit_plot_dir <- file.path(plot_dir, "powerlaw_fits")
result_dir <- "results/scaling"

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fit_plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

# Use raw degree
degree_dist <- df %>%
  group_by(region, degree) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(region) %>%
  mutate(pk = count / sum(count)) %>%
  filter(degree > 0, pk > 0)

p_loglog <- ggplot(degree_dist, aes(x = degree, y = pk)) +
  geom_point(alpha = 0.6) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~ region, scales = "free") +
  labs(
    title = "Log–Log Degree Distribution by Region",
    x = "Degree (log scale)",
    y = "P(k) (log scale)"
  )

print(p_loglog)

ggsave(
  filename = file.path(plot_dir, "loglog_degree_distribution_by_region.png"),
  plot = p_loglog,
  width = 16,
  height = 12,
  dpi = 300
)

# Fit power law and compare to lognormal/exponential
fit_powerlaw_region <- function(deg, region_name, n_boot = 100) {
  
  deg <- deg[is.finite(deg)]
  deg <- as.integer(round(deg))
  deg <- deg[deg > 0]
  
  if (length(deg) < 50 || length(unique(deg)) < 10) {
    return(tibble(
      n_nodes = length(deg),
      xmin = NA_real_,
      alpha = NA_real_,
      p_boot = NA_real_,
      R_pl_vs_lognorm = NA_real_,
      p_pl_vs_lognorm = NA_real_,
      R_pl_vs_exp = NA_real_,
      p_pl_vs_exp = NA_real_,
      note = "Too few nodes / unique degrees for reliable fit"
    ))
  }
  
  out <- tryCatch({
    
    m_pl <- displ$new(deg)
    
    est_xmin_pl <- estimate_xmin(m_pl)
    m_pl$setXmin(est_xmin_pl$xmin)
    
    est_pars_pl <- estimate_pars(m_pl)
    m_pl$setPars(est_pars_pl$pars)
    
    bs <- bootstrap_p(m_pl, no_of_sims = n_boot, threads = 1)
    p_val <- bs$p
    
    xmin <- m_pl$getXmin()
    
    m_ln <- dislnorm$new(deg)
    m_ln$setXmin(xmin)
    m_ln$setPars(estimate_pars(m_ln)$pars)
    
    m_exp <- disexp$new(deg)
    m_exp$setXmin(xmin)
    m_exp$setPars(estimate_pars(m_exp)$pars)
    
    comp_pl_ln <- compare_distributions(m_pl, m_ln)
    comp_pl_ex <- compare_distributions(m_pl, m_exp)
    
    k_vals <- sort(unique(deg))
    df_ccdf <- data.frame(
      k = k_vals,
      ccdf = sapply(k_vals, function(x) mean(deg >= x))
    )
    
    k_tail <- df_ccdf$k[df_ccdf$k >= xmin]
    ccdf_fit <- 1 - dist_cdf(m_pl, q = k_tail - 1)
    df_fit <- data.frame(k = k_tail, ccdf_fit = ccdf_fit)
    
    alpha <- m_pl$getPars()
    
    p <- ggplot(df_ccdf, aes(x = k, y = ccdf)) +
      geom_point(size = 1.2, alpha = 0.8) +
      geom_line(data = df_fit, aes(x = k, y = ccdf_fit), linewidth = 1) +
      scale_x_log10() +
      scale_y_log10() +
      labs(
        title = paste0("Degree CCDF + power-law fit — ", region_name),
        subtitle = paste0(
          "xmin=", xmin,
          ", alpha=", round(alpha, 3),
          ", bootstrap p=", signif(p_val, 3)
        ),
        x = "Degree (log scale)",
        y = "P(K ≥ k) (log scale)"
      )
    
    ggsave(
      filename = file.path(fit_plot_dir, paste0("powerlaw_ccdf_", region_name, ".png")),
      plot = p,
      width = 7.5,
      height = 5.5,
      dpi = 300
    )
    
    tibble(
      n_nodes = length(deg),
      xmin = xmin,
      alpha = alpha,
      p_boot = p_val,
      R_pl_vs_lognorm = comp_pl_ln$test_statistic,
      p_pl_vs_lognorm = comp_pl_ln$p_two_sided,
      R_pl_vs_exp = comp_pl_ex$test_statistic,
      p_pl_vs_exp = comp_pl_ex$p_two_sided,
      note = NA_character_
    )
    
  }, error = function(e) {
    tibble(
      n_nodes = length(deg),
      xmin = NA_real_,
      alpha = NA_real_,
      p_boot = NA_real_,
      R_pl_vs_lognorm = NA_real_,
      p_pl_vs_lognorm = NA_real_,
      R_pl_vs_exp = NA_real_,
      p_pl_vs_exp = NA_real_,
      note = paste("Fit failed:", conditionMessage(e))
    )
  })
  
  out
}

# Run per region
results <- df %>%
  group_by(region) %>%
  summarise(deg = list(degree), .groups = "drop") %>%
  mutate(fit = map2(deg, region, ~ fit_powerlaw_region(.x, .y, n_boot = 100))) %>%
  select(-deg) %>%
  unnest(fit)

# Save summary table
write_csv(
  results,
  file.path(result_dir, "powerlaw_fit_summary_by_region.csv")
)

summary(results)