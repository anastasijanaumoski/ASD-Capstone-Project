library(readr)

people_dir <- "data/peoplexgene_hv_donor"
genexgene_dir <- "data/genexgene_hv"

dir.create(genexgene_dir, showWarnings = FALSE, recursive = TRUE)

files <- list.files(people_dir,
                    pattern = "peoplexgene_.*_hv_donor\\.rds$",
                    full.names = TRUE)

for (f in files) {
  
  region_clean <- gsub("peoplexgene_|_hv_donor\\.rds", "", basename(f))
  
  cat("Processing:", region_clean, "\n")
  
  pxg <- readRDS(f)
  
  # remove genes with zero variance (prevents NA correlations)
  pxg <- pxg[, apply(pxg, 2, sd, na.rm = TRUE) > 0]
  
  # Pearson correlation across genes
  gxg <- cor(pxg, use = "pairwise.complete.obs", method = "pearson")
  
  saveRDS(gxg,
          file.path(genexgene_dir,
                    paste0("genexgene_", region_clean, "_hv.rds")))
}