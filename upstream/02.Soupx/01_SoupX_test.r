suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(SoupX)
  library(ggplot2)
  library(dplyr)
})

# =========================
# USER SETTINGS
# =========================
base_dir <- "/mnt/18T/chibao/pbmc/data"
out_base <- file.path(base_dir, "SoupX_outputs")
dir.create(out_base, showWarnings = FALSE, recursive = TRUE)

# Minimal QC (keep permissive; don't over-filter before SoupX)
qc_min_features <- 200
qc_max_features <- 8000
qc_max_mt <- 25

# Clustering (only to provide clusters + UMAP for SoupX)
n_pcs <- 30
cluster_resolution <- 0.8

# SCTransform options
use_glmGamPoi <- TRUE   # set TRUE if you have glmGamPoi installed
vars_to_regress <- "percent.mt"

set.seed(1234)

# =========================
# UTIL: step wrapper for easy debugging
# =========================
run_step <- function(step_name, expr) {
  message("\n--- STEP: ", step_name, " ---")
  tryCatch(
    expr,
    error = function(e) {
      message("❌ ERROR in step: ", step_name)
      message("   ", conditionMessage(e))
      stop(e)
    }
  )
}

# =========================
# HELPERS
# =========================

read_10x_h5_gex <- function(h5_path) {
  # [FIX] Added use.names/unique.features to ensure consistent gene naming
  x <- Seurat::Read10X_h5(h5_path, use.names = TRUE, unique.features = TRUE)
  if (is.list(x)) {
    if ("Gene Expression" %in% names(x)) return(x[["Gene Expression"]])
    return(x[[1]])
  }
  x
}

find_h5_pair <- function(sample_dir) {
  raw <- list.files(sample_dir, pattern = "raw_feature.*\\.h5$", full.names = TRUE, recursive = TRUE)
  fil <- list.files(sample_dir, pattern = "filtered_feature.*\\.h5$", full.names = TRUE, recursive = TRUE)
  
  # Handle typo case "filtered_featured"
  if (length(fil) == 0) {
    fil <- list.files(sample_dir, pattern = "filtered_featured.*\\.h5$", full.names = TRUE, recursive = TRUE)
  }
  
  if (length(raw) == 0 || length(fil) == 0) {
    return(list(raw = NA_character_, filtered = NA_character_))
  }
  
  # Pick shortest path (often avoids duplicates in subfolders)
  raw <- raw[which.min(nchar(raw))]
  fil <- fil[which.min(nchar(fil))]
  list(raw = raw, filtered = fil)
}

save_qc_plot <- function(seu, out_path, title_text) {
  p <- VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
               ncol = 3, pt.size = 0) +
    ggtitle(title_text)
  ggsave(out_path, plot = p, width = 12, height = 4, dpi = 150)
}

# =========================
# SECTION 1: Load matrices
# =========================
load_cellranger_h5 <- function(sample_dir) {
  h5s <- find_h5_pair(sample_dir)
  if (is.na(h5s$raw) || is.na(h5s$filtered)) {
    stop("Missing raw/filtered .h5 in: ", sample_dir)
  }
  tod <- read_10x_h5_gex(h5s$raw)       # raw droplets
  toc <- read_10x_h5_gex(h5s$filtered)  # filtered cells
  list(tod = tod, toc = toc, raw_h5 = h5s$raw, filtered_h5 = h5s$filtered)
}

# =========================
# SECTION 2: Build Seurat + QC + SCT clustering
# =========================
build_seurat_for_soupx <- function(toc, sample_id, out_dir) {
  seu <- CreateSeuratObject(counts = toc, project = sample_id, min.cells = 3, min.features = 0)
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

  # Save PRE-QC plot
  save_qc_plot(
    seu,
    out_path = file.path(out_dir, paste0(sample_id, "_qc_before_SoupX.png")),
    title_text = paste0(sample_id, " QC (before SoupX)")
  )

  # Filter
  seu <- subset(seu, subset = nFeature_RNA >= qc_min_features &
                      nFeature_RNA <= qc_max_features &
                      percent.mt <= qc_max_mt)

  # Normalization & Clustering
  if (use_glmGamPoi) {
    seu <- SCTransform(seu, vars.to.regress = vars_to_regress, method = "glmGamPoi", verbose = FALSE)
  } else {
    seu <- SCTransform(seu, vars.to.regress = vars_to_regress, verbose = FALSE)
  }

  seu <- RunPCA(seu, npcs = n_pcs, verbose = FALSE)
  seu <- FindNeighbors(seu, dims = 1:n_pcs, verbose = FALSE)
  seu <- FindClusters(seu, resolution = cluster_resolution, verbose = FALSE)
  seu <- RunUMAP(seu, dims = 1:n_pcs, verbose = FALSE)

  seu
}

# =========================
# SECTION 3: Make SoupX channel (FIXED SUBSETTING)
# =========================
build_soup_channel <- function(tod, toc, seu) {
  # [FIX] Intersect barcodes first
  common_cells <- intersect(colnames(toc), colnames(seu))
  
  if (length(common_cells) < 100) {
    stop("Critical Error: Low overlap between SoupX toc and Seurat cells: ", length(common_cells))
  }

  # [FIX] Subset the count matrix BEFORE creating SoupChannel
  # This prevents the 'incorrect dimensions' and 'invalid cluster' errors
  toc_subset <- toc[, common_cells]
  
  # Create object with subsetted matrix
  sc <- SoupChannel(tod = tod, toc = toc_subset)

  # Get metadata
  clusters <- seu$seurat_clusters
  names(clusters) <- colnames(seu)
  
  umap <- Embeddings(seu, "umap")

  # Attach clusters + 2D embedding (now guaranteed to match)
  sc <- setClusters(sc, clusters[common_cells])
  sc <- setDR(sc, umap[common_cells, ])

  list(sc = sc, common_cells = common_cells)
}

# =========================
# SECTION 4: Estimate contamination + correct counts (FIXED REDUNDANCY)
# =========================
estimate_and_correct <- function(sc, round_to_int = TRUE) {
  # [FIX] Check if soup profile exists before recalculating
  # This prevents the "rowSums array dimensions" crash
  if (is.null(sc$soupProfile)) {
    message("   Calculating soup profile...")
    sc <- estimateSoup(sc)
  } else {
    message("   Soup profile already exists. Skipping explicit estimateSoup().")
  }
  
  message("   Running autoEstCont...")
  sc <- autoEstCont(sc, forceAccept = TRUE)
  
  corrected <- adjustCounts(sc, roundToInt = round_to_int)
  list(sc = sc, corrected = corrected)
}

# =========================
# SECTION 5: Save outputs (FIXED DATAFRAME ERROR)
# =========================
save_outputs <- function(sample_id, out_dir, seu_pre, sc, corrected) {
  # [FIX] Use rownames(sc$metaData) to get barcodes reliably
  rho_vals <- sc$metaData$rho
  barcodes <- rownames(sc$metaData)
  
  rho_df <- data.frame(barcode = barcodes, rho = as.numeric(rho_vals))
  write.csv(rho_df, file.path(out_dir, paste0(sample_id, "_rho_per_cell.csv")), row.names = FALSE)

  # Histogram
  p_rho <- ggplot(rho_df, aes(x = rho)) +
    geom_histogram(bins = 50, fill = "steelblue", color = "white") +
    ggtitle(paste0(sample_id, " rho (contamination fraction)")) +
    xlab("rho") + ylab("Number of cells") + theme_minimal()
  ggsave(file.path(out_dir, paste0(sample_id, "_rho_hist.png")), p_rho, width = 7, height = 4, dpi = 150)

  # Build corrected Seurat object
  meta <- seu_pre@meta.data
  # [FIX] Ensure meta matches corrected counts
  meta <- meta[intersect(rownames(meta), colnames(corrected)), , drop = FALSE]

  seu_corr <- CreateSeuratObject(
    counts = corrected[, rownames(meta), drop = FALSE],
    meta.data = meta,
    project = paste0(sample_id, "_SoupX")
  )
  seu_corr[["percent.mt"]] <- PercentageFeatureSet(seu_corr, pattern = "^MT-")

  save_qc_plot(
    seu_corr,
    out_path = file.path(out_dir, paste0(sample_id, "_qc_after_SoupX.png")),
    title_text = paste0(sample_id, " QC (after SoupX)")
  )

  # Save RDS objects
  saveRDS(corrected, file.path(out_dir, paste0(sample_id, "_SoupX_corrected_counts.rds")))
  saveRDS(seu_corr,  file.path(out_dir, paste0(sample_id, "_SoupX_seurat.rds")))

  # Return numeric summary
  data.frame(
    sample_id = sample_id,
    n_cells_used_for_clusters = ncol(seu_pre),
    n_cells_corrected = ncol(corrected),
    rho_median = median(rho_vals, na.rm = TRUE),
    rho_mean = mean(rho_vals, na.rm = TRUE),
    rho_p95 = as.numeric(quantile(rho_vals, 0.95, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
}

# =========================
# MAIN: One sample pipeline
# =========================
run_soupx_one_sample <- function(sample_dir, sample_id, out_base_dir) {
  out_dir <- file.path(out_base_dir, sample_id)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  message("\n==============================")
  message("Sample: ", sample_id)
  message("Dir:    ", sample_dir)

  mats <- run_step("Load Cell Ranger H5 matrices", {
    load_cellranger_h5(sample_dir)
  })
  message("   Raw H5:      ", mats$raw_h5)
  message("   Filtered H5: ", mats$filtered_h5)

  seu <- run_step("Build Seurat (QC + SCTransform + clustering + UMAP)", {
    build_seurat_for_soupx(mats$toc, sample_id, out_dir)
  })

  soup <- run_step("Build SoupChannel + attach clusters/UMAP", {
    build_soup_channel(mats$tod, mats$toc, seu)
  })

  est <- run_step("Estimate contamination + adjust counts", {
    estimate_and_correct(soup$sc, round_to_int = TRUE)
  })

  smry <- run_step("Save outputs (rho, plots, corrected Seurat)", {
    save_outputs(sample_id, out_dir, seu, est$sc, est$corrected)
  })

  message("✅ ", sample_id, ": n_cells=", smry$n_cells_corrected,
          ", rho_median=", sprintf("%.3f", smry$rho_median),
          ", rho_p95=", sprintf("%.3f", smry$rho_p95))

  smry
}

# =========================
# RUN ALL SAMPLES
# =========================
# Get all subdirectories in base_dir
sample_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
# Exclude the output directory itself if it exists inside base_dir
sample_dirs <- sample_dirs[basename(sample_dirs) != basename(out_base)]

all_summaries <- list()

for (sd in sample_dirs) {
  sample_id <- basename(sd)
  
  sm <- tryCatch(
    run_soupx_one_sample(sd, sample_id, out_base),
    error = function(e) {
      message("⚠️ Skipping sample ", sample_id, " due to error: ", conditionMessage(e))
      NULL
    }
  )
  if (!is.null(sm)) all_summaries[[sample_id]] <- sm
}

# Final Summary CSV and Plot
summary_df <- do.call(rbind, all_summaries)
if (!is.null(summary_df) && nrow(summary_df) > 0) {
  write.csv(summary_df, file.path(out_base, "SoupX_summary_by_sample.csv"), row.names = FALSE)

  p_bar <- ggplot(summary_df, aes(x = reorder(sample_id, rho_median), y = rho_median)) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    ylab("Median rho per sample") + xlab("Sample") +
    ggtitle("SoupX contamination (rho) by sample (median)") +
    theme_minimal()
  
  ggsave(file.path(out_base, "SoupX_rho_median_barplot.png"), p_bar, width = 8, height = 6, dpi = 150)
}

message("\nDONE. Outputs written to: ", out_base)