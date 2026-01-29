suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(SoupX)
  library(ggplot2)
  library(dplyr)
})

# =========================
# USER SETTINGS (MAD VERSION)
# =========================
base_dir <- "/mnt/18T/chibao/pbmc/data_3"
out_base <- file.path(base_dir, "SoupX_outputs_MAD")
dir.create(out_base, showWarnings = FALSE, recursive = TRUE)

# MAD Settings (Dynamic QC)
nmads <- 3               # Number of MADs to define outliers (standard is 3 or 4)
min_cells_feature <- 200 # Hard floor for features (optional safety net)

# Clustering (only to provide clusters + UMAP for SoupX)
n_pcs <- 30
cluster_resolution <- 0.8

# SCTransform options
use_glmGamPoi <- TRUE
vars_to_regress <- "percent.mt"

set.seed(1234)

# =========================
# UTIL: step wrapper
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
  if (length(fil) == 0) fil <- list.files(sample_dir, pattern = "filtered_featured.*\\.h5$", full.names = TRUE, recursive = TRUE)
  if (length(raw) == 0 || length(fil) == 0) return(list(raw = NA_character_, filtered = NA_character_))
  
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

# Helper: Calculate MAD thresholds
# We log-transform features/counts because they are usually log-normally distributed
get_mad_bounds <- function(x, nmads, lower = TRUE, log_scale = FALSE) {
  if (log_scale) x <- log1p(x)
  med_val <- median(x, na.rm = TRUE)
  mad_val <- mad(x, constant = 1.4826, na.rm = TRUE)
  
  upper_bound <- med_val + (nmads * mad_val)
  lower_bound <- med_val - (nmads * mad_val)
  
  if (log_scale) {
    upper_bound <- expm1(upper_bound)
    lower_bound <- expm1(lower_bound)
  }
  
  list(lower = if(lower) lower_bound else -Inf, upper = upper_bound)
}

# =========================
# SECTION 1: Load matrices
# =========================
load_cellranger_h5 <- function(sample_dir) {
  h5s <- find_h5_pair(sample_dir)
  if (is.na(h5s$raw) || is.na(h5s$filtered)) stop("Missing raw/filtered .h5 in: ", sample_dir)
  tod <- read_10x_h5_gex(h5s$raw)
  toc <- read_10x_h5_gex(h5s$filtered)
  list(tod = tod, toc = toc, raw_h5 = h5s$raw, filtered_h5 = h5s$filtered)
}

# =========================
# SECTION 2: Build Seurat + QC (MAD)
# =========================
build_seurat_for_soupx <- function(toc, sample_id, out_dir) {
  seu <- CreateSeuratObject(counts = toc, project = sample_id, min.cells = 3, min.features = 0)
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

  # 1. Save Pre-QC Plot (Raw State)
  save_qc_plot(seu, file.path(out_dir, paste0(sample_id, "_qc_before_MAD.png")), paste0(sample_id, " Raw Data"))

  # 2. Calculate MAD Thresholds
  # Features (Low end is bad)
  feat_lims <- get_mad_bounds(seu$nFeature_RNA, nmads = nmads, lower = TRUE, log_scale = TRUE)
  # Mito (High end is bad)
  mt_lims   <- get_mad_bounds(seu$percent.mt, nmads = nmads, lower = FALSE, log_scale = FALSE)
  
  # Safety: Enforce a hard floor for features (e.g., don't keep cells with 10 genes even if MAD says ok)
  cut_feature_low <- max(feat_lims$lower, min_cells_feature)
  cut_mt_high     <- mt_lims$upper

  message(sprintf("   MAD Cutoffs -> nFeature < %.0f | percent.mt > %.2f", cut_feature_low, cut_mt_high))

  # 3. Apply Filter
  seu <- subset(seu, subset = nFeature_RNA >= cut_feature_low & percent.mt <= cut_mt_high)
  
  # 4. Save Post-QC Plot (To visualize what MAD kept)
  save_qc_plot(seu, file.path(out_dir, paste0(sample_id, "_qc_after_MAD.png")), paste0(sample_id, " After MAD Filter"))

  # 5. SCTransform & Clustering
  if (use_glmGamPoi) {
    seu <- SCTransform(seu, vars.to.regress = vars_to_regress, method = "glmGamPoi", verbose = FALSE)
  } else {
    seu <- SCTransform(seu, vars.to.regress = vars_to_regress, verbose = FALSE)
  }

  seu <- RunPCA(seu, npcs = n_pcs, verbose = FALSE)
  seu <- FindNeighbors(seu, dims = 1:n_pcs, verbose = FALSE)
  seu <- FindClusters(seu, resolution = cluster_resolution, verbose = FALSE)
  seu <- RunUMAP(seu, dims = 1:n_pcs, verbose = FALSE)

  list(seu = seu, cut_feature = cut_feature_low, cut_mt = cut_mt_high)
}

# =========================
# SECTION 3: Make SoupX channel
# =========================
build_soup_channel <- function(tod, toc, seu) {
  common_cells <- intersect(colnames(toc), colnames(seu))
  if (length(common_cells) < 100) stop("Critical Error: Low overlap after QC.")

  toc_subset <- toc[, common_cells]
  sc <- SoupChannel(tod = tod, toc = toc_subset)

  clusters <- seu$seurat_clusters
  names(clusters) <- colnames(seu)
  umap <- Embeddings(seu, "umap")

  sc <- setClusters(sc, clusters[common_cells])
  sc <- setDR(sc, umap[common_cells, ])

  list(sc = sc, common_cells = common_cells)
}

# =========================
# SECTION 4: Estimate & Correct
# =========================
estimate_and_correct <- function(sc, round_to_int = TRUE) {
  if (is.null(sc$soupProfile)) {
    message("   Calculating soup profile...")
    sc <- estimateSoup(sc)
  }
  message("   Running autoEstCont...")
  sc <- autoEstCont(sc, forceAccept = TRUE)
  corrected <- adjustCounts(sc, roundToInt = round_to_int)
  list(sc = sc, corrected = corrected)
}

# =========================
# SECTION 5: Save outputs
# =========================
save_outputs <- function(sample_id, out_dir, seu_pre, sc, corrected, mad_stats) {
  # Rho CSV
  rho_vals <- sc$metaData$rho
  rho_df <- data.frame(barcode = rownames(sc$metaData), rho = as.numeric(rho_vals))
  write.csv(rho_df, file.path(out_dir, paste0(sample_id, "_rho_per_cell.csv")), row.names = FALSE)

  # Rho Plot
  p_rho <- ggplot(rho_df, aes(x = rho)) +
    geom_histogram(bins = 50, fill = "forestgreen", color = "white") +
    ggtitle(paste0(sample_id, " Contamination (MAD Filtered)")) + theme_minimal()
  ggsave(file.path(out_dir, paste0(sample_id, "_rho_hist.png")), p_rho, width = 7, height = 4)

  # Corrected Seurat
  meta <- seu_pre@meta.data
  meta <- meta[intersect(rownames(meta), colnames(corrected)), , drop = FALSE]

  seu_corr <- CreateSeuratObject(counts = corrected[, rownames(meta)], meta.data = meta, project = paste0(sample_id, "_SoupX_MAD"))
  seu_corr[["percent.mt"]] <- PercentageFeatureSet(seu_corr, pattern = "^MT-")

  saveRDS(corrected, file.path(out_dir, paste0(sample_id, "_SoupX_corrected_counts.rds")))
  saveRDS(seu_corr,  file.path(out_dir, paste0(sample_id, "_SoupX_seurat.rds")))

  # Summary Stats
  data.frame(
    sample_id = sample_id,
    mad_cutoff_feature = mad_stats$cut_feature,
    mad_cutoff_mt = mad_stats$cut_mt,
    n_cells_final = ncol(corrected),
    rho_median = median(rho_vals, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

# =========================
# MAIN LOOP
# =========================
run_soupx_mad <- function(sample_dir, sample_id, out_base_dir) {
  out_dir <- file.path(out_base_dir, sample_id)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  message("\n==============================")
  message("MAD Pipeline: ", sample_id)

  mats <- run_step("Load Data", { load_cellranger_h5(sample_dir) })

  # Modified Step: Returns list(seu, cut_feature, cut_mt)
  seu_res <- run_step("Seurat + MAD QC", { build_seurat_for_soupx(mats$toc, sample_id, out_dir) })
  seu <- seu_res$seu # Extract the object

  soup <- run_step("SoupChannel", { build_soup_channel(mats$tod, mats$toc, seu) })
  
  est <- run_step("Soup Correction", { estimate_and_correct(soup$sc) })

  smry <- run_step("Save Results", { save_outputs(sample_id, out_dir, seu, est$sc, est$corrected, seu_res) })

  message("✅ ", sample_id, ": rho_median=", sprintf("%.3f", smry$rho_median), " | MT Cutoff used: ", sprintf("%.2f", seu_res$cut_mt), "%")
  smry
}

# Run
sample_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
sample_dirs <- sample_dirs[basename(sample_dirs) != basename(out_base)]

all_summaries <- list()
for (sd in sample_dirs) {
  sample_id <- basename(sd)
  sm <- tryCatch(run_soupx_mad(sd, sample_id, out_base), error = function(e) { message("⚠️ Error: ", conditionMessage(e)); NULL })
  if (!is.null(sm)) all_summaries[[sample_id]] <- sm
}

summary_df <- do.call(rbind, all_summaries)
if (!is.null(summary_df)) write.csv(summary_df, file.path(out_base, "SoupX_MAD_summary.csv"), row.names = FALSE)
message("DONE.")