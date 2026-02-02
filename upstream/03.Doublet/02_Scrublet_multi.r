suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(reticulate)
})

set.seed(1234)
options(future.globals.maxSize = 8 * 1024^3) # 8GB limit

# ==============================================================================
# USER SETTINGS
# ==============================================================================
base_dir <- "/mnt/18T/chibao/pbmc/data"
soupx_base_dir <- file.path(base_dir, "SoupX_outputs")
scrublet_out_base <- file.path(base_dir, "Scrublet_outputs")
dir.create(scrublet_out_base, recursive = TRUE, showWarnings = FALSE)

target_env <- "sc_bao"

# ==============================================================================
# UTIL: Step Wrapper (For nice logging)
# ==============================================================================
run_step <- function(step_name, expr) {
  message("\n--- STEP: ", step_name, " ---")
  tryCatch(expr, error = function(e) {
    message("❌ ERROR in step: ", step_name)
    stop(e)
  })
}

# ==============================================================================
# SECTION 1: Setup Python (Run ONCE)
# ==============================================================================
setup_python_env <- function(env_name) {
  message("--- Configuring Python Environment ---")
  tryCatch({
    reticulate::use_condaenv(env_name, required = TRUE)
  }, error = function(e) {
    message("⚠️ Python already initialized. Verifying active environment...")
  })
  
  # Check packages
  required <- c("scrublet", "numpy", "scipy", "scikit-learn", "matplotlib")
  missing <- c()
  for (p in required) {
    if (!reticulate::py_module_available(p)) missing <- c(missing, p)
  }
  
  if (length(missing) > 0) {
    message("Installing missing packages: ", paste(missing, collapse=", "))
    reticulate::py_install(missing, envname = env_name, pip = TRUE)
  }
  
  reticulate::py_run_string("import matplotlib; matplotlib.use('Agg')")
  message("Python ready: ", reticulate::py_config()$python)
}

# ==============================================================================
# SECTION 2: Load Data
# ==============================================================================
load_soupx_seurat <- function(sample_id, soupx_base) {
  rds_path <- file.path(soupx_base, sample_id, paste0(sample_id, "_SoupX_seurat.rds"))
  if (!file.exists(rds_path)) stop("File not found: ", rds_path)
  
  seu <- readRDS(rds_path)
  DefaultAssay(seu) <- "RNA"
  
  # Get counts safely
  if ("counts" %in% Layers(seu[["RNA"]])) {
    counts <- GetAssayData(seu, assay = "RNA", layer = "counts")
  } else {
    counts <- GetAssayData(seu, assay = "RNA", slot = "counts")
  }
  list(seu = seu, counts = counts)
}

# ==============================================================================
# SECTION 3: Estimate Doublet Rate
# ==============================================================================
calc_expected_rate <- function(n_cells) {
  # 10x Genomics Multiplet Rate Table
  tab <- data.frame(
    recovered = c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000),
    rate      = c(0.004,0.008,0.016,0.023,0.031,0.039,0.046,0.054,0.061,0.069,0.076)
  )
  n_clamped <- max(min(n_cells, max(tab$recovered)), min(tab$recovered))
  rate <- approx(tab$recovered, tab$rate, xout = n_clamped)$y
  rate
}

# ==============================================================================
# SECTION 4: Python Bridge (Export & Run)
# ==============================================================================
run_scrublet_py <- function(counts, expected_rate, out_dir, sample_id) {
  # 1. Export Matrix
  tmp_dir <- file.path(out_dir, "scrublet_tmp")
  dir.create(tmp_dir, showWarnings = FALSE)
  
  mtx_path <- file.path(tmp_dir, "matrix.mtx")
  Matrix::writeMM(counts, file = mtx_path)
  
  # 2. Run Python
  scipy_io <- reticulate::import("scipy.io", convert = FALSE) # [CRITICAL]
  scr      <- reticulate::import("scrublet")
  plt      <- reticulate::import("matplotlib.pyplot")
  
  # Load & Transpose (R is Genes x Cells, Python needs Cells x Genes)
  X_genes_cells <- scipy_io$mmread(mtx_path)
  X_cells_genes <- X_genes_cells$T$tocsr()
  
  scrub <- scr$Scrublet(X_cells_genes, expected_doublet_rate = expected_rate)
  res   <- scrub$scrub_doublets(
    min_counts = 2L, min_cells = 3L, 
    min_gene_variability_pctl = 85L, n_prin_comps = 30L
  )
  
  # Save Histogram
  scrub$plot_histogram()
  plt$savefig(file.path(out_dir, paste0(sample_id, "_Scrublet_histogram.png")), dpi=150L)
  plt$close("all")
  
  # Cleanup
  unlink(tmp_dir, recursive = TRUE)
  
  # Return Results
  list(
    scores = as.numeric(res[[1]]),
    preds  = as.logical(res[[2]]),
    thresh = as.numeric(scrub$threshold_)
  )
}

# ==============================================================================
# SECTION 5: Merge & Visualize
# ==============================================================================
merge_and_visualize <- function(seu, scrub_res, expected_rate, sample_id, out_dir) {
  # Add Metadata
  seu$Scrublet_score <- scrub_res$scores
  seu$Scrublet_pred  <- ifelse(scrub_res$preds, "Doublet", "Singlet")
  seu$Scrublet_threshold <- scrub_res$thresh
  seu$Scrublet_expected_rate <- expected_rate
  
  # Quick Visualization
  if (!"umap" %in% names(seu@reductions)) {
    seu <- SCTransform(seu, verbose = FALSE)
    seu <- RunPCA(seu, verbose = FALSE)
    seu <- RunUMAP(seu, dims = 1:30, verbose = FALSE)
  }
  
  p1 <- DimPlot(seu, reduction = "umap", group.by = "Scrublet_pred", cols = c("Singlet"="grey", "Doublet"="red")) +
    ggtitle(paste0(sample_id, " Prediction"))
  ggsave(file.path(out_dir, paste0(sample_id, "_UMAP_Scrublet_pred.png")), p1, width = 7, height = 5)
  
  p2 <- FeaturePlot(seu, features = "Scrublet_score") + scale_color_viridis_c(option = "magma")
  ggsave(file.path(out_dir, paste0(sample_id, "_UMAP_Scrublet_score.png")), p2, width = 7, height = 5)
  
  seu
}

# ==============================================================================
# SECTION 6: Save Outputs
# ==============================================================================
save_results <- function(seu, sample_id, out_dir) {
  # 1. Full Object
  saveRDS(seu, file.path(out_dir, paste0(sample_id, "_SoupX_Scrublet_annotated_seurat.rds")))
  
  # 2. Singlets Only
  seu_singlets <- subset(seu, subset = Scrublet_pred == "Singlet")
  saveRDS(seu_singlets, file.path(out_dir, paste0(sample_id, "_SoupX_Scrublet_singlets_seurat.rds")))
  
  # 3. CSV Stats
  write.csv(
    data.frame(barcode = colnames(seu), score = seu$Scrublet_score, pred = seu$Scrublet_pred),
    file.path(out_dir, paste0(sample_id, "_Scrublet_barcodes.csv")), 
    row.names = FALSE
  )
  
  # Return Summary Row
  data.frame(
    sample_id = sample_id,
    n_total = ncol(seu),
    n_singlets = ncol(seu_singlets),
    doublet_rate = round((1 - ncol(seu_singlets)/ncol(seu))*100, 2),
    threshold = seu$Scrublet_threshold[1]
  )
}

# ==============================================================================
# MAIN PIPELINE WRAPPER
# ==============================================================================
process_one_sample_scrublet <- function(sample_id, soupx_base, out_base) {
  out_dir <- file.path(out_base, sample_id)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  message("\n=== Processing: ", sample_id, " ===")
  
  # Step 1: Load
  data <- run_step("Load SoupX Data", { load_soupx_seurat(sample_id, soupx_base) })
  
  # Step 2: Calculate Params
  rate <- run_step("Calc Expected Rate", { calc_expected_rate(ncol(data$seu)) })
  message("   Expected Doublet Rate: ", round(rate, 4))
  
  # Step 3: Run Scrublet (Python)
  scrub_res <- run_step("Run Scrublet (Python)", { 
    run_scrublet_py(data$counts, rate, out_dir, sample_id) 
  })
  message("   Detected ", sum(scrub_res$preds), " doublets.")
  
  # Step 4: Merge & Plot
  seu_final <- run_step("Merge & Visualize", { 
    merge_and_visualize(data$seu, scrub_res, rate, sample_id, out_dir) 
  })
  
  # Step 5: Save
  stats <- run_step("Save Outputs", { save_results(seu_final, sample_id, out_dir) })
  
  stats
}

# ==============================================================================
# EXECUTION LOOP
# ==============================================================================
# 1. Setup Python
setup_python_env(target_env)

# 2. Find Samples
sample_dirs <- list.dirs(soupx_base_dir, full.names = TRUE, recursive = FALSE)
sample_dirs <- sample_dirs[basename(sample_dirs) != basename(scrublet_out_base)]

all_stats <- list()

# 3. Loop
for (sd in sample_dirs) {
  sid <- basename(sd)
  
  res <- tryCatch({
    process_one_sample_scrublet(sid, soupx_base_dir, scrublet_out_base)
  }, error = function(e) {
    message("⚠️ Skipped ", sid, ": ", conditionMessage(e))
    NULL
  })
  
  if (!is.null(res)) all_stats[[sid]] <- res
}

# 4. Final Summary
if (length(all_stats) > 0) {
  final_df <- do.call(rbind, all_stats)
  print(final_df)
  write.csv(final_df, file.path(scrublet_out_base, "Scrublet_Batch_Summary.csv"), row.names = FALSE)
  message("\n✅ Batch Processing Complete.")
}