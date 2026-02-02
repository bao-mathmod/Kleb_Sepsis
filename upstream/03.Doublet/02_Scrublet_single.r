suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(reticulate)
})

set.seed(1234)

# ==============================================================================
# USER SETTINGS (Edit this section)
# ==============================================================================
# The ID of the sample you want to test
target_sample_id <- "PBMC_data" 

# Directory containing your SoupX outputs
# Structure assumed: base_dir/SoupX_outputs/Sample_Name_01/Sample_Name_01_SoupX_seurat.rds
base_dir <- "/mnt/18T/chibao/pbmc/data"
soupx_dir <- file.path(base_dir, "SoupX_outputs")

# Where to save Scrublet results
out_dir <- file.path(base_dir, "Scrublet_outputs", target_sample_id)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Processing parameters
dims_use <- 1:30
scrublet_env_name <- "scrublet_env"

message("Processing Sample: ", target_sample_id)
message("Output Directory:  ", out_dir)

# ==============================================================================
# STEP 1: Configure Python & Scrublet Environment
# ==============================================================================
message("\n--- STEP 1: Setting up Python Environment ---")

# # Check if conda is available
# if (reticulate::conda_binary(conda = "auto") != "") {
#   # Create env if missing
#   if (!scrublet_env_name %in% reticulate::conda_list()$name) {
#     message("Creating conda env: ", scrublet_env_name)
#     reticulate::conda_create(scrublet_env_name, python_version = "3.10")
#     # Install packages
#     reticulate::py_install(c("scrublet", "numpy", "scipy", "scikit-learn", "matplotlib"), 
#                            envname = scrublet_env_name, pip = TRUE)
#   }
#   # Use the environment
#   reticulate::use_condaenv(scrublet_env_name, required = TRUE)
# } else {
#   message("Conda not found. Attempting to use default Python...")
#   # Fallback: install into default python
#   reticulate::py_install(c("scrublet", "numpy", "scipy", "scikit-learn", "matplotlib"), pip = TRUE)
# }

# # Force matplotlib to use a non-interactive backend (prevents crashes on servers)
# reticulate::py_run_string("import matplotlib; matplotlib.use('Agg')")
# message("Python configuration ready.")
# print(reticulate::py_config())

message("\n--- STEP 1: Checking Python Environment 'sc_bao' ---")

# 1. Define the environment and packages we need
target_env <- "sc_bao"
required_pkgs <- c("scrublet", "numpy", "scipy", "scikit-learn", "matplotlib")

# 2. Tell Reticulate to use your specific environment
#    Note: Since you started R inside the environment, reticulate usually picks it up automatically.
#    However, this command explicitly forces it to ensure consistency.
tryCatch({
  reticulate::use_condaenv(target_env, required = TRUE)
}, error = function(e) {
  # If Python was already initialized (e.g. by another library), this might warn.
  # We print the warning but continue, as likely it initialized to the right place.
  message("⚠️ Notice: Python already initialized. Verifying active environment...")
})

# 3. Verify which Python R is actually using
py_config <- reticulate::py_config()
message("Active Python Path: ", py_config$python)
message("Active Virtual Env: ", Sys.getenv("CONDA_PREFIX"))

# 4. Check for missing packages and install ONLY if needed
missing_pkgs <- c()

for (pkg in required_pkgs) {
  # Check if the package can be imported
  if (!reticulate::py_module_available(pkg)) {
    message("❌ Missing package: ", pkg)
    missing_pkgs <- c(missing_pkgs, pkg)
  } else {
    message("✅ Found package: ", pkg)
  }
}

# 5. Install missing packages into 'sc_bao'
if (length(missing_pkgs) > 0) {
  message("Installing missing packages: ", paste(missing_pkgs, collapse = ", "))
  
  # We use pip = TRUE to ensure compatibility with Scrublet
  reticulate::py_install(missing_pkgs, envname = target_env, pip = TRUE)
  
  message("Installation complete.")
} else {
  message("All required Python packages are present.")
}

# 6. Configure Matplotlib (Headless mode)
reticulate::py_run_string("import matplotlib; matplotlib.use('Agg')")
message("Python configuration ready.")

# ==============================================================================
# STEP 2: Load Data & Prepare Counts
# ==============================================================================
message("\n--- STEP 2: Loading SoupX Seurat Object ---")

# Construct path to the SoupX object
seu_path <- file.path(soupx_dir, target_sample_id, paste0(target_sample_id, "_SoupX_seurat.rds"))

if (!file.exists(seu_path)) stop("Could not find input file: ", seu_path)

seu <- readRDS(seu_path)
DefaultAssay(seu) <- "RNA"

# Extract Counts (Seurat v5 friendly)
if ("counts" %in% Layers(seu[["RNA"]])) {
  counts <- GetAssayData(seu, assay = "RNA", layer = "counts")
} else {
  counts <- GetAssayData(seu, assay = "RNA", slot = "counts")
}

message("Loaded dimensions: ", paste(dim(counts), collapse = " x "))


# ==============================================================================
# STEP 3: Calculate Expected Doublet Rate
# ==============================================================================
message("\n--- STEP 3: Calculating Expected Doublet Rate ---")

n_cells <- ncol(seu)

# 10x Genomics Multiplet Rate Table (linear interpolation)
# Source: 10x User Guide
tab <- data.frame(
  recovered = c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000),
  rate      = c(0.004,0.008,0.016,0.023,0.031,0.039,0.046,0.054,0.061,0.069,0.076)
)
n_clamped <- max(min(n_cells, max(tab$recovered)), min(tab$recovered))
expected_rate <- approx(tab$recovered, tab$rate, xout = n_clamped)$y

message("Number of cells: ", n_cells)
message("Estimated 10x Doublet Rate: ", round(expected_rate, 4), " (", round(expected_rate*100, 2), "%)")


# ==============================================================================
# STEP 4: Export Matrix for Python
# ==============================================================================
message("\n--- STEP 4: Transferring Matrix to Python ---")

# Create a temporary directory for data transfer
tmp_dir <- file.path(out_dir, "scrublet_tmp")
dir.create(tmp_dir, showWarnings = FALSE)

mtx_path      <- file.path(tmp_dir, "matrix.mtx")
genes_path    <- file.path(tmp_dir, "genes.tsv")
barcodes_path <- file.path(tmp_dir, "barcodes.tsv")

message("Writing MatrixMarket files to: ", tmp_dir)
Matrix::writeMM(counts, file = mtx_path)
writeLines(rownames(counts), con = genes_path)
writeLines(colnames(counts), con = barcodes_path)


# ==============================================================================
# STEP 5: Run Scrublet (Python Code)
# ==============================================================================
message("\n--- STEP 5: Running Scrublet in Python ---")

# Import Python libraries
# [FIX] We set convert = FALSE for scipy.io. 
# This ensures the loaded matrix stays a "Python Object" so we can run .T (transpose) on it.
scipy_io     <- reticulate::import("scipy.io", convert = FALSE) 
scr          <- reticulate::import("scrublet")
plt          <- reticulate::import("matplotlib.pyplot")

# Load data in Python
# Now 'X_genes_cells' will be a pointer to a Python object, NOT an R matrix.
X_genes_cells <- scipy_io$mmread(mtx_path)

# Now .T (transpose) and .tocsr() will work because it's still treated as Python code
X_cells_genes <- X_genes_cells$T$tocsr()

# Initialize Scrublet object
message("Initializing Scrublet...")
scrub <- scr$Scrublet(X_cells_genes, expected_doublet_rate = expected_rate)

# Run Doublet Detection
message("Scrubbing doublets...")
res <- scrub$scrub_doublets(
  min_counts = 2L,
  min_cells = 3L,
  min_gene_variability_pctl = 85L,
  n_prin_comps = 30L
)

# Extract results
# 'res' is a Python tuple. We access elements using standard R list syntax [[ ]]
doublet_scores     <- as.numeric(res[[1]])
predicted_doublets <- as.logical(res[[2]])
used_threshold     <- as.numeric(scrub$threshold_)

message("Scrublet Finished.")
message("Threshold chosen: ", round(used_threshold, 4))
message("Doublets detected: ", sum(predicted_doublets))


# ==============================================================================
# STEP 6: Save Scrublet Plots (Python)
# ==============================================================================
message("\n--- STEP 6: Saving Scrublet Histogram ---")

scrub$plot_histogram()
hist_filename <- file.path(out_dir, paste0(target_sample_id, "_Scrublet_histogram.png"))
plt$savefig(hist_filename, dpi = as.integer(150), bbox_inches = "tight")
plt$close("all") # Clear plot memory

message("Histogram saved to: ", hist_filename)


# ==============================================================================
# STEP 7: Import Results to Seurat & Cleanup
# ==============================================================================
message("\n--- STEP 7: Merging Results into Seurat ---")

# Verify barcodes align
# The matrix read by Python was in the same order as written by R, so safe to assume alignment.
if (length(doublet_scores) != ncol(seu)) {
  stop("Error: Number of scores returned by Python does not match Seurat object.")
}

# Add Metadata
seu$Scrublet_score         <- doublet_scores
seu$Scrublet_pred          <- ifelse(predicted_doublets, "Doublet", "Singlet")
seu$Scrublet_threshold     <- used_threshold
seu$Scrublet_expected_rate <- expected_rate

# Clean up temp files
unlink(tmp_dir, recursive = TRUE)


# ==============================================================================
# STEP 8: Visualization (R Plots)
# ==============================================================================
message("\n--- STEP 8: Generating QC Plots ---")
options(future.globals.maxSize = 8 * 1024^3)  # 8GB
# Check if UMAP exists. If not, run a quick SCT/PCA/UMAP
if (!"umap" %in% names(seu@reductions)) {
  message("No UMAP found. Running quick preprocessing...")
  seu <- SCTransform(seu, verbose = FALSE)
  seu <- RunPCA(seu, verbose = FALSE)
  seu <- RunUMAP(seu, dims = 1:30, verbose = FALSE)
}

# 1. UMAP: Prediction
p1 <- DimPlot(seu, reduction = "umap", group.by = "Scrublet_pred", cols = c("Singlet" = "grey", "Doublet" = "red")) +
  ggtitle(paste0(target_sample_id, " (Scrublet Prediction)"))
ggsave(file.path(out_dir, paste0(target_sample_id, "_UMAP_Scrublet_pred.png")), p1, width = 7, height = 5)

# 2. UMAP: Score (Continuous)
p2 <- FeaturePlot(seu, features = "Scrublet_score") +
  scale_color_viridis_c(option = "magma") +
  ggtitle(paste0(target_sample_id, " (Scrublet Score)"))
ggsave(file.path(out_dir, paste0(target_sample_id, "_UMAP_Scrublet_score.png")), p2, width = 7, height = 5)

# 3. Violin Plot: QC Metrics
p3 <- VlnPlot(seu, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), 
              group.by = "Scrublet_pred", pt.size = 0, ncol = 3)
ggsave(file.path(out_dir, paste0(target_sample_id, "_QC_by_Scrublet.png")), p3, width = 10, height = 4)


# ==============================================================================
# STEP 9: Save Outputs
# ==============================================================================
message("\n--- STEP 9: Saving Final Objects ---")

# 1. Save Full Object (Annotated)
saveRDS(seu, file.path(out_dir, paste0(target_sample_id, "_SoupX_Scrublet_annotated_seurat.rds")))

# 2. Save CSV Summary
anno_df <- data.frame(
  barcode = colnames(seu),
  Scrublet_score = seu$Scrublet_score,
  Scrublet_pred = seu$Scrublet_pred,
  stringsAsFactors = FALSE
)
write.csv(anno_df, file.path(out_dir, paste0(target_sample_id, "_Scrublet_barcodes.csv")), row.names = FALSE)

# 3. Save Singlets Only (Clean Object)
seu_singlets <- subset(seu, subset = Scrublet_pred == "Singlet")
saveRDS(seu_singlets, file.path(out_dir, paste0(target_sample_id, "_SoupX_Scrublet_singlets_seurat.rds")))

message("✅ DONE! Processing complete for: ", target_sample_id)
message("   Singlets saved: ", ncol(seu_singlets))
message("   Doublets removed: ", ncol(seu) - ncol(seu_singlets))