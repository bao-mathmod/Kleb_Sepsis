# ============================================================
# DoubletFinder on SoupX-corrected Seurat (Seurat v5.3.1)
# Input: *_SoupX_seurat.rds and *_SoupX_corrected_counts.rds
# Output: DF-annotated Seurat, singlet-only Seurat, plots + CSVs
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(remotes)
  library(DoubletFinder)
})


set.seed(1234)

# -------------------------
# 0) Paths
# -------------------------
sample_id   <- "PBMC_data"
seu_path    <- "/mnt/18T/chibao/pbmc/data_3/SoupX_outputs_MAD/PBMC_data/PBMC_data_SoupX_seurat.rds"
counts_path <- "/mnt/18T/chibao/pbmc/data_3/SoupX_outputs_MAD/PBMC_data/PBMC_data_SoupX_corrected_counts.rds"

out_dir <- file.path("/mnt/18T/chibao/pbmc/data_3/DoubletFinder_outputs", sample_id)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------
# 1) Load SoupX outputs
# -------------------------
seu <- readRDS(seu_path)
corrected_counts <- readRDS(counts_path)

# Optional: verify counts match Seurat RNA counts layer (sanity check)
# Seurat v5 stores data in layers. We'll compare dimensions + a few entries.
DefaultAssay(seu) <- "RNA"

rna_counts <- GetAssayData(seu, assay = "RNA", layer = "counts")
message("RNA counts dim: ", paste(dim(rna_counts), collapse=" x "),
        " | corrected_counts dim: ", paste(dim(corrected_counts), collapse=" x "))

# -------------------------
# 2) Seurat v5 "layers" safety
# -------------------------
# If you ever have multiple layers (e.g., from split/integration), JoinLayers can help.
# (In your SoupX object it's usually just one counts layer.)
if ("RNA" %in% names(seu@assays)) {
  # Check layer names if available
  # Layers() is the Seurat v5 way to see stored layers.
  if ("Layers" %in% ls(getNamespace("SeuratObject"))) {
    # not strictly necessary, but safe if you use split layers later
    # message("RNA layers: ", paste(Layers(seu[["RNA"]]), collapse=", "))
  }
}

# -------------------------
# 3) Preprocess (SCT -> PCA -> UMAP -> clustering)
#    DoubletFinder expects a "fully processed" object (normalization, PCA, etc.) :contentReference[oaicite:5]{index=5}
# -------------------------
options(future.globals.maxSize = 8 * 1024^3)  # 8GB
run_preprocess_sct <- function(seu, dims_use = 1:30, resolution = 0.8) {
  DefaultAssay(seu) <- "RNA"

  # SCT creates an SCT assay; use corrected counts in RNA as input.
  seu <- SCTransform(seu, verbose = FALSE)
  DefaultAssay(seu) <- "SCT"

  seu <- RunPCA(seu, npcs = max(dims_use), verbose = FALSE)
  seu <- RunUMAP(seu, dims = dims_use, reduction = "pca", verbose = FALSE)
  seu <- FindNeighbors(seu, dims = dims_use, reduction = "pca", verbose = FALSE)
  seu <- FindClusters(seu, resolution = resolution, verbose = FALSE)

  return(seu)
}

dims_use <- 1:30
seu <- run_preprocess_sct(seu, dims_use = dims_use, resolution = 0.8)

# Save a quick UMAP before DF
p_umap_pre <- DimPlot(seu, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle(paste0(sample_id, " (post-SoupX, pre-DoubletFinder clusters)"))
ggsave(file.path(out_dir, paste0(sample_id, "_UMAP_pre_DF_clusters.png")),
       p_umap_pre, width = 7, height = 5, dpi = 150)

# -------------------------
# 4) pK sweep (paramSweep -> summarizeSweep -> find.pK)
# -------------------------
# NOTE: with Seurat v5, DoubletFinder removed the _v3 naming and supports v5 now
# Main Purpose: Find optimal pK param for DoubletFinder
# For SCT workflow, set sct=TRUE.

sweep.res <- paramSweep(seu, PCs = dims_use, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

# Save sweep table + plot
write.csv(bcmvn, file.path(out_dir, paste0(sample_id, "_pK_sweep_BCmvn.csv")), row.names = FALSE)

bcmvn$pK_num <- as.numeric(as.character(bcmvn$pK))
p_pk <- ggplot(bcmvn, aes(x = pK_num, y = BCmetric)) +
  geom_line() + geom_point() +
  ggtitle(paste0(sample_id, " DoubletFinder pK sweep (BCmvn)")) +
  xlab("pK") + ylab("BCmvn metric")
ggsave(file.path(out_dir, paste0(sample_id, "_pK_sweep_plot.png")),
       p_pk, width = 7, height = 4, dpi = 150)

best_pK <- bcmvn$pK[which.max(bcmvn$BCmetric)]
best_pK <- as.numeric(as.character(best_pK))
message("Selected best pK = ", best_pK)

# -------------------------
# 5) Estimate expected doublets (nExp)
# -------------------------
# DoubletFinder recommends estimating nExp from platform loading densities and adjusting for homotypic doublets.

estimate_10x_multiplet_rate <- function(n_recovered) {
  # From 10x "Getting Started: Single Cell Immune Profiling" table (recovered vs multiplet rate).
  tab <- data.frame(
    recovered = c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000),
    rate      = c(0.004,0.008,0.016,0.023,0.031,0.039,0.046,0.054,0.061,0.069,0.076)
  ) # Source: https://cdn.10xgenomics.com/image/upload/v1660261285/support-documents/CG000361_GettingStartedImmuneProfiling_RevA.pdf
  # clamp + linear interpolation
  n <- max(min(n_recovered, max(tab$recovered)), min(tab$recovered))
  approx(tab$recovered, tab$rate, xout = n)$y
}

n_cells <- ncol(seu)
dbl_rate <- estimate_10x_multiplet_rate(n_cells)
nExp_poi <- round(dbl_rate * n_cells)

# Homotypic adjustment
homotypic.prop <- modelHomotypic(seu$seurat_clusters)
nExp_adj <- round(nExp_poi * (1 - homotypic.prop))

message("Cells=", n_cells,
        " | 10x multiplet rate~", round(dbl_rate*100, 2), "% ",
        " | nExp_poi=", nExp_poi,
        " | homotypic.prop~", round(homotypic.prop, 3),
        " | nExp_adj=", nExp_adj)

# Save expected doublet summary
exp_df <- data.frame(
  sample_id = sample_id,
  n_cells = n_cells,
  multiplet_rate_est = dbl_rate,
  nExp_poi = nExp_poi,
  homotypic_prop = homotypic.prop,
  nExp_adj = nExp_adj
)
write.csv(exp_df, file.path(out_dir, paste0(sample_id, "_expected_doublets_summary.csv")),
          row.names = FALSE)

# -------------------------
# 6) Run DoubletFinder (two-pass)
# -------------------------
# Pass 1: nExp_poi
seu <- doubletFinder(
  seu,
  PCs = dims_use,
  pN = 0.25,
  pK = best_pK,
  nExp = nExp_poi,
  reuse.pANN = NULL,   
  sct = TRUE
)

# Identify the pANN column created in pass 1 (pattern starts with "pANN_")
pANN_col <- tail(grep("^pANN_", colnames(seu@meta.data), value = TRUE), 1)
# pANN_col <- pANN_col[length(pANN_col)]  # take the newest one

# Pass 2: homotypic-adjusted nExp (reuse pANN from pass 1)
seu <- doubletFinder(
  seu,
  PCs = dims_use,
  pN = 0.25,
  pK = best_pK,
  nExp = nExp_adj,
  reuse.pANN = pANN_col, 
  sct = TRUE
)

# Find the newest DF classification column
df_class_cols <- grep("^DF.classifications_", colnames(seu@meta.data), value = TRUE)
df_class_col  <- df_class_cols[length(df_class_cols)]
message("Using DF classification column: ", df_class_col)

# -------------------------
# 7) QC plots (graphs) + counts (numbers)
# -------------------------
# Count predicted doublets/singlets (example numeric output)
df_counts <- table(seu@meta.data[[df_class_col]])
df_counts_df <- as.data.frame(df_counts)
colnames(df_counts_df) <- c("class", "n_cells")
write.csv(df_counts_df, file.path(out_dir, paste0(sample_id, "_DF_class_counts.csv")),
          row.names = FALSE)

# UMAP colored by DF class
p_umap_df <- DimPlot(seu, reduction = "umap", group.by = df_class_col) +
  ggtitle(paste0(sample_id, " DoubletFinder classification"))
ggsave(file.path(out_dir, paste0(sample_id, "_UMAP_DF_class.png")),
       p_umap_df, width = 7, height = 5, dpi = 150)

# Violin plots: typical QC metrics differ in doublets vs singlets
# (doublets often have higher nCount/nFeature)
p_vln <- VlnPlot(seu, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
                 group.by = df_class_col, pt.size = 0.1, ncol = 3) +
  ggtitle(paste0(sample_id, " QC metrics by DF class"))
ggsave(file.path(out_dir, paste0(sample_id, "_QC_by_DF_class.png")),
       p_vln, width = 10, height = 4, dpi = 150)

# -------------------------
# 8) Filter singlets + re-run preprocessing (recommended)
# -------------------------
singlet_cells <- rownames(seu@meta.data)[seu@meta.data[[df_class_col]] == "Singlet"]
seu_singlet <- subset(seu, cells = singlet_cells)

# Recompute SCT/PCA/UMAP after removing doublets
seu_singlet <- run_preprocess_sct(seu_singlet, dims_use = dims_use, resolution = 0.8)

p_umap_post <- DimPlot(seu_singlet, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle(paste0(sample_id, " (Singlets only, post-DF reclustering)"))
ggsave(file.path(out_dir, paste0(sample_id, "_UMAP_singlets_post_recluster.png")),
       p_umap_post, width = 7, height = 5, dpi = 150)

# -------------------------
# 9) Save outputs
# -------------------------
saveRDS(seu,         file.path(out_dir, paste0(sample_id, "_SoupX_DF_annotated_seurat.rds")))
saveRDS(seu_singlet, file.path(out_dir, paste0(sample_id, "_SoupX_DF_singlets_seurat.rds")))

# Save a barcode-level annotation table
anno <- data.frame(
  barcode = colnames(seu),
  DF_class = seu@meta.data[[df_class_col]],
  stringsAsFactors = FALSE
)
write.csv(anno, file.path(out_dir, paste0(sample_id, "_DF_barcode_annotations.csv")),
          row.names = FALSE)

message("DONE: DoubletFinder outputs saved to: ", out_dir)
