# Load libraries
library(Seurat)
library(tidyverse)
library(glmGamPoi)
library(presto)
library(future)
library(harmony)
options(Seurat.object.assay.version = "v5")


in_path <- "/mnt/18T/chibao/pbmc/data_3/DoubletFinder_outputs/MERGED/merged_singlets_seurat.rds" # Change path if needed
out_dir <- "/mnt/18T/chibao/pbmc/data_2/Integration_outputs/PBMC_data"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Load merged obj 
merged_obj <- readRDS(in_path)
merged_obj
merged_obj@meta.data |> head()

# NEEDED if have split layers 
# merged_obj <- JoinLayers(merged_obj, assay = 'RNA)
# merged_obj

# Split the RNA assay by sample_id 
merged_obj[["RNA"]] <- split(merged_obj[["RNA"]], f = merged_obj$sample_id)
merged_obj

# Perform Cell Cycle Scoring 
DefaultAssay(merged_obj) <- "RNA"
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
merged_obj <- NormalizeData(merged_obj, verbose = FALSE)
merged_obj <- CellCycleScoring(merged_obj, s.features = s.genes, g2m.features = g2m.genes)
merged_obj
merged_obj@meta.data |> head(2)

# Perform SCT Normalize
merged_obj <- SCTransform(merged_obj, method = "glmGamPoi", vars.to.regress = c("S.Score", "GS2M.Score"), verbose = FALSE)
merged_obj 

# Save for backup 
saveRDS(merged_obj, file.path(out_dir, "merged_SCT_backup.rds"))

# Run PCA and UMAP
merged_obj <- RunPCA(merged_obj, assay = "SCT", verbose = FALSE)
merged_obj <- RunUMAP(merged_obj, dims = 1:30, verbose = FALSE)
merged_obj

# Intergate with Harmony
harmony_obj <- merged_obj

harmony_obj <- IntegrateLayers(
    obj = harmony_obj,
    method = HarmonyIntegration,
    normalization.method = "SCT",
    orig.reduction = "pca",
    new.reduction = "harmony",
    dims = 1:30,
    verbose = TRUE
)

harmony_obj

# Run PCA, UMAP, and Clustering on harmony_obj
harmony_obj <- RunPCA(harmony_obj, assay = "SCT", npcs = 30, verbose = FALSE)
harmony_obj <- RunUMAP(harmony_obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
harmony_obj <- FindNeighbors(harmony_obj, reduction = "harmony", dims = 1:30)
harmony_obj <- FindClusters(harmony_obj, resolution = c(0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1,
                                                        0.2, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4), verbose = FALSE)

harmony_obj
harmony_obj@meta.data |> head(2)

# Visualize UMAP 
p <- DimPlot(harmony_obj, reduction = "umap.harmony", groupd.by = 'SCT_snn_res.0.4', label = TRUE)
ggsave(filename = file.path(out_dir, "harmony_umap_0.04.png"), plot = p)