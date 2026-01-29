# Load libraries
library(Seurat)
library(SoupX)
library(tidyverse)

# Setting 
base_dir <- "/mnt/18T/chibao/pbmc/data"
out_dir <- "/mnt/18T/chibao/pbmc/data/SoupX"
# Create dir 
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Load files 
filt.mat <- Read10X_h5(file.path(base_dir, "pbmc10k_filt.h5"), use.names = TRUE)
raw.mat <- Read10X_h5(file.path(base_dir, "pbmc10k_raw.h5"), use.names = TRUE)

# Check loaded data
filt.mat
raw.mat 

# Create Seurat object from filtered matrix
seurat_obj <- CreateSeuratObject(counts = filt.mat)
seurat_obj
seurat_obj@meta.data |> head(2)

# Create SoupX object 
soupx_obj <- SoupChannel(tod = raw.mat, toc = filt.mat)
soupx_obj

# Pre-Process Seurat obj 
seurat_obj <- SCTransform(seurat_obj, verbose = FALSE)
seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, verbose = FALSE)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, verbose = FALSE)
seurat_obj <- FindClusters(seurat_obj, resolution = c(0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 1, 1.2), verbose = FALSE)

seurat_obj
seurat_obj@meta.data |> head(2)

# Extract 
meta <- seurat_obj@meta.data
meta
umap <- seurat_obj$umap@cell.embeddings
umap

# Attach meta and umap to SoupX object 
soupx_obj <- setClusters(soupx_obj, setNames(meta$seurat_clusters, rownames(meta)))
soupx_obj <- setDR(soupx_obj, umap)
soupx_obj

# Calculate ambient 
soupx_obj <- autoEstCont(soupx_obj)
soupx_obj

head(soupx_obj$soupProfile[order(soupx_obj$soupProfile$est, decreasing = T), ], n = 20)

# Adjust counts 
adj_obj <- adjustCounts(soupx_obj, roundToInt = TRUE)
adj_obj

# Write adjusted counts to 10x format
DropletUtils:::write10xCounts("soupX_pbmc10k_filt", adj_obj)
