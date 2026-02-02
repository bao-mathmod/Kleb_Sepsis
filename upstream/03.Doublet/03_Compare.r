suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(Seurat)
  library(gridExtra)
  library(ggVennDiagram)
})

# ==============================================================================
# USER SETTINGS
# ==============================================================================
base_dir <- "/mnt/18T/chibao/pbmc/data"
sample_id <- "PBMC_data"  # Change this to loop if needed later

# Paths based on your structure
df_dir  <- file.path(base_dir, "DoubletFinder_outputs", sample_id)
scr_dir <- file.path(base_dir, "Scrublet_outputs", sample_id)
out_dir <- file.path(base_dir, "Doublet_Comparison_Benchmark", sample_id)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================
message("Loading annotations for: ", sample_id)

# Load DoubletFinder CSV
df_csv <- file.path(df_dir, paste0(sample_id, "_DF_barcode_annotations.csv"))
if (!file.exists(df_csv)) stop("Missing DF CSV: ", df_csv)
df_data <- read.csv(df_csv, stringsAsFactors = FALSE)
colnames(df_data)[colnames(df_data) == "DF_class"] <- "DF_pred" # Standardize

# Load Scrublet CSV
scr_csv <- file.path(scr_dir, paste0(sample_id, "_Scrublet_barcodes.csv"))
if (!file.exists(scr_csv)) stop("Missing Scrublet CSV: ", scr_csv)
scr_data <- read.csv(scr_csv, stringsAsFactors = FALSE)

# Merge by Barcode
# Inner join ensures we only compare cells present in both (should be identical)
merged <- inner_join(df_data, scr_data, by = "barcode")

# Standardize labels (Singlet/Doublet)
merged$DF_pred  <- ifelse(merged$DF_pred == "Doublet", "Doublet", "Singlet")
merged$Scrublet_pred <- ifelse(merged$Scrublet_pred == "Doublet", "Doublet", "Singlet")

# Create Consensus Column
merged$Consensus <- case_when(
  merged$DF_pred == "Doublet" & merged$Scrublet_pred == "Doublet" ~ "Consensus_Doublet",
  merged$DF_pred == "Doublet" & merged$Scrublet_pred == "Singlet" ~ "DF_Only",
  merged$DF_pred == "Singlet" & merged$Scrublet_pred == "Doublet" ~ "Scrublet_Only",
  TRUE ~ "Singlet"
)

# ==============================================================================
# 2. STATISTICAL BENCHMARKING
# ==============================================================================
# Confusion Matrix
conf_mat <- table(DoubletFinder = merged$DF_pred, Scrublet = merged$Scrublet_pred)
print(conf_mat)
write.csv(as.data.frame(conf_mat), file.path(out_dir, "Confusion_Matrix.csv"))

# Calculate Metrics
n_total <- nrow(merged)
n_df    <- sum(merged$DF_pred == "Doublet")
n_scr   <- sum(merged$Scrublet_pred == "Doublet")
n_both  <- sum(merged$Consensus == "Consensus_Doublet")
n_union <- sum(merged$Consensus != "Singlet")

stats_df <- data.frame(
  Metric = c("Total Cells", "DF Called Doublets", "Scrublet Called Doublets", 
             "Intersection (Both)", "Union (Either)", "Jaccard Index (Overlap)"),
  Value = c(n_total, n_df, n_scr, n_both, n_union, 
            round(n_both / n_union, 3))
)
write.csv(stats_df, file.path(out_dir, "Benchmark_Stats.csv"), row.names = FALSE)

message("\n--- BENCHMARK SUMMARY ---")
print(stats_df)

# ==============================================================================
# 3. VISUALIZATION (ggVennDiagram)
# ==============================================================================
# Prepare list for plotting
venn_list <- list(
  DoubletFinder = merged$barcode[merged$DF_pred == "Doublet"],
  Scrublet = merged$barcode[merged$Scrublet_pred == "Doublet"]
)

# Create Plot
p_venn <- ggVennDiagram(venn_list) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  ggtitle(paste0("Doublet Overlap: ", sample_id))

# Save
ggsave(file.path(out_dir, "Venn_Overlap_gg.png"), p_venn, width = 6, height = 6)

# ==============================================================================
# 4. UMAP CONSENSUS (Crucial Step)
# ==============================================================================
# Load the Seurat object (We can use either annotated one, they have the same UMAP)
# We'll use the Scrublet one as it's cleaner to load fresh
seu_path <- file.path(scr_dir, paste0(sample_id, "_SoupX_Scrublet_annotated_seurat.rds"))
seu <- readRDS(seu_path)

# Add the Consensus metadata
rownames(merged) <- merged$barcode
seu <- AddMetaData(seu, metadata = merged[colnames(seu), "Consensus", drop=FALSE], col.name = "Consensus_Class")

# 1. Consensus Plot
# Colors: Singlet=Grey, Both=Red (High Conf), DF_Only=Blue, Scrublet_Only=Green
cols <- c("Singlet" = "lightgrey", 
          "Consensus_Doublet" = "#D55E00",  # Dark Orange/Red
          "DF_Only" = "#0072B2",            # Blue
          "Scrublet_Only" = "#009E73")      # Green

p_con <- DimPlot(seu, group.by = "Consensus_Class", cols = cols, order = c("Consensus_Doublet", "DF_Only", "Scrublet_Only")) +
  ggtitle(paste0(sample_id, " Consensus Doublet Calls"))
ggsave(file.path(out_dir, "UMAP_Consensus_Comparison.png"), p_con, width = 8, height = 6)

# 2. Score Correlation Plot
# Do high DF scores correlate with high Scrublet scores?
# Note: DF doesn't output a raw score in the CSV by default, only classification.
# However, Scrublet does. We can visualize Scrublet Score vs Consensus Class.
p_box <- VlnPlot(seu, features = "Scrublet_score", group.by = "Consensus_Class", pt.size = 0.1) +
  ggtitle("Scrublet Scores split by Consensus Class")
ggsave(file.path(out_dir, "Violin_ScrubletScore_by_Consensus.png"), p_box, width = 8, height = 6)

message("✅ Comparison Complete. Check folder: ", out_dir)