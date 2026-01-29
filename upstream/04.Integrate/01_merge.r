suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

# =========================
# INPUTS
# =========================
df_base <- "/mnt/18T/chibao/pbmc/data_3/DoubletFinder_outputs" # Change path if needed 

out_dir <- file.path(df_base, "MERGED")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =========================
# 1) Find singlet RDS files
# =========================
rds_files <- list.files(
  df_base,
  pattern = "_SoupX_DF_singlets_seurat\\.rds$",
  full.names = TRUE,
  recursive = TRUE
)

stopifnot(length(rds_files) > 1)

# sample_id inferred from parent folder
sample_ids <- basename(dirname(rds_files))

message("Found ", length(rds_files), " samples:")
print(sample_ids)

# =========================
# 2) Read objects + enforce unique cell names
# =========================
objs <- vector("list", length(rds_files))
names(objs) <- sample_ids

for (i in seq_along(rds_files)) {
  sid <- sample_ids[i]
  obj <- readRDS(rds_files[i])

  # Ensure we can track sample ID later
  obj$sample_id <- sid

  # NEEDED if have split layers:
  # obj <- JoinLayers(obj, assay = "RNA")

  # Make cell names unique across samples (critical)
  # Example: AAAC...-1 becomes PBMC_data_AAAC...-1
  obj <- RenameCells(obj, add.cell.id = sid)

  objs[[sid]] <- obj
}

# Quick numeric check (example output: number of cells per sample)
cell_counts <- sapply(objs, ncol)
print(cell_counts)

# =========================
# 3) Merge
# =========================
merged <- merge(
  x = objs[[1]],
  y = objs[-1],
  project = "PBMC_sepsis_merged"
)

# Another quick check: cell totals
message("Merged object cells: ", ncol(merged),
        " | sum of per-sample cells: ", sum(cell_counts))

# =========================
# 4) Save merged object
# =========================
saveRDS(merged, file.path(out_dir, "PBMC_sepsis_merged_singlets.rds"))

# =========================
# 5) Save summary table + plot (numbers + graph)
# =========================
summary_df <- data.frame(
  sample_id = names(cell_counts),
  n_cells = as.integer(cell_counts)
) %>% arrange(desc(n_cells))

write.csv(summary_df, file.path(out_dir, "merged_cell_counts_by_sample.csv"),
          row.names = FALSE)

p <- ggplot(summary_df, aes(x = reorder(sample_id, n_cells), y = n_cells)) +
  geom_col() +
  coord_flip() +
  xlab("Sample") + ylab("Number of singlet cells") +
  ggtitle("Cells per sample after SoupX + DoubletFinder (singlets)")

ggsave(file.path(out_dir, "merged_cell_counts_by_sample.png"),
       p, width = 8, height = 6, dpi = 150)

message("DONE. Merged object + summaries saved in: ", out_dir)
