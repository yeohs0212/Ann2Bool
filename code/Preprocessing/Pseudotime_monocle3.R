library(Seurat)
library(dplyr)
library(ggplot2)
library(anndata)
library(tidyverse)
data = read_h5ad("/data/yeohs0212/MM/GSE207938/03_ADM_for_trajectory.h5ad")
library(monocle3)
seu <- CreateSeuratObject(counts = t(as.matrix(data$X)), meta.data = data$obs)
saveRDS(data,"/data/yeohs0212/MM/02_PDAC_subset_for_trajectory.rds")
custom_coords <- read.csv("/data/yeohs0212/MM/GSE207938/umap_coordinates.csv", row.names = 1)
exp_mat <- seu@assays$RNA$counts
cell_meta <- seu@meta.data
gene_meta <- data.frame(
  gene_short_name = rownames(exp_mat),
  row.names = rownames(exp_mat)
)


cds <- new_cell_data_set(exp_mat,
                         cell_metadata = cell_meta,
                         gene_metadata = gene_meta)

cds <- preprocess_cds(cds, num_dim = 50)

reducedDims(cds)$UMAP <- as.matrix(custom_coords[colnames(cds), ])
cds <- cluster_cells(cds, reduction_method = "UMAP")

# 5. 궤적 그래프 학습 (Learn Graph)
# 여기서 기존 좌표를 유지한 채 선을 긋습니다.
cds <- learn_graph(cds, use_partition = FALSE)
plot_cells(cds,color_cells_by = "cell_type")
cds <- order_cells(cds, root_cells = "235737665616812_DACD511_Kate_plus")
plot_cells(cds, 
           color_cells_by = "pseudotime", 
           label_cell_groups = FALSE,
           label_leaves = TRUE, 
           label_branch_points = TRUE,
           graph_label_size = 3)


pt_values <- pseudotime(cds)
pt_df <- as.data.frame(pt_values)
colnames(pt_df) <- "PT"

pt_df <- pt_df %>% 
  rownames_to_column(var = "cell_id")
write.csv(pt_df, "/data/yeohs0212/MM/GSE207938/pseudotime_results.csv",row.names = FALSE)
