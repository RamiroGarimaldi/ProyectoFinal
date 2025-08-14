library(data.table)
library(Seurat)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(RColorBrewer)

#Preparación de la matriz para el análisis
rep1 <- fread("DatosCrudos/Rep1.txt", sep = "\t")
rep2 <- fread("DatosCrudos/Rep2.txt", sep = "\t")
rep3 <- fread("DatosCrudos/Rep3.txt", sep = "\t")

rep1[, Name := paste0("rep1_", Name)]
rep1_df <- as.data.frame(rep1)
rownames(rep1_df) <- rep1_df$Name
rep1_df$Name <- NULL

rep2[, Name := paste0("rep2_", Name)]
rep2_df <- as.data.frame(rep2)
rownames(rep2_df) <- rep2_df$Name
rep2_df$Name <- NULL

rep3[, Name := paste0("rep3_", Name)]
rep3_df <- as.data.frame(rep3)
rownames(rep3_df) <- rep3_df$Name
rep3_df$Name <- NULL

total <- rbind(rep1_df, rep2_df, rep3_df)

genes_a_buscar <- c("lncRNA:roX1", "lncRNA:roX2", "lncRNA:CR40469")
genes_presentes <- genes_a_buscar[genes_a_buscar %in% colnames(total)]
genes_rRNA <- grep("rRNA", colnames(total), value = TRUE)
genes_a_filtrar <- c(genes_presentes, genes_rRNA)

total_df <- as.data.frame(total)
total_filt <- total_df[, !(colnames(total_df) %in% genes_a_filtrar), drop = FALSE]

total_filt_T <- t(total_filt)
counts_mat <- apply(total_filt_T, 2, function(x) as.numeric(trimws(x))) #Tarda mucho
rownames(counts_mat) <- rownames(total_filt_T)

#Una vez preparada la matriz, creo el objeto de Seurat
total_seurat <- CreateSeuratObject(counts = counts_mat)
total_seurat <- NormalizeData(total_seurat, normalization.method = "LogNormalize", scale.factor = 1e6)
total_seurat <- FindVariableFeatures(total_seurat, selection.method = "vst", nfeatures = 2000)

genes_ribosomales <- grep("^Rp[SL]", rownames(total_seurat), value = TRUE)
variable_features <- setdiff(VariableFeatures(total_seurat), genes_ribosomales)

total_seurat <- ScaleData(total_seurat, features = variable_features)
total_seurat <- RunPCA(total_seurat, features = variable_features)
ElbowPlot(total_seurat)

set.seed(1234) #Para asegurar la reproducibilidad, y que las semillas no sean siempre al azar

pcs_to_use <- 1:13
total_seurat <- FindNeighbors(total_seurat, dims = pcs_to_use, k.param = 50)  # 50 vecinos
total_seurat <- FindClusters(total_seurat, resolution = 0.4) #En mi codigo use 0.4, ver si da bien o volver a mi valor

#Grafico los resultados obtenidos para observar los distintos clusters
total_seurat <- RunTSNE(total_seurat, dims = pcs_to_use, perplexity = 50, reduction = "pca")
DimPlot(total_seurat, reduction = "tsne", label = TRUE) + coord_flip()

Idents(total_seurat)

#Renombro los clusters según los genes característicos
total_seurat <- RenameIdents(total_seurat, 
                             '0' = '0',
                             '1' = 'MZ2',
                             '2' = 'PL',
                             '3' = 'IZ',
                             '4' = 'CC',
                             '5' = '5',
                             '6' = 'MZ1',
                             '7' = 'PSC')
total_seurat@active.ident <- factor(total_seurat@active.ident, levels = c("PSC", "MZ1", "MZ2", "IZ", "0", "5", "PL", "CC"))

#Elaboración del heatmap para verificar que los genes marcadores sean los correspondientes para cada cluster
hallmark_genes <- c("Antp", "kn", "dlp", "Dl", "shg", "EGFP", "Tep4", "dome", "CG30090",    
                    "lectin-24A", "MFS3", "DsRed", "NimC1",  "vkg", "Col4a1", "Pxn",
                    "Hml", "PPO2", "PPO1", "lz", "peb") #Orden del paper

genes_disponibles <- rownames(GetAssayData(total_seurat, slot = "scale.data"))
hallmark_genes_filtrados <- intersect(hallmark_genes, genes_disponibles)
mat <- GetAssayData(total_seurat, slot = "scale.data")[hallmark_genes_filtrados, ]

cell_order <- order(Idents(total_seurat))
mat <- mat[, cell_order]
clusters <- as.character(Idents(total_seurat)[cell_order])

cluster_levels <- c("PSC", "MZ1", "MZ2", "IZ", "0", "5", "PL", "CC")
clusters <- factor(clusters, levels = cluster_levels)
clusters <- droplevels(clusters)
cluster_colors <- setNames(RColorBrewer::brewer.pal(length(cluster_levels), "Set2"), cluster_levels)

top_anno <- HeatmapAnnotation(
  Cluster = clusters,
  col = list(Cluster = cluster_colors),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

heatmap_cluster_propio <- Heatmap(
  mat,
  name = "Z-score",
  top_annotation = top_anno,
  col = colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = TRUE,
  row_names_side = "left",
  column_split = clusters,
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  column_title_rot = 90,
  column_names_rot = 90,
  row_names_gp = gpar(fontsize = 8),
  use_raster = FALSE
)

#Elaboración del heatmap a partir de la expresión génica promedio por cluster
expr_data_rna <- GetAssayData(total_seurat, slot = "data")[hallmark_genes_filtrados, ]
ident_vector <- Idents(total_seurat)

avg_expr_rna <- sapply(levels(ident_vector), function(cluster) {
  cells_in_cluster <- WhichCells(total_seurat, idents = cluster)
  rowMeans(expr_data_rna[, cells_in_cluster, drop = FALSE])
})

avg_expr_rna_zscore <- t(scale(t(avg_expr_rna)))

cluster_levels <- c("PSC", "MZ1", "MZ2", "IZ", "0", "5", "PL", "CC")
colnames(avg_expr_rna_zscore) <- factor(colnames(avg_expr_rna_zscore), levels = cluster_levels)
avg_expr_rna_zscore <- avg_expr_rna_zscore[, cluster_levels]
cluster_colors <- setNames(RColorBrewer::brewer.pal(length(cluster_levels), "Set2"), cluster_levels)

top_anno_avg <- HeatmapAnnotation(
  Cluster = colnames(avg_expr_rna_zscore),
  col = list(Cluster = cluster_colors),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

heatmap_rna_prom <- Heatmap(
  avg_expr_rna_zscore,
  name = "Z-score",
  top_annotation = top_anno_avg,
  col = colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = TRUE,
  row_names_side = "left",
  column_split = factor(cluster_levels, levels = cluster_levels),
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_title_rot = 90,
  row_names_gp = gpar(fontsize = 8),
  use_raster = FALSE
)


