library(tidyverse)
library(DropletUtils)
library(Seurat)
library(Matrix)
library(scales)
library(rjson)
library(R2HTML)
library(DT)
library(data.table)
library(ComplexHeatmap)
library(circlize)

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
dim(total_filt)

#Búsqueda genes a chequear----
Presentan_5HT1A <- total_filt[
  total_filt[, "5-HT1A"] > 0, ] #Hay 34

Presentan_5HT1B <- total_filt[
  total_filt[, "5-HT1B"] > 0, ] #Hay 201

Presentan_5HT2A <- total_filt[
  total_filt[, "5-HT2A"] > 0, ] #Hay 48

Presentan_5HT7 <- total_filt[
  total_filt[, "5-HT7"] > 0, ] #Hay 26

Presentan_Tg <- total_filt[
  total_filt[, "Tg"] > 0, ] #Hay 45

Presentan_Ddc <- total_filt[
  total_filt[, "Ddc"] > 0, ] #Hay 2133



#Caracterización de las poblaciones----
total_filt_T <- t(total_filt)
counts_mat <- apply(total_filt_T, 2, function(x) as.numeric(trimws(x))) #Tarda mucho
rownames(counts_mat) <- rownames(total_filt_T)

total_seurat <- CreateSeuratObject(counts = counts_mat)
total_seurat <- NormalizeData(total_seurat, normalization.method = "LogNormalize", scale.factor = 1e6)
total_seurat <- FindVariableFeatures(total_seurat, selection.method = "vst", nfeatures = 2000)

#Filtro genes ribosomales para que no esten en las features más variables
genes_ribosomales <- grep("^Rp[SL]", rownames(total_seurat), value = TRUE)
variable_features <- VariableFeatures(total_seurat)
filtered_features <- setdiff(variable_features, genes_ribosomales)

total_seurat <- ScaleData(total_seurat, features = variable_features)
total_seurat <- RunPCA(total_seurat, features = variable_features)
ElbowPlot(total_seurat)

pcs_to_use <- 1:13
total_seurat <- FindNeighbors(total_seurat, dims = pcs_to_use, k.param = 50)  # 50 vecinos
total_seurat <- FindClusters(total_seurat, resolution = 0.4)  # resolución que da 8 clusters en paper es 0.19

#Nombres clusters
Idents(total_seurat)
total_seurat <- RenameIdents(total_seurat, 
                             '0' = '0',
                             '1' = 'MZ2',
                             '2' = 'PL',
                             '3' = '3',
                             '4' = 'CC',
                             '5' = '5',
                             '6' = 'MZ1',
                             '7' = 'PSC')
total_seurat@active.ident <- factor(total_seurat@active.ident, levels = c("PSC", "MZ1", "MZ2", "0", "3", "5", "PL", "CC"))

total_seurat <- RunTSNE(total_seurat, dims = pcs_to_use, perplexity = 50, reduction = "pca", nfeatures = 2000)
total_seurat <- RunUMAP(total_seurat, reduction = "pca", dims = pcs_to_use)
DimPlot(total_seurat, reduction = "tsne", label = TRUE)
DimPlot(total_seurat, reduction = "umap", label = TRUE)

ggsave("Clusters_tsne_TodasLasReplicas_nombradas.png", plot = DimPlot(total_seurat, reduction = "tsne", label = TRUE), width = 8, height = 6, dpi = 300)

marc_PSC <- c("Antp", "kn", "dlp") #Cluster 7
marc_MZ <- c("shg", "EGFP", "Tep4", "dome") #Cluster 1 y 6 (MZ1 y MZ2?)
marc_IZ <- c("CG30090", "lectin-24A", "MFS3") #Cluster 3
marc_CZ <- c("DsRed", "vkg", "Col4a1", "Pxn", "Hml")# Mezclado entre 0,2,3,5 no muy claro
marc_PLs <- c("NimC1") #Principalmente en cluster 2
marc_CCs <- c("PPO2", "PPO1", "lz", "peb")  #Cluster 4

#Búsqueda de genes en el cluster
FeaturePlot(total_seurat, 
            reduction = "tsne", 
            features = "Ddc",
            pt.size = 0.4, 
            order = TRUE,
            #split.by = "orig.ident",
            min.cutoff = 'q10',
            label = FALSE)

#Genes a chequear: 5-HT1A; 5-HT1B; 5-HT2A; 5-HT2B (no encontrado); 5-HT7; Thrn (no encontrado)
# Hn (no encontrado); SerT (no encontrado); Tg; Ddc

hallmark_genes <- c(marc_PSC, marc_MZ, marc_IZ, marc_CZ, marc_PLs, marc_CCs)
heatmap_marcadores_clusters <- DoHeatmap(total_seurat, features = hallmark_genes, group.by = "ident") +
  scale_fill_gradientn(colors = c("navy", "white", "firebrick3"))

ggsave("heatmap_marcadores_clusters_TodasLasReplicas.png", plot = heatmap_marcadores_clusters, width = 10, height = 8)

#Prueba Heatmap más prolijo
genes_disponibles <- rownames(GetAssayData(total_seurat, slot = "scale.data"))
hallmark_genes_filtrados <- intersect(hallmark_genes, genes_disponibles)
mat <- GetAssayData(total_seurat, slot = "scale.data")[hallmark_genes_filtrados, ]

cell_order <- order(Idents(total_seurat))
mat <- mat[, cell_order]
clusters <- as.character(Idents(total_seurat)[cell_order])

cluster_levels <- c("PSC", "MZ1", "MZ2", "0", "3", "5", "PL", "CC")  # O el orden que quieras
clusters <- factor(clusters, levels = cluster_levels)
clusters <- droplevels(clusters)
cluster_colors <- setNames(
  RColorBrewer::brewer.pal(length(cluster_levels), "Set2"),
  cluster_levels
)

top_anno <- HeatmapAnnotation(
  Cluster = clusters,
  col = list(Cluster = cluster_colors),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

heatmap_nombrado <- Heatmap(
  mat,
  name = "Z-score",
  top_annotation = top_anno,
  col = colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = TRUE,
  row_names_side = "left",
  column_split = clusters,  # separa las columnas por cluster y pone título encima
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  column_title_rot = 90,
  column_names_rot = 90,
  row_names_gp = gpar(fontsize = 8),
  use_raster = FALSE
)
png("heatmap_marcadores_clusters_TodasLasReplicas_nombradas.png", width = 10, height = 8, units = "in", res = 300)
draw(heatmap_nombrado)
dev.off()

#Subclusterización células cristal----
cc_cells <- subset(total_seurat, idents = "CC")

cc_cells <- NormalizeData(cc_cells)
cc_cells <- FindVariableFeatures(cc_cells, selection.method = "vst", nfeatures = 2000)
cc_cells <- ScaleData(cc_cells)
cc_cells <- RunPCA(cc_cells)

cc_cells <- FindNeighbors(cc_cells, dims = 1:10)
cc_cells <- FindClusters(cc_cells, resolution = 0.17)
cc_cells <- RunTSNE(cc_cells, dims = 1:10)

DimPlot(cc_cells, reduction = "tsne", label = TRUE) + ggtitle("Subcluster de células cristal (CC)")

CC_hallmark_genes <- c("PPO1", "PPO2", "lz", "peb", "Hml")

new_idents <- Idents(cc_cells)
levels(new_idents)
levels(new_idents) <- c("iCC", "mCC")
Idents(cc_cells) <- new_idents

genes_disponibles_cc <- rownames(GetAssayData(cc_cells, slot = "scale.data"))
genes_cc_filtrados <- intersect(CC_hallmark_genes, genes_disponibles_cc)

mat_cc <- GetAssayData(cc_cells, slot = "scale.data")[genes_cc_filtrados, ]
mat_cc <- mat_cc[, order(Idents(cc_cells))]

clusters_cc <- as.character(Idents(cc_cells)[colnames(mat_cc)])
clusters_cc <- factor(clusters_cc, levels = c("iCC", "mCC"))
cluster_colors_cc <- c(iCC = "#1f78b4", mCC = "#33a02c")

top_anno_cc <- HeatmapAnnotation(
  Cluster = clusters_cc,
  col = list(Cluster = cluster_colors_cc),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

heatmap_cc <- Heatmap(
  mat_cc,
  name = "Z-score",
  top_annotation = top_anno_cc,
  col = colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = TRUE,
  row_names_side = "left",
  column_split = clusters_cc,
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  row_names_gp = gpar(fontsize = 10),
  use_raster = FALSE
)

png("heatmap_CC.png", width = 10, height = 8, units = "in", res = 300)
draw(heatmap_cc)
dev.off()

#Observar grupo de células en particular----
genes_pobPrevia <- FetchData(total_seurat, vars = c("dome", "Tep4", "Hml"))
cells_interes <- rownames(genes_pobPrevia)[
  genes_pobPrevia$dome > 0 &
  genes_pobPrevia$Tep4 > 0 &
  genes_pobPrevia$Hml == 0]

total_seurat$marker_combo <- ifelse(colnames(total_seurat) %in% cells_interes, "Previa", "Otros")
p <- DimPlot(
  total_seurat,
  group.by = "marker_combo",
  reduction = "tsne",
  pt.size = 0.4,
  cols = c("lightgray", "red"))

p + ggtitle("Ubicación de células dome+, Tep4+, Hml-")


#Genes que forman parte de “Regulation of intracellular signal transduction”(GO:1902531)----
genes_go1902531 <- c("sgg", "Mo25", "sra", "cact", "Hipk", "msn", "CG42674", "Myd88", "Cdc37", 
                     "Crag", "CG12290", "Atg1", "RhoGEF3", "Smurf", "M1BP", "Csk", "wcy", 
                     "Cka", "Patronin", "raw", "Nup44A", "Clbn", "Dhx15", "Dad", "pll", "dos", 
                     "step", "wrd", "Cbl", "aru", "CG18659", "alph", "Sac1", "Src42A", "RhoGEF2", 
                     "Pak", "Sur-8", "Stlk", "Slmap", "armi", "Pde8", "PDZ-GEF", "key", "Lst", 
                     "LRR", "PEK", "Cul4", "Pp2C1", "AMPKalpha", "CycD", "sigmar", "RhoGAP1A", 
                     "jub", "Kcmf1", "spri", "RagA-B", "trbl", "Pdk1", "Mnn1", "PRAS40", "Egfr", 
                     "Duba", "conu", "RtGEF", "garz", "rictor", "Mer", "Tnks", "mbt", "RhoGAP5A", 
                     "Socs16D", "hppy", "Src64B", "ex", "Ced-12", "CG4853", "cnk", "tefu", "sl", 
                     "CG5521", "Rgl", "Art4", "CG15547", "cno", "Rlip", "chico", "uri", "hpo", 
                     "puc", "trr", "Tsc1", "Lpt", "put")

genes_disponibles <- rownames(GetAssayData(total_seurat, slot = "scale.data"))
genes_go1902531_filtrados <- intersect(genes_go1902531, genes_disponibles) #Genes de go1902531 presentes entre los 2000 más variados

mat_go1902531 <- GetAssayData(total_seurat, slot = "scale.data")[genes_go1902531_filtrados, ]
mat_go1902531 <- mat_go1902531[, order(Idents(total_seurat))]

clusters_go1902531 <- as.character(Idents(total_seurat)[colnames(mat_go1902531)])
clusters_go1902531 <- factor(clusters_go1902531, levels = cluster_levels)
clusters_go1902531 <- droplevels(clusters_go1902531)

top_anno_go1902531 <- HeatmapAnnotation(
  Cluster = clusters_go1902531,
  col = list(Cluster = cluster_colors),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

heatmap_go1902531 <- Heatmap(
  mat_go1902531,
  name = "Z-score",
  top_annotation = top_anno_go1902531,
  col = colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3")),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = TRUE,
  row_names_side = "left",
  column_split = clusters_go1902531,
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  column_title_rot = 90,
  row_names_gp = gpar(fontsize = 6),
  use_raster = FALSE
)

png("heatmap_genes_GO1902531.png", width = 10, height = 8, units = "in", res = 300)
draw(heatmap_go1902531)
dev.off()
