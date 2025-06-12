#Paquetes utilizados----
library(data.table)
library(edgeR)
library(limma)
library(clusterProfiler)
library(org.Dm.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(stringr)

#Lectura de archivos----
list.files("DatosCrudos")
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

genes_sexo <- c("lncRNA:roX1", "lncRNA:roX2", "lncRNA:CR40469")
genes_rRNA <- grep("rRNA", colnames(total), value = TRUE)
genes_a_filtrar <- c(genes_sexo, genes_rRNA)
total_filt <-  total[, !colnames(total) %in% genes_a_filtrar]

Poblacion_interes <- total_filt[
  total_filt[, "dome"] > 0 & 
    total_filt[, "Tep4"] == 0 & 
    total_filt[, "Hml"] == 0, 
] #Hay 1

Poblacion_previa <- total_filt[
  total_filt[, "dome"] > 0 & 
    total_filt[, "Tep4"] > 0 & 
    total_filt[, "Hml"] == 0, 
] #Hay 171

Poblacion_post <- total_filt[
  total_filt[, "dome"] > 0 & 
    total_filt[, "Tep4"]== 0 & 
    total_filt[, "Hml"] > 0, 
] #Hay 90

#Análisis de las poblaciones----

#Análisis población de interés----
matriz_interes <- total_filt[rownames(total_filt) %in% rownames(Poblacion_interes), ]
matriz_interes_t <- t(matriz_interes)
expresion_promedio <- rowMeans(matriz_interes_t)
top_genes_interes_sin_contraste <- names(sort(expresion_promedio, decreasing = TRUE))[1:200]
top_genes_entrez <- bitr(top_genes_interes_sin_contraste, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Dm.eg.db)$ENTREZID

ego_interes_sin_contraste <- enrichGO(
  gene = top_genes_entrez,
  OrgDb = org.Dm.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

dotplot(ego_interes_sin_contraste, showCategory = 15) +
        ggtitle("GO - Población de Interés") +
        theme(axis.text.y = element_text(size = 8))

ggsave("GO_PoblacionInteres.png", 
        plot = dotplot(ego_interes_sin_contraste, showCategory = 15) + 
        ggtitle("GO - Población de Interés") +
        theme(axis.text.y = element_text(size = 8)),
        width = 10, height = 6, dpi = 300)

ego_interes_df_sin_contraste <- as.data.frame(ego_interes_sin_contraste)

tabla_interes_sin_contraste <- ego_interes_df_sin_contraste %>%
  mutate(
    N_genes_enriquecidos = str_count(geneID, "/") + 1,
    Genes_representativos = sapply(strsplit(geneID, "/"), function(x) paste(head(x, 5), collapse = ", "))
  ) %>%
  select(Description, N_genes_enriquecidos, Count, p.adjust, Genes_representativos)

write.csv(tabla_interes_sin_contraste, "Tabla_GOPoblacionInteres.csv", row.names = FALSE)
Tabla_GOPoblacionesInteres <- fread("Tabla_GOPoblacionInteres.csv")


#Análisis población previa----
matriz_previa <- total_filt[rownames(total_filt) %in% rownames(Poblacion_previa), ]
matriz_previa_t <- t(matriz_previa)
expresion_promedio_previa <- rowMeans(matriz_previa_t)
top_genes_previa_sin_contraste <- names(sort(expresion_promedio_previa, decreasing = TRUE))[1:200]
top_genes_previa_entrez <- bitr(top_genes_previa_sin_contraste, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Dm.eg.db)$ENTREZID

ego_previa_sin_contraste <- enrichGO(
  gene = top_genes_previa_entrez,
  OrgDb = org.Dm.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

dotplot(ego_previa_sin_contraste, showCategory = 15) +
  ggtitle("GO - Población previa") +
  theme(axis.text.y = element_text(size = 8))

ggsave("GO_PoblacionPrevia.png", 
       plot = dotplot(ego_previa_sin_contraste, showCategory = 15) + 
         ggtitle("GO - Población previa") +
         theme(axis.text.y = element_text(size = 8)),
       width = 10, height = 6, dpi = 300)

ego_previa_df_sin_contraste <- as.data.frame(ego_previa_sin_contraste)

tabla_previa_sin_contraste <- ego_previa_df_sin_contraste %>%
  mutate(
    N_genes_enriquecidos_previa = str_count(geneID, "/") + 1,
    Genes_representativos_previa = sapply(strsplit(geneID, "/"), function(x) paste(head(x, 5), collapse = ", "))
  ) %>%
  select(Description, N_genes_enriquecidos_previa, Count, p.adjust, Genes_representativos_previa)

write.csv(tabla_previa_sin_contraste, "Tabla_GOPoblacionPrevia.csv", row.names = FALSE)
Tabla_GOPoblacionPrevia <- fread("Tabla_GOPoblacionPrevia.csv")

#Análisis población post----
matriz_post <- total_filt[rownames(total_filt) %in% rownames(Poblacion_post), ]
matriz_post_t <- t(matriz_post)
expresion_promedio_post <- rowMeans(matriz_post_t)
top_genes_post_sin_contraste <- names(sort(expresion_promedio_post, decreasing = TRUE))[1:200]
top_genes_post_entrez <- bitr(top_genes_post_sin_contraste, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Dm.eg.db)$ENTREZID

ego_post_sin_contraste <- enrichGO(
  gene = top_genes_post_entrez,
  OrgDb = org.Dm.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

dotplot(ego_post_sin_contraste, showCategory = 15) +
  ggtitle("GO - Población post") +
  theme(axis.text.y = element_text(size = 8))

ggsave("GO_PoblacionPost.png", 
       plot = dotplot(ego_post_sin_contraste, showCategory = 15) + 
         ggtitle("GO - Población post") +
         theme(axis.text.y = element_text(size = 8)),
       width = 10, height = 6, dpi = 300)

ego_post_df_sin_contraste <- as.data.frame(ego_post_sin_contraste)

tabla_post_sin_contraste <- ego_post_df_sin_contraste %>%
  mutate(
    N_genes_enriquecidos_post = str_count(geneID, "/") + 1,
    Genes_representativos_post = sapply(strsplit(geneID, "/"), function(x) paste(head(x, 5), collapse = ", "))
  ) %>%
  select(Description, N_genes_enriquecidos_post, Count, p.adjust, Genes_representativos_post)

write.csv(tabla_post_sin_contraste, "Tabla_GOPoblacionPost.csv", row.names = FALSE)
Tabla_GOPoblacionPost <- fread("Tabla_GOPoblacionPost.csv")


#Comparación entre poblaciones----
poblaciones_a_estudiar <- data.frame(
  celda = rownames(total_filt),
  grupo = NA_character_
)

poblaciones_a_estudiar$grupo[total_filt$dome > 0 & total_filt$Tep4 == 0 & total_filt$Hml == 0] <- "interes"
poblaciones_a_estudiar$grupo[total_filt$dome > 0 & total_filt$Tep4 > 0 & total_filt$Hml == 0] <- "previa"
poblaciones_a_estudiar$grupo[total_filt$dome > 0 & total_filt$Tep4 == 0 & total_filt$Hml > 0] <- "post"

poblaciones_a_estudiar <- poblaciones_a_estudiar[!is.na(poblaciones_a_estudiar$grupo), ]
matriz_poblaciones <- total_filt[rownames(total_filt) %in% poblaciones_a_estudiar$celda, ]
poblaciones_a_estudiar <- poblaciones_a_estudiar[match(rownames(matriz_poblaciones), poblaciones_a_estudiar$celda), ]

grupo_factor <- factor(poblaciones_a_estudiar$grupo)
design <- model.matrix(~0 + grupo_factor)
colnames(design) <- levels(grupo_factor)

if (nrow(matriz_poblaciones) == nrow(poblaciones_a_estudiar)) {
  matriz_poblaciones <- t(matriz_poblaciones)
}

dge <- DGEList(counts = matriz_poblaciones)
dge <- calcNormFactors(dge)
v <- voom(dge, design)
fit <- lmFit(v, design)

# Interes vs Previa----
contraste <- makeContrasts(interes_vs_previa = interes - previa, levels = design)
fit2 <- contrasts.fit(fit, contraste)
fit2 <- eBayes(fit2)
top_genes <- topTable(fit2, coef = "interes_vs_previa", number = Inf, adjust = "fdr", sort.by = "P")

genes_top_interes <- rownames(top_genes)[top_genes$logFC > 1 & top_genes$adj.P.Val < 0.05]
genes_top_previa <- rownames(top_genes)[top_genes$logFC < -1 & top_genes$adj.P.Val < 0.05]
#Conversion a ID para flybase
genes_top_interes_entrez <- bitr(genes_top_interes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Dm.eg.db)$ENTREZID
genes_top_previa_entrez <- bitr(genes_top_previa, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Dm.eg.db)$ENTREZID

#Analisis excluyendo IDs faltantes
#Genes enriquecidos en población de interés frente a previa
ego_interes <- enrichGO(
  gene = genes_top_interes_entrez,
  OrgDb = org.Dm.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)
dotplot(ego_interes, showCategory = 15) + ggtitle("GO enriquecido - Población de Interés")

ggsave("GO_Interes_InteresvsPrevia.png", 
        plot = dotplot(ego_interes, showCategory = 15) + 
        ggtitle("GO enriquecido - Población de Interés"),
        width = 10, height = 6, dpi = 300)

ego_interes_df <- as.data.frame(ego_interes)
Tabla_GOInteres_interesvsprevia <- ego_interes_df %>%
  mutate(
    N_genes_enriquecidos = str_count(geneID, "/") + 1,  # Cuenta los genes en geneID
    Genes_representativos = sapply(strsplit(geneID, "/"), function(x) paste(x, collapse = ", "))
  ) %>%
  dplyr::select(Description, N_genes_enriquecidos, Count, p.adjust, Genes_representativos)

write.csv(Tabla_GOInteres_interesvsprevia, "Tabla_GOInteres_interesvsprevia.csv", row.names = FALSE)

#Genes enriquecidos en población previa frente a la de interés
ego_previa <- enrichGO(
  gene = genes_top_previa_entrez,
  OrgDb = org.Dm.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)
dotplot(ego_previa, showCategory = 15) + ggtitle("GO enriquecido - Población Previa")

ggsave("GO_Previa_InteresvsPrevia.png", 
        plot = dotplot(ego_previa, showCategory = 15) + 
        ggtitle("GO enriquecido - Población previa"),
        width = 10, height = 6, dpi = 300)

ego_previa_df <- as.data.frame(ego_previa)
Tabla_GOPrevia_interesvsprevia <- ego_previa_df %>%
  mutate(
    N_genes_enriquecidos = str_count(geneID, "/") + 1,  # Cuenta los genes en geneID
    Genes_representativos = sapply(strsplit(geneID, "/"), function(x) paste(x, collapse = ", "))
  ) %>%
  dplyr::select(Description, N_genes_enriquecidos, Count, p.adjust, Genes_representativos)

write.csv(Tabla_GOPrevia_interesvsprevia, "Tabla_GOPrevia_interesvsprevia.csv", row.names = FALSE)

#Agregar IDs faltantes
mapped <- bitr(genes_top_interes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Dm.eg.db)
no_mapped <- setdiff(genes_top_interes, mapped$SYMBOL)
print(no_mapped)


# Interes vs Post----
contrasteIvPost <- makeContrasts(interes_vs_post = interes - post, levels = design)
fit3 <- contrasts.fit(fit, contrasteIvPost)
fit3 <- eBayes(fit3)
top_genes_IvPost <- topTable(fit3, coef = "interes_vs_post", number = Inf, adjust = "fdr", sort.by = "P")

genes_top_interes2 <- rownames(top_genes_IvPost)[top_genes_IvPost$logFC > 1 & top_genes_IvPost$adj.P.Val < 0.05]
genes_top_post <- rownames(top_genes_IvPost)[top_genes_IvPost$logFC < -1 & top_genes_IvPost$adj.P.Val < 0.05]
#Conversion a ID para flybase
genes_top_interes2_entrez <- bitr(genes_top_interes2, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Dm.eg.db)$ENTREZID
genes_top_post_entrez <- bitr(genes_top_post, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Dm.eg.db)$ENTREZID

#Analisis excluyendo IDs faltantes
#Genes enriquecidos en población de interés frente a previa
ego_interes2 <- enrichGO(
  gene = genes_top_interes2_entrez,
  OrgDb = org.Dm.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)
dotplot(ego_interes2, showCategory = 15) + ggtitle("GO enriquecido - Población de Interés")

ggsave("GO_Interes_InteresvsPost.png", 
       plot = dotplot(ego_interes2, showCategory = 15) + 
         ggtitle("GO enriquecido - Población de Interés"),
       width = 10, height = 6, dpi = 300)

ego_interes2_df <- as.data.frame(ego_interes2)
Tabla_GOInteres_interesvspost <- ego_interes2_df %>%
  mutate(
    N_genes_enriquecidos = str_count(geneID, "/") + 1,  # Cuenta los genes en geneID
    Genes_representativos = sapply(strsplit(geneID, "/"), function(x) paste(head(x, 5), collapse = ", "))
  ) %>%
  select(Description, N_genes_enriquecidos, Count, p.adjust, Genes_representativos)

write.csv(Tabla_GOInteres_interesvspost, "Tabla_GOInteres_interesvspost.csv", row.names = FALSE)

#Genes enriquecidos en población post frente a la de interés
ego_post <- enrichGO(
  gene = genes_top_post_entrez,
  OrgDb = org.Dm.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)
dotplot(ego_post, showCategory = 15) + ggtitle("GO enriquecido - Población Post")

ggsave("GO_Post_InteresvsPost.png", 
       plot = dotplot(ego_post, showCategory = 15) + 
         ggtitle("GO enriquecido - Población Post") +
         theme(axis.text.y = element_text(size = 8)),
       width = 10, height = 6, dpi = 300)

ego_post_df <- as.data.frame(ego_post)
Tabla_GOPost_interesvspost <- ego_post_df %>%
  mutate(
    N_genes_enriquecidos = str_count(geneID, "/") + 1,  # Cuenta los genes en geneID
    Genes_representativos = sapply(strsplit(geneID, "/"), function(x) paste(head(x, 5), collapse = ", "))
  ) %>%
  select(Description, N_genes_enriquecidos, Count, p.adjust, Genes_representativos)

write.csv(Tabla_GOPost_interesvspost, "Tabla_GOPost_interesvspost.csv", row.names = FALSE)

# Previa vs Post----
contrastePrevPost <- makeContrasts(pre_vs_post = previa - post, levels = design)
fit4 <- contrasts.fit(fit, contrastePrevPost)
fit4 <- eBayes(fit4)
top_genes_PrevPost <- topTable(fit4, coef = "pre_vs_post", number = Inf, adjust = "fdr", sort.by = "P")

genes_top_pre2 <- rownames(top_genes_PrevPost)[top_genes_PrevPost$logFC > 1 & top_genes_PrevPost$adj.P.Val < 0.05]
genes_top_post2 <- rownames(top_genes_PrevPost)[top_genes_PrevPost$logFC < -1 & top_genes_PrevPost$adj.P.Val < 0.05]
#Conversion a ID para flybase
genes_top_pre2_entrez <- bitr(genes_top_pre2, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Dm.eg.db)$ENTREZID
genes_top_post2_entrez <- bitr(genes_top_post2, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Dm.eg.db)$ENTREZID

#Analisis excluyendo IDs faltantes
#Genes enriquecidos en población de interés frente a previa
ego_pre2 <- enrichGO(
  gene = genes_top_pre2_entrez,
  OrgDb = org.Dm.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)
dotplot(ego_pre2, showCategory = 15) + ggtitle("GO enriquecido - Población Previa")

ggsave("GO_Previa_PrevsPost.png", 
       plot = dotplot(ego_pre2, showCategory = 15) + 
         ggtitle("GO enriquecido - Población Previa"),
       width = 10, height = 6, dpi = 300)

ego_pre2_df <- as.data.frame(ego_pre2)
Tabla_GOPrevia_prevspost <- ego_pre2_df %>%
  mutate(
    N_genes_enriquecidos = str_count(geneID, "/") + 1,  # Cuenta los genes en geneID
    Genes_representativos = sapply(strsplit(geneID, "/"), function(x) paste(head(x, 5), collapse = ", "))
  ) %>%
  select(Description, N_genes_enriquecidos, Count, p.adjust, Genes_representativos)

write.csv(Tabla_GOPrevia_prevspost, "Tabla_GOPrevia_prevspost.csv", row.names = FALSE)

#Genes enriquecidos en población post frente a la de interés
ego_post2 <- enrichGO(
  gene = genes_top_post2_entrez,
  OrgDb = org.Dm.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)
dotplot(ego_post2, showCategory = 15) + ggtitle("GO enriquecido - Población Post")

ggsave("GO_Post_PrevsPost.png", 
       plot = dotplot(ego_post2, showCategory = 15) + 
         ggtitle("GO enriquecido - Población Post"),
       width = 10, height = 6, dpi = 300)

ego_post2_df <- as.data.frame(ego_post2)
Tabla_GOPost_prevspost <- ego_post2_df %>%
  mutate(
    N_genes_enriquecidos = str_count(geneID, "/") + 1,  # Cuenta los genes en geneID
    Genes_representativos = sapply(strsplit(geneID, "/"), function(x) paste(head(x, 5), collapse = ", "))
  ) %>%
  select(Description, N_genes_enriquecidos, Count, p.adjust, Genes_representativos)

write.csv(Tabla_GOPost_prevspost, "Tabla_GOPost_prevspost.csv", row.names = FALSE)
