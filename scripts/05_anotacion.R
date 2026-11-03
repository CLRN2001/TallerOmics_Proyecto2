# ============================================================
# 05_anotacion.R — Anotación celular
# Replicación: Schafflick et al. 2020 (GSE138266)
# ============================================================
library(Seurat)
library(ggplot2)
library(patchwork)

ms <- readRDS("data/processed/04_ms_clusterizado.rds")

col_elegida <- grep("RNA_snn_res.0.8", colnames(ms@meta.data), value = TRUE)[1]
Idents(ms)  <- ms@meta.data[[col_elegida]]
cat("Niveles de Idents:", levels(Idents(ms)), "\n")   # confirmar qué valores tiene

# --- 1. Genes marcadores ---
# Solo genes que existen en el dataset
genes_dataset <- rownames(ms)

marcadores_candidatos <- list(
  "CD4 T"  = c("CD3D", "CD4", "IL7R"),
  "CD8 T"  = c("CD3D", "CD8A", "CD8B"),
  "Tregs"  = c("FOXP3", "IL2RA", "CTLA4"),
  "Tdg"    = c("TRDC", "TRGC1", "TRGC2"),
  "NK"     = c("NKG7", "GNLY", "KLRD1"),
  "B"      = c("MS4A1", "CD79A"),
  "plasma" = c("MZB1", "JCHAIN"),
  "Mono"   = c("CD14", "LYZ", "S100A8"),
  "ncMono" = c("FCGR3A", "MS4A7", "CX3CR1"),
  "mDC"    = c("FCER1A", "CD1C", "CLEC10A"),
  "pDC"    = c("CLEC4C", "IL3RA", "LILRA4"),
  "Gran"   = c("S100A8", "S100A9", "CXCR1", "CSF3R"),
  "MegaK"  = c("ITGA2B", "GP9", "VWF")
)

# Filtrar solo los que existen en el dataset
marcadores <- lapply(marcadores_candidatos, function(g) g[g %in% genes_dataset])
marcadores <- marcadores[sapply(marcadores, length) > 0]

genes_no_encontrados <- setdiff(unlist(marcadores_candidatos), genes_dataset)
cat("Genes no disponibles en el dataset:", paste(genes_no_encontrados, collapse = ", "), "\n")

genes_plot <- unique(unlist(marcadores))
cat("Genes disponibles para el dotplot:", length(genes_plot), "\n")

p_dot <- DotPlot(ms, features = genes_plot, group.by = col_elegida) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7)) +
  labs(title = "Expresión de marcadores por cluster")
ggsave("results/clustering/dotplot_marcadores.png", p_dot, width = 14, height = 10)

# --- 2. Asignación de nombres ---
# Ajusta según lo que mostraron el heatmap (Paso 4) y el dotplot
etiquetas <- c(
  "0"  = "CD4",
  "1"  = "CD4/Tregs",
  "2"  = "CD8a",
  "3"  = "Gran",
  "4"  = "CD8a/Tdg",
  "5"  = "NK1",
  "6"  = "CD8n",
  "7"  = "Tregs",
  "8"  = "B1",
  "9"  = "Mono",
  "10" = "B2",
  "11" = "ncMono",
  "12" = "Mono",
  "13" = "NK2",
  "14" = "MegaK",
  "15" = "pDC",
  "16" = "plasma/mDC2",
  "17" = "Gran",
  "18" = "Gran/Mono"
)

# Asignación segura: nombrar por código de barras de célula
cluster_por_celula  <- as.character(Idents(ms))
etiqueta_por_celula <- etiquetas[cluster_por_celula]
names(etiqueta_por_celula) <- colnames(ms)   # clave: nombres = códigos de barras
ms <- AddMetaData(ms, metadata = etiqueta_por_celula, col.name = "celltype_nuestro")

cat("Distribución de tipos anotados:\n")
print(table(ms$celltype_nuestro))
cat("NAs en anotación:", sum(is.na(ms$celltype_nuestro)), "\n")

# --- 3. UMAP final ---
p_final <- DimPlot(ms, reduction = "umap", group.by = "celltype_nuestro",
                   label = TRUE, repel = TRUE) + ggtitle("Nuestra anotación") + NoLegend() +
  DimPlot(ms, reduction = "umap", group.by = "celltype_paper",
          label = TRUE, repel = TRUE) + ggtitle("Anotación del paper") + NoLegend()
ggsave("figures/umap_anotacion_final.png", p_final, width = 14, height = 6)

saveRDS(ms, "data/processed/05_ms_anotado.rds")
message("Paso 5 terminado: data/processed/05_ms_anotado.rds")