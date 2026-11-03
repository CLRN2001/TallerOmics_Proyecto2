# ============================================================
# 04_clustering.R — Clustering propio y comparación con el paper
# Replicación: Schafflick et al. 2020 (GSE138266)
# ============================================================
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

ms <- readRDS("data/processed/03_ms_integrado.rds")

# --- 1. Grafo de vecinos (k=30, igual que el paper) ---
ms <- FindNeighbors(ms, reduction = "harmony", dims = 1:20,
                    k.param = 30, verbose = FALSE)

# --- 2. Clustering a varias resoluciones ---
for (res in c(0.3, 0.5, 0.8, 1.2)) {
  ms <- FindClusters(ms, algorithm = 1, resolution = res, verbose = FALSE)
}

# Ver qué columnas se crearon realmente y cuántos clusters produce cada una
cols_res <- grep("RNA_snn_res", colnames(ms@meta.data), value = TRUE)
cat("Columnas creadas:\n"); print(cols_res)
for (col in cols_res) {
  cat("  ", col, "->", nlevels(ms@meta.data[[col]]), "clusters\n")
}

# --- 3. UMAP coloreado por las primeras 3 resoluciones ---
p_res <- lapply(cols_res[1:min(3, length(cols_res))], function(col) {
  DimPlot(ms, reduction = "umap", group.by = col, label = TRUE) +
    ggtitle(col) + NoLegend()
})
ggsave("results/clustering/umap_resoluciones.png",
       wrap_plots(p_res, ncol = 3), width = 18, height = 6)

# --- 4. Elegir la resolución con más clusters cercano a 21 (tipos del paper) ---
# Cambiar "0.8" por la resolución que elijas según la tabla de arriba
col_elegida <- cols_res[grep("0.8", cols_res)[1]]
cat("Resolución elegida:", col_elegida,
    "->", nlevels(ms@meta.data[[col_elegida]]), "clusters\n")
Idents(ms) <- ms@meta.data[[col_elegida]]

# --- 5. Tabla de confusión: clusters propios vs. tipos del paper ---
tabla <- table(cluster_nuestro = Idents(ms),
               tipo_paper      = ms$celltype_paper)

write.csv(as.data.frame.matrix(tabla),
          "results/clustering/confusion_nuestro_vs_paper.csv")

# Heatmap proporcional (cada columna suma 1)
prop <- prop.table(tabla, margin = 2)
df   <- as.data.frame(prop)
colnames(df) <- c("cluster_nuestro", "tipo_paper", "Freq")

p_conf <- ggplot(df, aes(tipo_paper, cluster_nuestro, fill = Freq)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "#2166ac", name = "Proporción") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  labs(title = "Correspondencia: clusters propios vs. tipos del paper",
       x = "Tipo celular (paper)", y = "Cluster propio")
ggsave("results/clustering/heatmap_confusion.png", p_conf, width = 12, height = 8)

saveRDS(ms, "data/processed/04_ms_clusterizado.rds")
message("Paso 4 terminado: data/processed/04_ms_clusterizado.rds")