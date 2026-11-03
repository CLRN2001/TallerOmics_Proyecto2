# ============================================================
# 03_normalizacion_integracion.R — Normalización, PCA, Harmony y UMAP
# Replicación: Schafflick et al. 2020 (GSE138266)
# Nota: el paper integró con scVI (Python); aquí se usa Harmony (R).
# ============================================================
library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)

ms <- readRDS("data/processed/02_ms_paper.rds")

# --- 1. Quedarse con las células válidas según los autores ---
ms <- subset(ms, subset = valid_paper == "True")
cat("Células válidas:", ncol(ms), "(esperado: 61051)\n")
print(table(ms$tejido, ms$condicion))

# --- 2. Normalización (log de conteos por 10.000) ---
ms <- NormalizeData(ms, normalization.method = "LogNormalize", scale.factor = 1e4)

# --- 3. Genes variables, escalado y PCA ---
ms <- FindVariableFeatures(ms, selection.method = "vst", nfeatures = 2000)
ms <- ScaleData(ms)                       # solo los 2.000 genes variables
ms <- RunPCA(ms, npcs = 50, verbose = FALSE)

p_elbow <- ElbowPlot(ms, ndims = 50) + ggtitle("Varianza explicada por PC")
ggsave("results/clustering/elbow_plot.png", p_elbow, width = 6, height = 4)

# --- 4. UMAP SIN integrar (para mostrar el efecto de lote) ---
n_dims <- 20     # AJUSTAR según el elbow plot
ms <- RunUMAP(ms, reduction = "pca", dims = 1:n_dims,
              reduction.name = "umap_sin_integrar", verbose = FALSE)

# --- 5. Integración con Harmony por lote (= muestra) ---
ms <- RunHarmony(ms, group.by.vars = "lote", reduction.use = "pca",
                 dims.use = 1:n_dims, verbose = FALSE)
ms <- RunUMAP(ms, reduction = "harmony", dims = 1:n_dims,
              reduction.name = "umap", verbose = FALSE)

# --- 6. Figuras: antes vs. después de integrar ---
p_lote <- DimPlot(ms, reduction = "umap_sin_integrar", group.by = "lote") +
  ggtitle("Sin integrar (por lote)") +
  DimPlot(ms, reduction = "umap", group.by = "lote") +
  ggtitle("Con Harmony (por lote)")
ggsave("results/clustering/umap_efecto_lote.png", p_lote, width = 14, height = 6)

# --- 7. Figuras: tu UMAP vs. el del paper ---
p_comp <- DimPlot(ms, reduction = "umap", group.by = "celltype_paper",
                  label = TRUE, repel = TRUE) + ggtitle("Nuestro UMAP (Harmony)") + NoLegend() +
  DimPlot(ms, reduction = "umap_paper", group.by = "celltype_paper",
          label = TRUE, repel = TRUE) + ggtitle("UMAP del paper (scVI)") + NoLegend()
ggsave("figures/umap_nuestro_vs_paper.png", p_comp, width = 14, height = 6)

p_meta <- DimPlot(ms, reduction = "umap", group.by = "tejido") +
  DimPlot(ms, reduction = "umap", group.by = "condicion")
ggsave("figures/umap_tejido_condicion.png", p_meta, width = 14, height = 6)

saveRDS(ms, "data/processed/03_ms_integrado.rds")
message("Paso 3 terminado: data/processed/03_ms_integrado.rds")