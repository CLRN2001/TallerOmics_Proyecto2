# ============================================================
# 02_cargar_datos_paper.R — Objeto Seurat desde el .h5ad + metadatos del GitHub
# Replicación: Schafflick et al. 2020 (GSE138266)
# ============================================================
library(rhdf5)
library(Matrix)
library(Seurat)
library(RcppCNPy)
library(ggplot2)

f    <- "data/reference/MS_CSF.h5ad"
meta <- "github_paper/MSscRNAseq2019/Notebooks/meta"

# --- 1. Leer una matriz dispersa del .h5ad (células x genes en CSR = genes x células en CSC) ---
leer_matriz <- function(grupo, n_genes, n_cells) {
  sparseMatrix(i = as.integer(h5read(f, paste0(grupo, "/indices"))) + 1,
               p = as.integer(h5read(f, paste0(grupo, "/indptr"))),
               x = as.numeric(h5read(f, paste0(grupo, "/data"))),
               dims = c(n_genes, n_cells))
}

celulas <- paste0("cell_", h5read(f, "/obs/_index"))
genes   <- as.character(h5read(f, "/var/_index"))
X <- leer_matriz("/X", length(genes), length(celulas))
dimnames(X) <- list(genes, celulas)

# Si /raw tiene conteos, se usan esos
raw_ok <- tryCatch({
  rv <- h5read(f, "/raw/var")
  genes_raw <- as.character(if (is.data.frame(rv)) rv[[1]] else rv[["_index"]])
  R <- leer_matriz("/raw/X", length(genes_raw), length(celulas))
  dimnames(R) <- list(make.unique(genes_raw), celulas)
  TRUE
}, error = function(e) FALSE)

counts <- if (raw_ok) R else X
rm(X); if (raw_ok) rm(R); invisible(gc())
cat("Matriz usada:", ifelse(raw_ok, "/raw/X", "/X"),
    "| genes:", nrow(counts), "| células:", ncol(counts),
    "| ¿conteos enteros?:", all(counts@x == round(counts@x)), "\n")

# --- 1b. Recuperar conteos enteros: cada célula = conteo × factor ---
col_de_x <- rep(seq_len(ncol(counts)), diff(counts@p))   # columna de cada valor no cero

min_por_col <- tapply(counts@x, col_de_x, min)
factor_cel  <- rep(NA_real_, ncol(counts))
factor_cel[as.integer(names(min_por_col))] <- as.numeric(min_por_col)
cat("Células sin valores:", sum(is.na(factor_cel)), "\n")   # esperado: 0

counts@x <- as.numeric(counts@x / factor_cel[col_de_x])
cat("Factores por célula:\n"); print(table(round(factor_cel, 3)))
cat("¿Conteos enteros tras corregir?:", all(abs(counts@x - round(counts@x)) < 1e-6), "\n")
counts@x <- round(counts@x)

# --- 2. Metadatos del .h5ad ---
leer_cat <- function(col) {
  codigos <- as.vector(h5read(f, paste0("/obs/", col)))
  niveles <- tryCatch(as.character(h5read(f, paste0("/obs/__categories/", col))),
                      error = function(e) NULL)
  if (is.null(niveles)) codigos else niveles[codigos + 1]
}
md <- data.frame(
  tejido         = ifelse(as.vector(h5read(f, "/obs/CSF")) == 1, "CSF", "blood"),
  condicion      = ifelse(as.vector(h5read(f, "/obs/MS"))  == 1, "MS", "control"),
  celltype_paper = leer_cat("labels"),
  cd4_paper      = leer_cat("CD4_labels"),
  valid_paper    = leer_cat("valid_clusters"),
  row.names      = celulas
)

# --- 3. Crear objeto Seurat ---
ms <- CreateSeuratObject(counts = counts, meta.data = md, project = "Schafflick2020")

# --- 3b. Metadatos del GitHub (mismo orden de células que el .h5ad) ---
leer_npy <- function(archivo) {
  ruta <- file.path(meta, archivo)
  if (!file.exists(ruta)) stop("No existe el archivo: ", normalizePath(ruta, mustWork = FALSE))
  x <- as.vector(npyLoad(ruta, type = "numeric"))
  if (length(x) != ncol(ms)) stop(archivo, " tiene ", length(x), " valores; se esperaban ", ncol(ms))
  # Enteros de 64 bits leídos como decimales diminutos -> recuperar el entero
  if (all(abs(x) < 1e-300)) x <- round(x / 4.94065645841247e-324)
  x
}
ms$lote          <- factor(leer_npy("batchid.npy"))
ms$louvain_paper <- factor(leer_npy("louvain.k30.allgenes.npy"))
ms$doublet_score <- leer_npy("doublet_score.npy")

# --- 4. UMAP del paper ---
umap <- t(h5read(f, "/obsm/X_umap"))
dimnames(umap) <- list(celulas, c("umapPaper_1", "umapPaper_2"))
ms[["umap_paper"]] <- CreateDimReducObject(embeddings = umap, key = "umapPaper_", assay = "RNA")

# --- 5. Comprobaciones ---
print(table(ms$tejido, ms$condicion))              # esperado: 25831 / 17138 / 9652 / 12705
print(table(ms$lote, ms$tejido))                   # lotes 0 a 17
print(table(ms$celltype_paper, ms$valid_paper))    # ¿doublet/RBC/contamination = no válidas?
print(summary(ms$nCount_RNA)); print(summary(ms$nFeature_RNA))
print(summary(ms$doublet_score))

p <- DimPlot(ms, reduction = "umap_paper", group.by = "celltype_paper", label = TRUE, repel = TRUE) +
  DimPlot(ms, reduction = "umap_paper", group.by = "tejido")
ggsave("figures/umap_paper_referencia.png", p, width = 14, height = 6)

saveRDS(ms, "data/processed/02_ms_paper.rds")
message("Paso 2 terminado: data/processed/02_ms_paper.rds")