# ============================================================
# PASO 8 — GSEA (Gene Set Enrichment Analysis)
# Input:  todos_DEG_CSF.csv  y  todos_DEG_blood.csv
# Output: gsea_CSF.csv, gsea_blood.csv,
#         dotplot_gsea_CSF.png, dotplot_gsea_blood.png
# ============================================================

library(dplyr)
library(fgsea)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

# Indicar ruta deseada para los resultados setwd("")

# ── 1. Cargar DEGs ──────────────────────────────────────────
deg_csf   <- read.csv("todos_DEG_CSF.csv",   stringsAsFactors = FALSE)
deg_blood <- read.csv("todos_DEG_blood.csv", stringsAsFactors = FALSE)

# ── 2. Función principal ─────────────────────────────────────
run_gsea <- function(deg_df, tejido) {
  
  tipos <- unique(deg_df$celltype)
  resultados <- list()
  
  for (ct in tipos) {
    
    df <- deg_df %>% filter(celltype == ct)
    
    # Estadístico de ranking: sign(log2FC) × -log10(FDR)
    df <- df %>%
      mutate(rank_stat = sign(avg_log2FC) * (-log10(p_val_adj + 1e-300))) %>%
      arrange(desc(rank_stat))
    
    # Convertir símbolos a Entrez ID
    genes_entrez <- bitr(df$gene,
                         fromType = "SYMBOL",
                         toType   = "ENTREZID",
                         OrgDb    = org.Hs.eg.db)
    
    df <- df %>%
      inner_join(genes_entrez, by = c("gene" = "SYMBOL"))
    
    # Vector nombrado para fgsea
    ranked <- setNames(df$rank_stat, df$ENTREZID)
    
    # Cargar gene sets desde archivo local GMT
    go_sets_df <- fgsea::gmtPathways("D:/Lelu/Documents/Chopo/MS_scRNAseq_replication/datos/c5.go.bp.v2023.2.Hs.symbols.gmt")
    
    # Vector nombrado con símbolos directamente
    ranked <- setNames(df$rank_stat, df$gene)
    
    # Correr fgsea con gene sets locales
    set.seed(42)
    res <- fgsea(pathways = go_sets_df,
                 stats    = ranked,
                 minSize  = 15,
                 maxSize  = 500,
                 nPermSimple = 1000)
    
    res$celltype <- ct
    res$tejido   <- tejido
    resultados[[ct]] <- res %>% filter(padj < 0.05) %>% arrange(padj)
  }
  
  bind_rows(resultados)
}

# ── 3. Correr para cada tejido ───────────────────────────────
cat("Corriendo GSEA en LCR...\n")
gsea_csf   <- run_gsea(deg_csf,   "CSF")

cat("Corriendo GSEA en sangre...\n")
gsea_blood <- run_gsea(deg_blood, "Blood")

# ── 4. Guardar CSVs ──────────────────────────────────────────
# Convertir leadingEdge de lista a texto antes de guardar
gsea_csf$leadingEdge <- sapply(gsea_csf$leadingEdge, paste, collapse = ";")
gsea_blood$leadingEdge <- sapply(gsea_blood$leadingEdge, paste, collapse = ";")

write.csv(gsea_csf,   "gsea_CSF.csv",   row.names = FALSE)
write.csv(gsea_blood, "gsea_blood.csv", row.names = FALSE) 

# ── 5. Dotplot LCR ───────────────────────────────────────────
top_csf <- gsea_csf %>%
  group_by(celltype) %>%
  slice_min(padj, n = 5) %>%
  ungroup() %>%
  mutate(pathway = stringr::str_trunc(pathway, 40))

p1 <- ggplot(top_csf, aes(x = celltype, y = pathway,
                          size = -log10(padj),
                          color = NES)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red",
                        midpoint = 0, name = "NES") +
  scale_size_continuous(name = "-log10(FDR)") +
  labs(title = "GSEA — LCR (Top 5 vías por tipo celular)",
       x = "Tipo celular", y = "Vía biológica (GO:BP)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 7))

ggsave("dotplot_gsea_CSF.png", p1, width = 14, height = 8, dpi = 150)

# ── 6. Dotplot sangre ────────────────────────────────────────
top_blood <- gsea_blood %>%
  group_by(celltype) %>%
  slice_min(padj, n = 5) %>%
  ungroup() %>%
  mutate(pathway = stringr::str_trunc(pathway, 40))

p2 <- ggplot(top_blood, aes(x = celltype, y = pathway,
                            size = -log10(padj),
                            color = NES)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red",
                        midpoint = 0, name = "NES") +
  scale_size_continuous(name = "-log10(FDR)") +
  labs(title = "GSEA — Sangre (Top 5 vías por tipo celular)",
       x = "Tipo celular", y = "Vía biológica (GO:BP)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 7))

ggsave("dotplot_gsea_blood.png", p2, width = 14, height = 8, dpi = 150)

cat("Paso 8 completado.\n")