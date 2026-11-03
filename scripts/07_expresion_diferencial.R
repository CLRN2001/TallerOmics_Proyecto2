# =============================================================================
# Paso 7: Expresión diferencial (DEG) — MS vs. control por tipo celular
# Método: MAST (como el paper Schafflick et al. 2020)
# Proyecto: MS_scRNAseq_replication
# =============================================================================

# NOTA: no se puede probar este código antes de entregarlo; puede requerir
# pequeños ajustes según el estado exacto del objeto guardado.

suppressPackageStartupMessages({
  library(Seurat)
  library(MAST)        
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

# ── Parámetros ────────────────────────────────────────────────────────────────
PVAL_CORTE <- 0.05   # FDR máximo para DEG significativo
FC_CORTE   <- 0.25   # |log2FC| mínimo (umbral estándar de Seurat)
MIN_PCT    <- 0.10   # gen debe expresarse en ≥10% de las células de algún grupo

# ── Directorios de salida ─────────────────────────────────────────────────────
dir.create("results/DEG/CSF",   recursive = TRUE, showWarnings = FALSE)
dir.create("results/DEG/blood", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/DEG",       showWarnings = FALSE)

# ── 1. Cargar objeto ──────────────────────────────────────────────────────────
cat("Cargando objeto anotado...\n")
ms <- readRDS("data/processed/05_ms_anotado.rds")

DefaultAssay(ms) <- "RNA"
Idents(ms)       <- "celltype_nuestro"

cat("Células totales:", ncol(ms), "\n")
cat("Tipos celulares:", paste(sort(unique(ms$celltype_nuestro)), collapse = ", "), "\n\n")

# ── 2. Tipos celulares elegibles ──────────────────────────────────────────────
# LCR: excluidos B1 (control=7), Gran/Mono (control=1), MegaK (control=4)
tipos_CSF <- c("B2", "CD4", "CD4/Tregs", "CD8a", "CD8a/Tdg", "CD8n",
               "Gran", "Mono", "ncMono", "NK1", "NK2", "pDC",
               "plasma/mDC2", "Tregs")

# Sangre: excluido ncMono (control=4); Gran/Mono incluido con advertencia
tipos_blood <- c("B1", "B2", "CD4", "CD4/Tregs", "CD8a", "CD8a/Tdg", "CD8n",
                 "Gran", "Gran/Mono", "MegaK", "Mono", "NK1", "NK2", "pDC",
                 "plasma/mDC2", "Tregs")

# ── 3. Función principal de DEG ───────────────────────────────────────────────
correr_DEG <- function(ms_tejido, tejido, tipos, outdir) {
  
  resultados <- list()
  
  for (ct in tipos) {
    cat(sprintf("  [%s] Procesando: %s\n", tejido, ct))
    
    # Subset a este tipo celular
    cells_use <- rownames(ms_tejido@meta.data)[ms_tejido$celltype_nuestro == ct]
    
    if (length(cells_use) < 20) {
      cat(sprintf("    -> Saltando (solo %d células totales)\n", length(cells_use)))
      next
    }
    
    obj_ct <- subset(ms_tejido, cells = cells_use)
    
    n_ctrl <- sum(obj_ct$condicion == "control")
    n_ms   <- sum(obj_ct$condicion == "MS")
    cat(sprintf("    -> control=%d, MS=%d\n", n_ctrl, n_ms))
    
    if (n_ctrl < 10 || n_ms < 10) {
      cat(sprintf("    -> Saltando (grupo demasiado pequeño)\n"))
      next
    }
    
    # Correr MAST
    tryCatch({
      Idents(obj_ct) <- "condicion"
      
      deg <- FindMarkers(
        object          = obj_ct,
        ident.1         = "MS",
        ident.2         = "control",
        test.use        = "MAST",
        min.pct         = MIN_PCT,
        logfc.threshold = FC_CORTE,
        assay           = "RNA",
        latent.vars     = "nCount_RNA"  # controla profundidad de secuenciación
      )
      
      # Añadir columnas informativas
      deg$gene      <- rownames(deg)
      deg$celltype  <- ct
      deg$tejido    <- tejido
      deg$n_control <- n_ctrl
      deg$n_MS      <- n_ms
      
      # Ordenar por FDR
      deg <- deg[order(deg$p_val_adj), ]
      
      # Guardar CSV individual
      nombre_limpio <- gsub("[/]", "_", ct)
      outfile <- file.path(outdir, paste0("DEG_", nombre_limpio, ".csv"))
      write.csv(deg, file = outfile, row.names = FALSE)
      
      n_sig <- sum(deg$p_val_adj < PVAL_CORTE, na.rm = TRUE)
      n_up  <- sum(deg$p_val_adj < PVAL_CORTE & deg$avg_log2FC > 0, na.rm = TRUE)
      n_dn  <- sum(deg$p_val_adj < PVAL_CORTE & deg$avg_log2FC < 0, na.rm = TRUE)
      cat(sprintf("    -> %d DEGs significativos (%d↑MS, %d↓MS)\n",
                  n_sig, n_up, n_dn))
      
      resultados[[ct]] <- deg
      
    }, error = function(e) {
      cat(sprintf("    -> ERROR en %s: %s\n", ct, conditionMessage(e)))
    })
  }
  
  return(resultados)
}

# ── 4. Correr DEG en LCR ──────────────────────────────────────────────────────
cat("=== DEG en LCR (CSF) — puede tardar 2-3 horas ===\n")
ms_CSF  <- subset(ms, subset = tejido == "CSF")
res_CSF <- correr_DEG(ms_CSF, "CSF", tipos_CSF, "results/DEG/CSF")
cat("LCR terminado.\n\n")

# ── 5. Correr DEG en sangre ───────────────────────────────────────────────────
cat("=== DEG en sangre (blood) — puede tardar 2-3 horas ===\n")
ms_blood  <- subset(ms, subset = tejido == "blood")
res_blood <- correr_DEG(ms_blood, "blood", tipos_blood, "results/DEG/blood")
cat("Sangre terminada.\n\n")

# ── 6. Consolidar resultados ──────────────────────────────────────────────────
consolidar <- function(lista_res, outfile_top, outfile_all) {
  if (length(lista_res) == 0) {
    cat("  (sin resultados para consolidar)\n")
    return(invisible(NULL))
  }
  
  df_all <- do.call(rbind, lista_res)
  rownames(df_all) <- NULL
  
  # Guardar tabla completa
  write.csv(df_all, file = outfile_all, row.names = FALSE)
  cat("  Tabla completa guardada:", outfile_all, "\n")
  
  # Top 10 DEGs por tipo celular (por |log2FC|, solo significativos)
  top_genes <- df_all %>%
    filter(p_val_adj < PVAL_CORTE) %>%
    group_by(celltype) %>%
    slice_max(order_by = abs(avg_log2FC), n = 10) %>%
    select(celltype, gene, avg_log2FC, p_val_adj, pct.1, pct.2,
           n_MS, n_control) %>%
    arrange(celltype, desc(abs(avg_log2FC)))
  
  write.csv(top_genes, file = outfile_top, row.names = FALSE)
  cat("  Top genes guardados:    ", outfile_top, "\n")
  
  return(df_all)
}

cat("Consolidando resultados...\n")
df_CSF   <- consolidar(res_CSF,
                       "results/DEG/top_DEG_CSF.csv",
                       "results/DEG/todos_DEG_CSF.csv")
df_blood <- consolidar(res_blood,
                       "results/DEG/top_DEG_blood.csv",
                       "results/DEG/todos_DEG_blood.csv")

# ── 7. Volcano plots ──────────────────────────────────────────────────────────
hacer_volcano <- function(df_all, tejido_label, outfile, ncol_facet = 4) {
  if (is.null(df_all) || nrow(df_all) == 0) {
    cat("  (sin datos para volcano de", tejido_label, ")\n")
    return(invisible(NULL))
  }
  
  df_all <- df_all %>%
    mutate(
      p_plot = pmin(p_val_adj, 1e-300),   # evitar -log10(0)
      sig = p_val_adj < PVAL_CORTE & abs(avg_log2FC) > FC_CORTE,
      direccion = case_when(
        sig & avg_log2FC > 0  ~ "Subido en MS",
        sig & avg_log2FC < 0  ~ "Bajado en MS",
        TRUE                  ~ "NS"
      )
    )
  
  # Etiquetar top 5 genes por tipo celular (mayor |FC|)
  genes_label <- df_all %>%
    filter(sig) %>%
    group_by(celltype) %>%
    slice_max(order_by = abs(avg_log2FC), n = 5) %>%
    pull(gene)
  
  df_all$label <- ifelse(df_all$gene %in% genes_label, df_all$gene, "")
  
  p <- ggplot(df_all,
              aes(x = avg_log2FC,
                  y = -log10(p_plot),
                  color = direccion,
                  label = label)) +
    geom_point(alpha = 0.35, size = 0.7) +
    geom_text_repel(
      data       = subset(df_all, label != ""),
      size       = 2.5,
      max.overlaps = 15,
      segment.size = 0.3
    ) +
    scale_color_manual(
      values = c("Subido en MS"  = "#E41A1C",
                 "Bajado en MS" = "#377EB8",
                 "NS"           = "grey75")
    ) +
    geom_vline(xintercept = c(-FC_CORTE, FC_CORTE),
               linetype = "dashed", color = "grey40", linewidth = 0.4) +
    geom_hline(yintercept = -log10(PVAL_CORTE),
               linetype = "dashed", color = "grey40", linewidth = 0.4) +
    facet_wrap(~ celltype, scales = "free_y", ncol = ncol_facet) +
    labs(
      title    = paste("DEG —", tejido_label, "— MS vs. control (MAST)"),
      subtitle = paste("FDR <", PVAL_CORTE, "  |log2FC| >", FC_CORTE),
      x        = "log2 Fold Change (MS / control)",
      y        = "-log10(FDR)",
      color    = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "bottom",
      strip.text      = element_text(size = 8, face = "bold"),
      panel.grid.minor = element_blank()
    )
  
  ggsave(outfile, p, width = 16, height = 14, dpi = 150)
  cat("  Volcano guardado:", outfile, "\n")
}

cat("Generando volcano plots...\n")
hacer_volcano(df_CSF,   "LCR",    "figures/DEG/volcano_CSF.png",   ncol_facet = 4)
hacer_volcano(df_blood, "Sangre", "figures/DEG/volcano_blood.png", ncol_facet = 4)

# ── 8. Resumen por consola ────────────────────────────────────────────────────
cat("\n")
cat("========================================\n")
cat("RESUMEN FINAL — Paso 7 (DEG)\n")
cat("========================================\n")

for (info in list(list(res = res_CSF, label = "LCR"),
                  list(res = res_blood, label = "Sangre"))) {
  cat(sprintf("\n[%s]\n", info$label))
  if (length(info$res) == 0) {
    cat("  (sin resultados)\n")
    next
  }
  for (ct in names(info$res)) {
    d <- info$res[[ct]]
    n_sig <- sum(d$p_val_adj < PVAL_CORTE, na.rm = TRUE)
    n_up  <- sum(d$p_val_adj < PVAL_CORTE & d$avg_log2FC > 0, na.rm = TRUE)
    n_dn  <- sum(d$p_val_adj < PVAL_CORTE & d$avg_log2FC < 0, na.rm = TRUE)
    cat(sprintf("  %-15s: %4d DEGs (%3d↑MS  %3d↓MS)\n",
                ct, n_sig, n_up, n_dn))
  }
}

cat("\n=== Paso 7 completado ===\n")
cat("Archivos en results/DEG/ y figures/DEG/\n")