# ============================================================
# 06_abundancia.R — Abundancia celular diferencial (Beta-Binomial)
# Replicación: Schafflick et al. 2020 (GSE138266)
# Usa BetaBinom.R del GitHub del paper
# ============================================================
library(Seurat)
library(ggplot2)
library(dplyr)
library(aod)
library(VGAM)

ms <- readRDS("data/processed/05_ms_anotado.rds")

# --- 1. Definir orden de muestras igual que el paper ---
# BetaBinom.R espera: primero 8 LCR (4 MS + 4 ctrl), luego 10 sangre (5 MS + 5 ctrl)
orden_lotes <- data.frame(
  lote      = factor(c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17)),
  tejido    = c(rep("CSF",8), rep("blood",10)),
  condicion = c("MS","MS","MS","MS","control","control","control","control",
                "MS","MS","MS","MS","MS","control","control","control","control","control")
)

# Verificar que coincide con tus datos
print(table(ms$lote, ms$tejido))

# --- 2. Construir la tabla de conteos ---
# Filas = tipos celulares, columnas = muestras en el orden del paper
tipos <- sort(unique(ms$celltype_nuestro))

count <- sapply(orden_lotes$lote, function(l) {
  celulas_lote <- ms$celltype_nuestro[ms$lote == l]
  sapply(tipos, function(t) sum(celulas_lote == t))
})
rownames(count) <- tipos
colnames(count) <- paste0("lote_", orden_lotes$lote)

cat("Tabla de conteos (tipos x muestras):\n")
print(count)
cat("\nTotal por muestra:\n")
print(colSums(count))

# --- 3. Función BetaBinomialReg del paper (extraída de BetaBinom.R) ---
# No se usa source() porque el script original requiere argumentos de línea de comandos
library(VGAM)
library(aod)

BetaBinomialReg <- function(count, N, C, filename) {
  res <- lapply(seq_len(nrow(count)), function(i) {
    x    <- as.numeric(count[i, ])
    data <- data.frame(x = x, N = N, C = as.factor(C))
    fm1  <- betabin(cbind(x, N - x) ~ C, ~ 1, data)
    r    <- summary(fm1)@Coef
    p    <- r[2, 4]
    fc   <- mean(x[C == TRUE] / N[C == TRUE]) /
      mean(x[C == FALSE] / N[C == FALSE])
    c(p, fc)
  })
  res <- do.call(cbind, res)
  write.table(t(res), col.names = FALSE, row.names = FALSE,
              quote = FALSE, sep = ",", file = filename)
  return(res)
}

# --- 4. Definir vectores de tejido y estado (mismo orden que count) ---
N      <- as.numeric(colSums(count))
tissue <- orden_lotes$tejido
states <- orden_lotes$condicion

# --- 5. Correr las 4 comparaciones relevantes ---
dir.create("results/DEG", showWarnings = FALSE)

# Comparación 3: MS vs. control en LCR (la más importante del paper)
res_CSF <- BetaBinomialReg(
  count  = count[, tissue == "CSF"],
  N      = N[tissue == "CSF"],
  C      = states[tissue == "CSF"] == "MS",
  filename = "results/DEG/betabinom_CSF_MSvsCtrl.csv"
)

# Comparación 4: MS vs. control en sangre
res_blood <- BetaBinomialReg(
  count  = count[, tissue == "blood"],
  N      = N[tissue == "blood"],
  C      = states[tissue == "blood"] == "MS",
  filename = "results/DEG/betabinom_blood_MSvsCtrl.csv"
)

# --- 6. Formatear resultados ---
formatear <- function(res, tipos, tejido_label) {
  df <- data.frame(
    tipo_celular = tipos,
    p_valor      = as.numeric(res[1, ]),
    fold_change  = as.numeric(res[2, ]),
    tejido       = tejido_label
  )
  df$p_adj     <- p.adjust(df$p_valor, method = "BH")
  df$log2FC    <- log2(df$fold_change)
  df$sig       <- df$p_adj < 0.05
  df[order(df$p_adj), ]
}

res_CSF_df   <- formatear(res_CSF,   tipos, "CSF")
res_blood_df <- formatear(res_blood, tipos, "blood")

print("--- LCR: MS vs. control ---")
print(res_CSF_df)
print("--- Sangre: MS vs. control ---")
print(res_blood_df)

write.csv(res_CSF_df,   "results/DEG/abundancia_CSF.csv",   row.names = FALSE)
write.csv(res_blood_df, "results/DEG/abundancia_blood.csv", row.names = FALSE)

# --- 7. Figura: volcano plot de abundancia ---
df_plot <- rbind(res_CSF_df, res_blood_df)

p_volc <- ggplot(df_plot, aes(x = log2FC, y = -log10(p_adj + 1e-10),
                              color = sig, label = tipo_celular)) +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(data = subset(df_plot, sig),
                           size = 3, max.overlaps = 20) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("grey70", "#d7191c"),
                     labels = c("No sig.", "FDR < 5 %")) +
  facet_wrap(~tejido) +
  labs(title = "Abundancia celular diferencial: MS vs. control",
       x = "log2(Fold Change)", y = "-log10(p ajustado)",
       color = NULL) +
  theme_bw()
ggsave("figures/abundancia_volcano.png", p_volc, width = 12, height = 6)

message("Paso 6 terminado")