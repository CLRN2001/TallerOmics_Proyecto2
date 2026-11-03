# ==============================================================================
# Script:      01_descarga.R
# Propósito:   Descargar y organizar los datos crudos (matrices 10x) de GEO y el
#              objeto integrado de referencia de Figshare.
# Paper:       Schafflick D, Xu CA, et al. (2020). "Integrated single cell
#              analysis of blood and cerebrospinal fluid leukocytes in multiple
#              sclerosis". Nature Communications 11:247.
#              doi:10.1038/s41467-019-14118-w
# Datos:       GEO GSE138266 (22 muestras: 12 LCR + 10 sangre; 6 EM + 6 IIH)
#              Figshare 14356661 (MS_CSF.h5ad, ~1 GB)
# Fecha:       2026-10-13
# Paquetes:    GEOquery (Bioconductor), tools (base R), here (opcional)
# Salida:      data/raw/<GSM>_<tejido>_<condicion>_<paciente>/
#                  barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz
#              data/raw/metadata_muestras.csv
#              data/reference/MS_CSF.h5ad
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. Configuración general
# ------------------------------------------------------------------------------
# En Windows, download.file corta las descargas a los 60 s por defecto.
# Subimos el límite a 1 hora para el archivo de Figshare (~1 GB).
options(timeout = 3600)

suppressPackageStartupMessages({
  library(GEOquery)
  library(tools)     # md5sum() para verificar la integridad de los archivos
})

# Convierte bytes a un formato legible (KB, MB, GB)
formato_tamano <- function(bytes) {
  unidades <- c("B", "KB", "MB", "GB", "TB")
  i <- if (bytes > 0) min(floor(log(bytes, 1024)), length(unidades) - 1) else 0
  sprintf("%.2f %s", bytes / 1024^i, unidades[i + 1])
}

mensaje <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)


# ------------------------------------------------------------------------------
# 1. Rutas del proyecto (relativas a la raíz del proyecto de RStudio)
# ------------------------------------------------------------------------------
# here::here() encuentra la raíz del proyecto (.Rproj); si 'here' no está
# instalado, se usa el directorio de trabajo (en un proyecto RStudio es la raíz).
raiz <- if (requireNamespace("here", quietly = TRUE)) here::here() else getwd()

dir_raw   <- file.path(raiz, "data", "raw")
dir_ref   <- file.path(raiz, "data", "reference")
dir_tmp   <- file.path(dir_raw, "_descarga_geo")       # descarga del .tar
dir_extr  <- file.path(dir_tmp, "extraido")            # extracción temporal

for (d in c(dir_raw, dir_ref, dir_tmp)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
mensaje("Raíz del proyecto: ", raiz)


# ------------------------------------------------------------------------------
# 2. Metadatos de las muestras (verificados en GEO, 22 GSM)
# ------------------------------------------------------------------------------
# Condición en GEO: "Multiple Sclerosis" (EM) o
# "idiopathic intracranial hypertension" (IIH = control del paper).
GSE <- "GSE138266"

pacientes <- c(
  # LCR (CSF): GSM4104122 - GSM4104133
  "MS19270", "MS58637", "MS71658", "MS49131", "MS60249", "MS74594",
  "PST83775", "PTC32190", "PST95809", "PTC41540", "PTC45044", "PTC85037",
  # Sangre (PBMC): GSM4104134 - GSM4104143
  "MS19270", "MS71658", "MS49131", "MS60249", "MS74594",
  "PST83775", "PTC32190", "PST95809", "PTC41540", "PTC85037"
)

muestras <- data.frame(
  gsm       = paste0("GSM", 4104122:4104143),
  paciente  = pacientes,
  tejido    = rep(c("CSF", "blood"), times = c(12, 10)),
  condicion = ifelse(startsWith(pacientes, "MS"), "MS", "IIH"),
  stringsAsFactors = FALSE
)
muestras$carpeta <- with(muestras, paste(gsm, tejido, condicion, paciente, sep = "_"))

write.csv(muestras, file.path(dir_raw, "metadata_muestras.csv"), row.names = FALSE)
mensaje("Metadatos guardados en data/raw/metadata_muestras.csv")
print(table(Tejido = muestras$tejido, Condicion = muestras$condicion))


# ------------------------------------------------------------------------------
# 3. Descarga de los archivos suplementarios de GEO (GSE138266_RAW.tar, 252 MB)
# ------------------------------------------------------------------------------
tar_path       <- file.path(dir_tmp, GSE, paste0(GSE, "_RAW.tar"))
tamano_tar_geo <- 264058880   # bytes, según filelist.txt de GEO

if (file.exists(tar_path) && file.size(tar_path) == tamano_tar_geo) {
  mensaje("El .tar de GEO ya existe y está completo; se omite la descarga.")
} else {
  mensaje("Descargando ", GSE, "_RAW.tar desde GEO (~252 MB)...")
  tryCatch(
    getGEOSuppFiles(GSE, baseDir = dir_tmp, makeDirectory = TRUE,
                    filter_regex = "RAW\\.tar$"),
    error = function(e) stop("Falló la descarga desde GEO: ", conditionMessage(e),
                             "\nRevisa tu conexión y vuelve a correr el script.")
  )
  if (!file.exists(tar_path) || file.size(tar_path) != tamano_tar_geo) {
    stop("El .tar descargado está incompleto (", formato_tamano(file.size(tar_path)),
         "). Bórralo y vuelve a correr el script.")
  }
  mensaje("Descarga de GEO completa.")
}


# ------------------------------------------------------------------------------
# 4. Descomprimir y organizar: una subcarpeta por muestra
# ------------------------------------------------------------------------------
# GEO nombra los archivos como GSM..._<paciente>_<tejido>_GRCh38_<tipo>.tsv.gz.
# Los renombramos al formato que espera Seurat::Read10X():
#   barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz
# (genes.tsv.gz de Cell Ranger v2 -> features.tsv.gz; tiene 2 columnas
#  ID Ensembl + símbolo, compatible con Read10X(gene.column = 2)).
archivos_10x <- c(barcodes = "barcodes.tsv.gz",
                  genes    = "features.tsv.gz",
                  matrix   = "matrix.mtx.gz")

muestra_completa <- function(carpeta) {
  all(file.exists(file.path(dir_raw, carpeta, archivos_10x)))
}

pendientes <- muestras$carpeta[!vapply(muestras$carpeta, muestra_completa, logical(1))]

if (length(pendientes) == 0) {
  mensaje("Las 22 muestras ya están organizadas; se omite la extracción.")
} else {
  mensaje("Extrayendo el .tar (", length(pendientes), " muestras pendientes)...")
  dir.create(dir_extr, showWarnings = FALSE)
  untar(tar_path, exdir = dir_extr)
  extraidos <- list.files(dir_extr)
  
  for (i in which(muestras$carpeta %in% pendientes)) {
    gsm     <- muestras$gsm[i]
    destino <- file.path(dir_raw, muestras$carpeta[i])
    dir.create(destino, showWarnings = FALSE)
    
    for (tipo in names(archivos_10x)) {
      # Se busca por el GSM y no por el nombre del paciente, porque GEO tiene
      # inconsistencias (ej. GSM4104132: título PTC45044, archivos PST45044).
      origen <- grep(paste0("^", gsm, "_.*_", tipo, "\\.(tsv|mtx)\\.gz$"),
                     extraidos, value = TRUE)
      if (length(origen) != 1) {
        stop("No se encontró un archivo único '", tipo, "' para ", gsm)
      }
      file.copy(file.path(dir_extr, origen), file.path(destino, archivos_10x[[tipo]]),
                overwrite = TRUE)
    }
    mensaje("  OK -> ", muestras$carpeta[i])
  }
  # Se borran los archivos extraídos temporales (se conserva el .tar para no
  # tener que descargarlo de nuevo).
  unlink(dir_extr, recursive = TRUE)
}


# ------------------------------------------------------------------------------
# 5. Descarga del objeto integrado de Figshare (MS_CSF.h5ad, ~1 GB)
# ------------------------------------------------------------------------------
url_figshare <- "https://ndownloader.figshare.com/files/27405182"
h5ad_path    <- file.path(dir_ref, "MS_CSF.h5ad")
md5_figshare <- "6ca095a0b0e1b8a41ad380aef71f54b2"   # MD5 publicado por Figshare

h5ad_valido <- function() {
  file.exists(h5ad_path) && unname(md5sum(h5ad_path)) == md5_figshare
}

if (h5ad_valido()) {
  mensaje("MS_CSF.h5ad ya existe y su MD5 es correcto; se omite la descarga.")
} else {
  mensaje("Descargando MS_CSF.h5ad desde Figshare (~1 GB, puede tardar)...")
  tmp_h5ad <- paste0(h5ad_path, ".part")
  estado <- tryCatch(
    download.file(url_figshare, destfile = tmp_h5ad, mode = "wb", quiet = FALSE),
    error = function(e) { message("Error: ", conditionMessage(e)); 1 }
  )
  if (estado != 0) stop("Falló la descarga de Figshare. Vuelve a correr el script.")
  file.rename(tmp_h5ad, h5ad_path)
  if (!h5ad_valido()) {
    stop("El MD5 de MS_CSF.h5ad no coincide (archivo corrupto). Bórralo y reintenta.")
  }
  mensaje("Descarga de Figshare completa y verificada (MD5 OK).")
}


# ------------------------------------------------------------------------------
# 6. Resumen final
# ------------------------------------------------------------------------------
cat("\n==================== RESUMEN DE LA DESCARGA ====================\n")

ok <- vapply(muestras$carpeta, muestra_completa, logical(1))
cat(sprintf("Muestras organizadas: %d de %d\n", sum(ok), nrow(muestras)))
if (any(!ok)) cat("  INCOMPLETAS:", muestras$carpeta[!ok], sep = "\n  ")

cat("\nArchivos por carpeta (data/raw/):\n")
tamano_muestras <- 0
for (carpeta in muestras$carpeta) {
  archivos <- list.files(file.path(dir_raw, carpeta), full.names = TRUE)
  tam      <- sum(file.size(archivos))
  tamano_muestras <- tamano_muestras + tam
  cat(sprintf("  %-38s %s  (%s)\n", carpeta,
              paste(basename(archivos), collapse = ", "), formato_tamano(tam)))
}

tamano_tar  <- if (file.exists(tar_path))  file.size(tar_path)  else 0
tamano_h5ad <- if (file.exists(h5ad_path)) file.size(h5ad_path) else 0

cat("\nTamaños:\n")
cat("  Matrices organizadas (22 muestras):", formato_tamano(tamano_muestras), "\n")
cat("  GSE138266_RAW.tar (respaldo):      ", formato_tamano(tamano_tar), "\n")
cat("  MS_CSF.h5ad (Figshare):            ", formato_tamano(tamano_h5ad), "\n")
cat("  TOTAL en disco:                    ",
    formato_tamano(tamano_muestras + tamano_tar + tamano_h5ad), "\n")

cat("\nArchivo de referencia Figshare:",
    if (h5ad_valido()) "PRESENTE y verificado (MD5 OK)" else "FALTA o CORRUPTO", "\n")
cat("================================================================\n")
cat("Listo para el Paso 2 (QC): scripts/02_QC.R\n")