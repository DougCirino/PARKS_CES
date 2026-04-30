# =========================================================
# Bivariate (4x4) Supply x Demand from Shapefile
#  - Inputs: Shapefile with PMR, DMR, PAR, DAR
#  - Outputs: GeoTIFFs (4x4 classes), x11 biplots, map previews
#  - Palette: Orange↔Green (Supply=Green on X; Demand=Orange on Y)
#  - Exports: CSV (value,label,hex) for the 16 categories
#  - No QGIS style files; no mini-graph inset in plots/maps
#  - Extra: opens a separate x11 with just the 4×4 color scale
# =========================================================

suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(dplyr)
  library(grDevices)
})

# ------------------------------
# 0) User inputs
# ------------------------------
shp_path   <- "G:/Meu Drive/PESQUISA - Ecologia Urbana e Serviços Ecossistêmicos/Doutorado/Dados/CES_Parques/Grids/DADOS_CES/grid_data_CES.shp"
setwd("G:/Meu Drive/PESQUISA - Ecologia Urbana e Serviços Ecossistêmicos/Doutorado/Dados/CES_Parques/Grids/DADOS_CES")

target_epsg    <- 31983
res_m          <- 50
n_classes      <- 4
palette_name   <- "OrangeGrn"   # <— padrão agora é Orange↔Green
draw_map_frame <- FALSE
SHOW_PERCENT   <- FALSE         # (sem efeito nos plots atuais; mantido p/ compat.)

out_pmr_tif <- "BiClass_PMRxDMR_50m_4x4.tif"
out_par_tif <- "BiClass_PARxDAR_50m_4x4.tif"
out_palette_csv <- "bivar_4x4_palette.csv"

# ------------------------------
# 1) Read shapefile and ensure projected CRS
# ------------------------------
g <- st_read(shp_path, quiet = TRUE)

req_cols <- c("PMR","DMR","PAR","DAR")
miss <- setdiff(req_cols, names(g))
if (length(miss)) stop("Missing required columns: ", paste(miss, collapse = ", "))

# Use representative points if input is polygonal
if (any(st_geometry_type(g) %in% c("POLYGON","MULTIPOLYGON"))) {
  g_pts <- st_point_on_surface(g)
} else {
  g_pts <- g
}

if (is.na(st_crs(g_pts))) stop("Shapefile has no CRS; define it before running.")
if (sf::st_is_longlat(g_pts)) {
  message("Input in degrees. Transforming to EPSG:", target_epsg, " ...")
  g_pts <- st_transform(g_pts, target_epsg)
} else {
  crs_epsg <- tryCatch(st_crs(g_pts)$epsg, error = function(e) NA_integer_)
  if (!is.na(crs_epsg) && crs_epsg != target_epsg) {
    g_pts <- st_transform(g_pts, target_epsg)
  }
}

# ------------------------------
# 2) Helpers
# ------------------------------
range01 <- function(x) {
  if (all(is.na(x))) return(x)
  r <- range(x, na.rm = TRUE)
  if (diff(r) == 0) ifelse(is.na(x), NA_real_, 0) else (x - r[1])/(r[2]-r[1])
}

qcut <- function(x, n = 4) {
  x0 <- range01(x)
  qs <- quantile(x0, probs = seq(0,1,length.out = n+1), na.rm = TRUE, type = 7)
  if (length(unique(qs)) < (n+1)) {
    x0 <- jitter(x0, factor = 1e-7)
    qs <- quantile(x0, probs = seq(0,1,length.out = n+1), na.rm = TRUE, type = 7)
  }
  cut(x0, breaks = qs, include.lowest = TRUE, labels = FALSE)
}

# Palette: Supply→Green (X), Demand→Orange (Y); blend row-major
bivar_palette_scheme <- function(n = 4, scheme = "OrangeGrn") {
  low_col <- "#f0f0f0"
  if (scheme == "PinkGrn")    { high_x <- "#008000"; high_y <- "#d95f9d"
  } else if (scheme == "OrangeGrn") { high_x <- "#008000"; high_y <- "#FF5C00"
  } else                       { high_x <- "#008000"; high_y <- "#6c83b5" }
  ramp_x <- colorRampPalette(c(low_col, high_x))(n)
  ramp_y <- colorRampPalette(c(low_col, high_y))(n)
  pal <- matrix(NA_character_, nrow = n, ncol = n)
  for (iy in 1:n) for (ix in 1:n) {
    cx <- col2rgb(ramp_x[ix]); cy <- col2rgb(ramp_y[iy])
    mix <- (cx + cy) / 2
    pal[iy, ix] <- rgb(mix[1]/255, mix[2]/255, mix[3]/255)
  }
  pal
}
make_palette_vector <- function(pal_mat) as.vector(t(pal_mat))  # row-major flatten

bivar_class <- function(x, y, n = 4) {
  xb <- qcut(x, n); yb <- qcut(y, n)
  code  <- as.integer((yb - 1L) * n + xb)  # 1..n^2 (row-major): (y-1)*n + x
  label <- paste0(xb, "-", yb)
  list(code = code, label = label, x_n = range01(x), y_n = range01(y), xb = xb, yb = yb)
}

code_to_color <- function(code, n = 4, pal_mat = NULL) {
  if (is.null(pal_mat)) pal_mat <- bivar_palette_scheme(n, scheme = palette_name)
  x <- ((code - 1L) %% n) + 1L
  y <- ((code - 1L) %/% n) + 1L
  unname(pal_mat[cbind(y, x)])
}

# Mini-graph (legend-only) — called in a separate x11
draw_bivar_minigraph <- function(n = 4,
                                 xlab = "Supply (low → high)",
                                 ylab = "Demand (low → high)",
                                 pal_mat = NULL) {
  if (is.null(pal_mat)) pal_mat <- bivar_palette_scheme(n, scheme = palette_name)
  par(xpd = NA, bty = "n")
  plot.new()
  rect(0, 0, 1, 1, col = NA, border = NA)
  for (iy in 1:n) for (ix in 1:n) {
    x0 <- (ix-1)/n; y0 <- (iy-1)/n
    rect(x0, y0, x0+1/n, y0+1/n, col = pal_mat[iy, ix], border = "white", lwd = 0.6)
  }
  mtext(xlab, side = 1, line = 2, cex = 0.9)
  mtext(ylab, side = 2, line = 2, cex = 0.9)
}

# ------------------------------
# 3) Compute 4x4 classes at feature level
# ------------------------------
dat <- g_pts |> st_drop_geometry() |> transmute(PMR, DMR, PAR, DAR)

b_PMR <- bivar_class(dat$PMR, dat$DMR, n = n_classes)
b_PAR <- bivar_class(dat$PAR, dat$DAR, n = n_classes)

g_pts$BiCode_PMRxDMR  <- b_PMR$code
g_pts$BiLabel_PMRxDMR <- b_PMR$label
g_pts$BiCode_PARxDAR  <- b_PAR$code
g_pts$BiLabel_PARxDAR <- b_PAR$label

# ------------------------------
# 4) Base raster and rasterize codes (force 1..16)
# ------------------------------
bb   <- st_bbox(g_pts)
pad  <- res_m / 2
rref <- rast(xmin = bb["xmin"] - pad, xmax = bb["xmax"] + pad,
             ymin = bb["ymin"] - pad, ymax = bb["ymax"] + pad,
             resolution = res_m, crs = st_crs(g_pts)$wkt)

v <- vect(g_pts)
r_PMR <- rasterize(v, rref, field = "BiCode_PMRxDMR", fun = "first")
r_PAR <- rasterize(v, rref, field = "BiCode_PARxDAR",  fun = "first")

# Enforce discrete classes 1..16 exactly; anything else -> NA; store as factor
r_PMR <- round(r_PMR); r_PAR <- round(r_PAR)
r_PMR[r_PMR < 1 | r_PMR > n_classes^2] <- NA
r_PAR[r_PAR < 1 | r_PAR > n_classes^2] <- NA
r_PMR <- as.factor(r_PMR); r_PAR <- as.factor(r_PAR)

# Write GeoTIFFs (Byte, LZW)
#writeRaster(r_PMR, out_pmr_tif, overwrite = TRUE, datatype = "INT1U", gdal = "COMPRESS=LZW")
#writeRaster(r_PAR, out_par_tif, overwrite = TRUE, datatype = "INT1U", gdal = "COMPRESS=LZW")
#cat("Saved TIFFs:\n  -", normalizePath(out_pmr_tif), "\n  -", normalizePath(out_par_tif), "\n")

# ------------------------------
# 5) Export palette CSV (value, label, hex)
# ------------------------------
pal_mat  <- bivar_palette_scheme(n_classes, scheme = palette_name)
pal_vec  <- make_palette_vector(pal_mat)  # length n^2, row-major

labels16 <- as.vector(t(matrix(
  paste0(rep(1:n_classes, each = n_classes), "-", rep(1:n_classes, times = n_classes)),
  nrow = n_classes, byrow = TRUE
)))

palette_tbl <- data.frame(
  value = 1:(n_classes * n_classes),
  label = labels16,
  hex   = pal_vec,
  stringsAsFactors = FALSE
)
#write.csv(palette_tbl, out_palette_csv, row.names = FALSE)
#cat("Palette CSV:\n  -", normalizePath(out_palette_csv), "\n")

# ------------------------------
# 6) POINT BIPLOTS (no inset)
# ------------------------------
cols_PMR <- code_to_color(g_pts$BiCode_PMRxDMR, n = n_classes, pal_mat = pal_mat)
cols_PAR <- code_to_color(g_pts$BiCode_PARxDAR,  n = n_classes, pal_mat = pal_mat)

# PMR x DMR
x11(width = 7, height = 6)
par(mar = c(4,4,2,1), bty = "n")
plot(b_PMR$x_n, b_PMR$y_n,
     pch = 16, cex = 0.65, col = cols_PMR,
     xlab = "PMR (normalized 0–1)",
     ylab = "DMR (normalized 0–1)",
     main = sprintf("Bivariate %dx%d: PMR × DMR", n_classes, n_classes))
grid()

# PAR x DAR
x11(width = 7, height = 6)
par(mar = c(4,4,2,1), bty = "n")
plot(b_PAR$x_n, b_PAR$y_n,
     pch = 16, cex = 0.65, col = cols_PAR,
     xlab = "PAR (normalized 0–1)",
     ylab = "DAR (normalized 0–1)",
     main = sprintf("Bivariate %dx%d: PAR × DAR", n_classes, n_classes))
grid()

# ------------------------------
# 7) MAP previews (no inset; no frame)
# ------------------------------
plot_bivar_raster <- function(r, main = "", n = 4, pal_mat = NULL, frame = FALSE) {
  if (is.null(pal_mat)) pal_mat <- bivar_palette_scheme(n, scheme = palette_name)
  pal_vec <- make_palette_vector(pal_mat)
  x11(width = 7, height = 6)
  par(mar = c(2,2,3,6), bty = "n")
  plot(r, col = pal_vec, axes = FALSE, legend = FALSE, main = main, box = FALSE)
  if (frame) box()
}
plot_bivar_raster(r_PMR, main = sprintf("Raster bi-class (%dx%d): PMR × DMR", n_classes, n_classes),
                  n = n_classes, pal_mat = pal_mat, frame = draw_map_frame)
plot_bivar_raster(r_PAR, main = sprintf("Raster bi-class (%dx%d): PAR × DAR", n_classes, n_classes),
                  n = n_classes, pal_mat = pal_mat, frame = draw_map_frame)

# ------------------------------
# 8) Separate x11 with ONLY the 4×4 color scale (legend)
# ------------------------------
x11(width = 5, height = 5)
par(mar = c(3,3,1,1))
draw_bivar_minigraph(n = n_classes,
                     xlab = "Supply (low → high)",
                     ylab = "Demand (low → high)",
                     pal_mat = pal_mat)

cat("Done.\n")

