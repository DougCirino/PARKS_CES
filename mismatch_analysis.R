# =========================================================
# Spatial (mis)match analysis from the original shapefile
#   Inputs (columns in grid_data_CES.shp):
#     - PMR, PAR  = supply+flow (two CES)
#     - DMR, DAR  = demand       (two CES)
#     - POP       = normalized population per grid cell (0..1; 0 = no people)
#     - park_dist = distance to nearest park (not used here, but kept available)
#
#   A) "No autocorrelation between supply and demand?"
#      - Global *bivariate* Moran’s I for PMR~DMR and PAR~DAR (permutation test)
#      - Moran scatterplots (z(S) vs lag_W(z(D))) with linear fit
#
#   B) Continuous mismatch maps: M = z(Demand) – z(Supply)  (for PMR/DMR and PAR/DAR)
#      - x11 map for each CES + GeoTIFF export
#
#   C) Four-category map of positive mismatch (Dz − Sz > 0)
#      - 0=None positive (forced where POP==0 or M<=0)
#      - 1=PMR>0 only
#      - 2=PAR>0 only
#      - 3=Both >0
#      - x11 map + GeoTIFF export
#
#   Notes:
#     - All rasters built on a 50 m grid in EPSG:31983
#     - x11() opened for each plot (so you can save from the device)
# =========================================================

suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(spdep)
})

# ------------------------------
# 0) User inputs
# ------------------------------
shp_path <- "G:/Meu Drive/PESQUISA - Ecologia Urbana e Serviços Ecossistêmicos/Doutorado/Dados/CES_Parques/Grids/DADOS_CES/grid_data_CES.shp"
setwd("G:/Meu Drive/PESQUISA - Ecologia Urbana e Serviços Ecossistêmicos/Doutorado/Dados/CES_Parques/Grids/DADOS_CES")

target_epsg <- 31983
res_m       <- 50

# ------------------------------
# 1) Read shapefile (point-on-surface if polygons) & reproject to metric CRS
# ------------------------------
g <- st_read(shp_path, quiet = TRUE)

req_cols <- c("PMR","DMR","PAR","DAR","POP","park_dist")
miss <- setdiff(req_cols, names(g))
if (length(miss)) stop("Missing required columns in shapefile: ", paste(miss, collapse = ", "))

if (any(st_geometry_type(g) %in% c("POLYGON","MULTIPOLYGON"))) g <- st_point_on_surface(g)
if (is.na(st_crs(g))) stop("Shapefile has no CRS; define it before running.")
if (sf::st_is_longlat(g)) {
  message("Input is in degrees; transforming to EPSG:", target_epsg)
  g <- st_transform(g, target_epsg)
} else {
  crs_epsg <- tryCatch(st_crs(g)$epsg, error = function(e) NA_integer_)
  if (!is.na(crs_epsg) && crs_epsg != target_epsg) g <- st_transform(g, target_epsg)
}

# ------------------------------
# 2) Build 50 m raster template and rasterize variables
# ------------------------------
bb  <- st_bbox(g)
pad <- res_m/2
rref <- rast(xmin=bb["xmin"]-pad, xmax=bb["xmax"]+pad,
             ymin=bb["ymin"]-pad, ymax=bb["ymax"]+pad,
             resolution=res_m, crs=st_crs(g)$wkt)

v      <- vect(g)
r_S_PMR <- rasterize(v, rref, field="PMR", fun="mean")
r_D_PMR <- rasterize(v, rref, field="DMR", fun="mean")
r_S_PAR <- rasterize(v, rref, field="PAR", fun="mean")
r_D_PAR <- rasterize(v, rref, field="DAR", fun="mean")
r_POP   <- rasterize(v, rref, field="POP", fun="mean")  # normalized pop already

# Reference mask: anywhere we have any CES value
ref_mask <- (!is.na(r_S_PMR)) | (!is.na(r_D_PMR)) | (!is.na(r_S_PAR)) | (!is.na(r_D_PAR))
r_POP <- mask(r_POP, ref_mask)

# ------------------------------
# 3) Helpers: z-score & neighbors from grid
# ------------------------------
zfun <- function(x) (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
z_stack <- function(S, D) {
  Sz <- app(S, zfun); Dz <- app(D, zfun); c(Sz=Sz, Dz=Dz)
}
pmr_z <- z_stack(r_S_PMR, r_D_PMR)
par_z <- z_stack(r_S_PAR, r_D_PAR)

# Queen neighbors for the grid (derived from rref dimensions)
d  <- dim(rref)  # c(nrow, ncol, nlyr)
nb_all <- spdep::cell2nb(d[1], d[2], type="queen")

# Robust listw + spatial lag (avoid class/attr pitfalls on Windows)
safe_listw <- function(nb) {
  lw <- spdep::nb2listw(nb, style="W")
  lw$weights <- lapply(lw$weights, function(w) as.numeric(w))
  lw$style <- "W"
  lw
}
lagW <- function(lw, x) {
  nb <- lw$neighbours; ww <- lw$weights
  x  <- as.numeric(x)
  out <- numeric(length(nb))
  for (i in seq_along(nb)) {
    if (length(nb[[i]]) == 0L) out[i] <- 0 else out[i] <- sum(x[nb[[i]]] * ww[[i]])
  }
  out
}
.prepare_for_moran <- function(r_layer, nb_all) {
  v    <- as.vector(r_layer[])
  keep <- is.finite(v)
  nbk  <- spdep::subset.nb(nb_all, keep)
  has  <- lengths(nbk) > 0
  nb   <- spdep::subset.nb(nbk, has)
  idx  <- which(keep)[has]
  x    <- as.numeric(v[idx])
  list(x=x, nb=nb, idx=idx)
}

# ------------------------------
# 4) A) Global bivariate Moran’s I + Moran scatter (PMR~DMR; PAR~DAR)
#     - If I ≈ 0 and p>0.05 ⇒ no spatial autocorrelation between S and D
# ------------------------------
global_bv_moran <- function(Sz, Dz, nb_all, label, nsim=999) {
  # restrict to cells where both Sz & Dz are finite
  vS <- as.vector(Sz[]); vD <- as.vector(Dz[])
  keep <- is.finite(vS) & is.finite(vD)
  nbk  <- spdep::subset.nb(nb_all, keep)
  has  <- lengths(nbk) > 0
  nb   <- spdep::subset.nb(nbk, has)
  idx  <- which(keep)[has]
  zS   <- as.numeric(vS[idx]); zD <- as.numeric(vD[idx])
  lw   <- safe_listw(nb)
  lag_zD <- lagW(lw, zD)
  
  # local cross-product and global mean (bivariate Moran’s I)
  Ii    <- zS * lag_zD
  Iglob <- mean(Ii)
  
  # permutation test (shuffle zD)
  seeds <- sample.int(1e9, nsim)
  perm_fun <- function(seed) {
    set.seed(seed)
    zDp <- sample(zD, replace=FALSE)
    mean(zS * lagW(lw, zDp))
  }
  Iperm <- vapply(seeds, perm_fun, numeric(1))
  pval  <- (sum(abs(Iperm) >= abs(Iglob)) + 1) / (nsim + 1)
  
  # Moran scatter (zS vs lag(zD))
  x11(width=7, height=6)
  gg <- ggplot(data.frame(zS=zS, lag_zD=lag_zD), aes(zS, lag_zD)) +
    geom_point(alpha=.35, size=.6, color="grey30") +
    geom_smooth(method="lm", se=FALSE, color="black") +
    geom_vline(xintercept=0, linetype=2) +
    geom_hline(yintercept=0, linetype=2) +
    labs(title=paste0(label, " — Moran scatter (bivariate)"),
         x="Supply (z)", y="Spatial lag of Demand (z)") +
    theme_minimal()
  print(gg)
  
  message(sprintf("%s — Global bivariate Moran's I: %.3f (p=%.3f)", label, Iglob, pval))
  invisible(list(I=Iglob, p=pval))
}

cat("\n=== A) Bivariate Moran’s I (Supply vs Demand) ===\n")
pmr_bv <- global_bv_moran(pmr_z[["Sz"]], pmr_z[["Dz"]], nb_all, "PMR ~ DMR", nsim=999)
par_bv <- global_bv_moran(par_z[["Sz"]], par_z[["Dz"]], nb_all, "PAR ~ DAR", nsim=999)

# Interpretation quick guide:
#  - If I ≈ 0 and p > 0.05: no spatial autocorrelation between supply and demand ⇒ supports "local mismatch".
#  - If I > 0 and p <= 0.05: positive spatial coupling (S high where neighboring D high).
#  - If I < 0 and p <= 0.05: negative association.

# ------------------------------
# 5) B) Continuous mismatch maps: M = z(Demand) – z(Supply)
# ------------------------------
make_mismatch <- function(Sz, Dz) {
  M <- Dz - Sz
  names(M) <- "Mismatch"
  M
}
pmr_M <- make_mismatch(pmr_z[["Sz"]], pmr_z[["Dz"]])
par_M <- make_mismatch(par_z[["Sz"]], par_z[["Dz"]])

# pretty, symmetric color clip at 99th percentile of |M|
sym_clip <- function(rM, q=0.99) as.numeric(quantile(abs(values(rM)), q, na.rm=TRUE))

plot_mismatch <- function(rM, title) {
  lim <- sym_clip(rM, 0.99)
  df  <- as.data.frame(rM, xy=TRUE, na.rm=TRUE)
  names(df) <- c("x","y","M")
  df$M <- pmin(pmax(df$M, -lim), lim)
  
  x11(width=7, height=6)
  gg <- ggplot(df, aes(x, y, fill=M)) +
    geom_raster() + coord_equal() +
    scale_fill_gradient2(low="#2F8F2F", mid="white", high="#FF8A33",
                         limits=c(-lim, lim), oob=squish,
                         breaks=pretty(c(-lim, lim), 5),
                         labels=label_number(accuracy=0.1)) +
    labs(title=title, fill="Mismatch\n(Dz − Sz)") +
    theme_void() + theme(legend.position="right")
  print(gg)
}

cat("\n=== B) Continuous mismatch maps (Dz − Sz) ===\n")
plot_mismatch(pmr_M, "PMR — Mismatch map (Dz − Sz)")
plot_mismatch(par_M, "PAR — Mismatch map (Dz − Sz)")

# GeoTIFF exports
writeRaster(pmr_M, "PMR_mismatch_zD_minus_zS_50m.tif",
            overwrite=TRUE, datatype="FLT4S", gdal="COMPRESS=LZW")
writeRaster(par_M, "PAR_mismatch_zD_minus_zS_50m.tif",
            overwrite=TRUE, datatype="FLT4S", gdal="COMPRESS=LZW")

# ------------------------------
# 6) C) Four-category map of positive mismatch, forcing POP==0 to "None"
#     Codes:
#       0 = None positive (POP==0 OR M<=0 for both CES)
#       1 = PMR>0 only
#       2 = PAR>0 only
#       3 = Both >0
# ------------------------------
pop_has_people <- (!is.na(r_POP)) & (r_POP > 0)

pmr_pos <- (pmr_M > 0) & pop_has_people
par_pos <- (par_M > 0) & pop_has_people

combo <- (pmr_pos * 1) + (par_pos * 2)  # 0,1,2,3 with pop==0 already 0 via mask above
combo <- classify(combo, rcl = matrix(c(0,0, 1,1, 2,2, 3,3), ncol=2, byrow=TRUE))
combo <- as.factor(combo)
levels(combo) <- data.frame(
  ID    = 0:3,
  label = c("None positive","PMR>0 only","PAR>0 only","Both >0")
)

# x11 map
x11(width=7, height=6)
plot(combo,
     main="Cells with positive mismatch (Dz − Sz > 0)\n(population==0 forced to None)",
     col=c("grey90", "#FFA866", "#66B2FF", "#FF4D4D"),
     plg=list(title="Positive mismatch"))

# Export GeoTIFF
writeRaster(combo, "CES_positive_mismatch_4cats_50m.tif",
            overwrite=TRUE, datatype="INT1U", gdal="COMPRESS=LZW")

# Quick area shares (by grid cell count)
cat("\nArea shares by class (grid-cell based):\n")
print( as.data.frame(freq(combo, useNA="no")) |>
         transform(share = count / sum(count)) )


# =========================================================
# D) Same-cell association: PMR×DMR and PAR×DAR
#     - Pearson & Spearman correlations (cell-wise)
#     - Hex-scatter plots with linear fit (x11 windows)
#     - Uses existing rasters: r_S_PMR, r_D_PMR, r_S_PAR, r_D_PAR
# =========================================================

suppressPackageStartupMessages({ library(terra); library(ggplot2); library(broom) })

# Helper: build paired df from two rasters on the same grid
make_pair_df <- function(Sr, Dr, s_name = "Supply", d_name = "Demand") {
  df <- as.data.frame(c(Sr, Dr), na.rm = TRUE)
  names(df) <- c(s_name, d_name)
  # keep only finite pairs
  df <- df[is.finite(df[[1]]) & is.finite(df[[2]]), , drop = FALSE]
  df
}

# Helper: correlation tests + neat print
corr_summary <- function(df, s_name, d_name, label) {
  x <- df[[s_name]]; y <- df[[d_name]]
  pear <- suppressWarnings(cor.test(x, y, method = "pearson"))
  spear <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  cat(sprintf("\n[%s] Same-cell association\n", label))
  cat(sprintf("  Pearson:  r = %.3f,  95%% CI [%.3f, %.3f],  p = %.3g,  n = %d\n",
              pear$estimate, pear$conf.int[1], pear$conf.int[2], pear$p.value, length(x)))
  cat(sprintf("  Spearman:  ρ = %.3f,              p = %.3g,  n = %d\n",
              spear$estimate, spear$p.value, length(x)))
  # return a tidy list in case you want to save
  list(
    pearson  = broom::tidy(pear),
    spearman = broom::tidy(spear),
    n = length(x)
  )
}

# Helper: hex-scatter with linear fit
plot_hex_scatter <- function(df, s_name, d_name, title) {
  x11(width = 7, height = 6)
  gg <- ggplot(df, aes_string(s_name, d_name)) +
    geom_bin2d(bins = 60) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.9, color = "black") +
    scale_fill_continuous(type = "viridis") +
    labs(title = title, x = s_name, y = d_name, fill = "Count") +
    theme_minimal()
  print(gg)
}

# ---------------- PMR × DMR ----------------
pmr_df <- make_pair_df(r_S_PMR, r_D_PMR, s_name = "PMR", d_name = "DMR")
pmr_corr <- corr_summary(pmr_df, "PMR", "DMR", "PMR × DMR")
plot_hex_scatter(pmr_df, "PMR", "DMR", "Same-cell association: PMR × DMR")

# ---------------- PAR × DAR ----------------
par_df <- make_pair_df(r_S_PAR, r_D_PAR, s_name = "PAR", d_name = "DAR")
par_corr <- corr_summary(par_df, "PAR", "DAR", "PAR × DAR")
plot_hex_scatter(par_df, "PAR", "DAR", "Same-cell association: PAR × DAR")

# (Optional) Save correlation tables
# write.csv(pmr_corr$pearson,  "PMR_DMR_Pearson.csv",  row.names = FALSE)
# write.csv(pmr_corr$spearman, "PMR_DMR_Spearman.csv", row.names = FALSE)
# write.csv(par_corr$pearson,  "PAR_DAR_Pearson.csv",  row.names = FALSE)
# write.csv(par_corr$spearman, "PAR_DAR_Spearman.csv", row.names = FALSE)


# ---- Dependencies ----
suppressPackageStartupMessages({ library(ggplot2) })

# ---- Helper for nicer p-values ----
p_fmt <- function(p, digits = 3, eps = 1e-16) {
  ifelse(p < eps, paste0("< ", format(eps, scientific = TRUE)),
         format.pval(p, digits = digits, eps = eps))
}

# ---- Build paired vectors from rasters (same-cell) ----
pair_from_rasters <- function(Sr, Dr, n_sample = 100000L, seed = 42) {
  S <- as.numeric(terra::values(Sr))
  D <- as.numeric(terra::values(Dr))
  keep <- is.finite(S) & is.finite(D)
  S <- S[keep]; D <- D[keep]
  # optional: z-score both to stabilize
  S <- as.numeric(scale(S)); D <- as.numeric(scale(D))
  set.seed(seed)
  if (length(S) > n_sample) {
    idx <- sample.int(length(S), n_sample)
    S <- S[idx]; D <- D[idx]
  }
  data.frame(S = S, D = D)
}

# ---- Compute correlations with CI for Pearson ----
corr_summary <- function(df) {
  # Pearson (with CI)
  ct_p <- cor.test(df$S, df$D, method = "pearson", alternative = "two.sided")
  # Spearman (no CI by default in base R)
  ct_s <- cor.test(df$S, df$D, method = "spearman", alternative = "two.sided", exact = FALSE)
  list(
    pearson = list(r = unname(ct_p$estimate),
                   ci = unname(ct_p$conf.int),
                   p  = ct_p$p.value,
                   n  = length(df$S)),
    spearman = list(rho = unname(ct_s$estimate),
                    p   = ct_s$p.value,
                    n   = length(df$S))
  )
}

# ---- Make a clean scatter (points, not pixels) ----
scatter_plot <- function(df, title, subtitle = NULL) {
  x11(width = 7, height = 6)
  gg <- ggplot(df, aes(S, D)) +
    geom_point(alpha = 0.1, size = 0.4) +   # points with transparency
    geom_smooth(method = "lm", se = FALSE) +
    labs(title = title, subtitle = subtitle,
         x = "Supply+Flow (z)", y = "Demand (z)") +
    theme_minimal()
  print(gg)
}

# ==========================
# PMR × DMR
# ==========================
df_pmr <- pair_from_rasters(r_S_PMR, r_D_PMR, n_sample = 100000)
cs_pmr <- corr_summary(df_pmr)

pmr_sub <- sprintf(
  "Pearson r = %.3f (95%% CI [%.3f, %.3f], p %s); Spearman ρ = %.3f (p %s); n = %s",
  cs_pmr$pearson$r, cs_pmr$pearson$ci[1], cs_pmr$pearson$ci[2], p_fmt(cs_pmr$pearson$p, digits=4, eps=1e-16),
  cs_pmr$spearman$rho, p_fmt(cs_pmr$spearman$p, digits=4, eps=1e-16),
  format(cs_pmr$pearson$n, big.mark = ",")
)
scatter_plot(df_pmr, "Same-cell association: PMR × DMR", pmr_sub)

# ==========================
# PAR × DAR
# ==========================
df_par <- pair_from_rasters(r_S_PAR, r_D_PAR, n_sample = 100000)
cs_par <- corr_summary(df_par)

par_sub <- sprintf(
  "Pearson r = %.3f (95%% CI [%.3f, %.3f], p %s); Spearman ρ = %.3f (p %s); n = %s",
  cs_par$pearson$r, cs_par$pearson$ci[1], cs_par$pearson$ci[2], p_fmt(cs_par$pearson$p, digits=4, eps=1e-16),
  cs_par$spearman$rho, p_fmt(cs_par$spearman$p, digits=4, eps=1e-16),
  format(cs_par$pearson$n, big.mark = ",")
)
scatter_plot(df_par, "Same-cell association: PAR × DAR", par_sub)

