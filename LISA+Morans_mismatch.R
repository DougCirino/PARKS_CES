# =========================================================
# LISA + Moran (POP>0) with parallel permutations (Windows-safe)
#   - Reuses existing rasters (r_S_*, r_D_*, r_POP, rref, nb_all, pmr_z/par_z)
#     or rebuilds them from grid_data_CES.shp if missing.
#   - Opens x11 windows for:
#       * Bivariate LISA map (POP>0)
#       * Moran scatter colored by LISA quadrants
#   - Key fix: export lagW/zS/zD/lw to PSOCK workers.
# =========================================================

suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(spdep)
  library(ggplot2)
})

# ---------- user inputs (only used if we need to rebuild) ----------
shp_path   <- "G:/Meu Drive/PESQUISA - Ecologia Urbana e Serviços Ecossistêmicos/Doutorado/Dados/CES_Parques/Grids/DADOS_CES/grid_data_CES.shp"
target_epsg <- 31983
res_m       <- 50

# ---------- 0) Build rasters / neighbors if not in memory ----------
need_build <- !all(c("r_S_PMR","r_D_PMR","r_S_PAR","r_D_PAR","r_POP","rref","nb_all") %in% ls())

if (need_build) {
  g <- sf::st_read(shp_path, quiet = TRUE)
  req_cols <- c("PMR","DMR","PAR","DAR","POP")
  miss <- setdiff(req_cols, names(g)); if (length(miss)) stop("Missing: ", paste(miss, collapse=", "))
  if (any(st_geometry_type(g) %in% c("POLYGON","MULTIPOLYGON"))) g <- st_point_on_surface(g)
  if (is.na(st_crs(g))) stop("Shapefile has no CRS.")
  if (sf::st_is_longlat(g)) g <- st_transform(g, target_epsg)
  
  bb  <- st_bbox(g); pad <- res_m/2
  rref <- rast(xmin=bb["xmin"]-pad, xmax=bb["xmax"]+pad,
               ymin=bb["ymin"]-pad, ymax=bb["ymax"]+pad,
               resolution=res_m, crs=st_crs(g)$wkt)
  
  v       <- vect(g)
  r_S_PMR <- rasterize(v, rref, field="PMR", fun="mean")
  r_D_PMR <- rasterize(v, rref, field="DMR", fun="mean")
  r_S_PAR <- rasterize(v, rref, field="PAR", fun="mean")
  r_D_PAR <- rasterize(v, rref, field="DAR", fun="mean")
  r_POP   <- rasterize(v, rref, field="POP", fun="mean")
  
  # neighbors for the full grid (queen)
  d <- dim(rref) # c(nrow, ncol, nlyr)
  nb_all <- spdep::cell2nb(d[1], d[2], type="queen")
}

# ---------- 1) Z-score stacks if not present ----------
zfun   <- function(x) (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
z_stack <- function(S, D) { Sz <- app(S, zfun); Dz <- app(D, zfun); c(Sz=Sz, Dz=Dz) }

if (!exists("pmr_z")) pmr_z <- z_stack(r_S_PMR, r_D_PMR)
if (!exists("par_z")) par_z <- z_stack(r_S_PAR, r_D_PAR)

# ---------- 2) Listw and spatial lag helpers (robust on Windows) ----------
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

# ---------- 3) PARALLEL helper (REPLACEMENT) ----------
# Exports named objects to workers before par*apply; robust on PSOCK.
permute_parallel <- function(nsim, FUN, export = character(), cores = max(1, parallel::detectCores()-1)) {
  if (cores <= 1) return(sapply(seq_len(nsim), function(i) FUN()))
  cl <- parallel::makeCluster(cores)
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  if (length(export)) parallel::clusterExport(cl, varlist = export, envir = parent.frame())
  parallel::parSapplyLB(cl, seq_len(nsim), function(i) FUN())
}

# ---------- 4) RUN LISA (REPLACEMENT) ----------
# POP>0 only; returns global I, local p, raster of quadrants; draws x11 plots.
n_cores   <- max(1, parallel::detectCores()-1)
nsim_perm <- 999

run_bv_lisa <- function(Sz, Dz, label, nsim = nsim_perm){
  # Keep only cells where both z-layers are finite AND POP > 0
  vS <- as.vector(Sz[]); vD <- as.vector(Dz[]); vP <- as.vector(r_POP[])
  keep <- is.finite(vS) & is.finite(vD) & is.finite(vP) & (vP > 0)
  
  # Neighbors restricted to kept cells; drop zero-neighbor cells
  nbk  <- spdep::subset.nb(nb_all, keep)
  has  <- lengths(nbk) > 0
  nb   <- spdep::subset.nb(nbk, has)
  idx  <- which(keep)[has]
  
  zS <- as.numeric(vS[idx]); zD <- as.numeric(vD[idx])
  lw <- safe_listw(nb)
  
  # Local bivariate Ii and global I
  lag_zD <- lagW(lw, zD)
  Ii     <- zS * lag_zD
  Iglob  <- mean(Ii)
  
  # ----- Global p-value via permutations (shuffle zD), parallel -----
  perm_fun <- function(){
    zDp <- sample(zD, replace = FALSE)
    mean(zS * lagW(lw, zDp))
  }
  Iperm <- permute_parallel(
    nsim,
    perm_fun,
    export = c("zS","zD","lw","lagW"),
    cores  = n_cores
  )
  p_global <- (sum(abs(Iperm) >= abs(Iglob)) + 1) / (nsim + 1)
  
  # ----- Local p-values (same Ii definition), parallel -----
  perm_local_fun <- function(){
    zDp <- sample(zD, replace = FALSE)
    zS * lagW(lw, zDp)
  }
  Ii_perm_mat <- permute_parallel(
    nsim,
    perm_local_fun,
    export = c("zS","zD","lw","lagW"),
    cores  = n_cores
  )
  if (is.null(dim(Ii_perm_mat))) Ii_perm_mat <- matrix(Ii_perm_mat, nrow = length(Ii))
  if (nrow(Ii_perm_mat) != length(Ii))  Ii_perm_mat <- matrix(Ii_perm_mat, nrow = length(Ii))
  
  ge_mat  <- abs(Ii_perm_mat) >= abs(Ii)
  p_local <- (rowSums(ge_mat) + 1) / (nsim + 1)
  
  # ----- Quadrants (p<=0.05) -----
  quad <- ifelse(zS >= 0 & lag_zD >= 0 & p_local <= 0.05, "HH",
                 ifelse(zS <  0 & lag_zD <  0 & p_local <= 0.05, "LL",
                        ifelse(zS <  0 & lag_zD >= 0 & p_local <= 0.05, "HL_deficit",
                               ifelse(zS >= 0 & lag_zD <  0 & p_local <= 0.05, "LH_surplus", "NS"))))
  lab_levels <- c("HL_deficit","LH_surplus","HH","LL","NS")
  
  # Map quadrants back to raster space (NA outside keep/has)
  quad_full <- rep(NA_character_, length(keep))
  quad_full[idx] <- quad
  r_quad <- Sz * NA
  terra::values(r_quad) <- match(quad_full, lab_levels)
  r_quad <- as.factor(r_quad)
  levels(r_quad) <- data.frame(ID = 1:length(lab_levels), label = lab_levels)
  
  # ----- Plot 1: LISA map (x11) -----
  x11(width=7, height=6)
  plot(r_quad, main = paste0(label, " — Bivariate LISA (POP>0)"),
       col = c("#FF8A33","#7BC67B","#2F8F2F","#C9A56A","grey85"),
       plg = list(title = "Cluster"))
  
  # ----- Plot 2: Moran scatter (colored by quadrants) -----
  x11(width=7, height=6)
  df <- data.frame(zS=zS, lag_zD=lag_zD, quad=factor(quad, levels=lab_levels))
  gg <- ggplot(df, aes(zS, lag_zD, color = quad)) +
    geom_point(alpha=.35, size=.6) +
    geom_smooth(method="lm", se=FALSE, color="black") +
    geom_vline(xintercept=0, linetype=2) +
    geom_hline(yintercept=0, linetype=2) +
    scale_color_manual(values = c(
      "HL_deficit"="#FF8A33","LH_surplus"="#7BC67B","HH"="#2F8F2F","LL"="#C9A56A","NS"="grey60"
    ), name="") +
    labs(title=paste0(label, " — Moran scatter (POP>0)"),
         x="Supply+flow (z)", y="Spatial lag of Demand (z)") +
    theme_minimal()
  print(gg)
  
  message(sprintf("%s — Global bivariate Moran’s I (POP>0): %.3f (p=%.3f)", label, Iglob, p_global))
  invisible(list(global=list(I=Iglob, p=p_global),
                 local=data.frame(Ii=Ii, p=p_local, quad=quad),
                 quad_r = r_quad))
}

# ---------- 5) RUN for both CES ----------
cat("\n=== A) Bivariate Moran’s I + LISA (POP>0 only) ===\n")
pmr_lisa <- run_bv_lisa(pmr_z[["Sz"]], pmr_z[["Dz"]], "PMR ~ DMR", nsim = nsim_perm)
par_lisa <- run_bv_lisa(par_z[["Sz"]], par_z[["Dz"]], "PAR ~ DAR", nsim = nsim_perm)

################
# =========================================================
# B) Continuous mismatch maps: M = z(Demand) − z(Supply)
#    (masked to POP>0, coherent with LISA)
# =========================================================
suppressPackageStartupMessages({ library(terra); library(ggplot2); library(scales) })

# Mask of exposed population
pop_has_people <- (!is.na(r_POP)) & (r_POP > 0)

# M = Dz - Sz (compute on full grid, then mask to POP>0 for mapping/exports)
make_mismatch <- function(Sz, Dz) { out <- Dz - Sz; names(out) <- "Mismatch"; out }
pmr_M_full <- make_mismatch(pmr_z[["Sz"]], pmr_z[["Dz"]])
par_M_full <- make_mismatch(par_z[["Sz"]], par_z[["Dz"]])

pmr_M <- mask(pmr_M_full, pop_has_people)
par_M <- mask(par_M_full, pop_has_people)

# symmetric clip for colors
.sym_clip <- function(rM, q=0.99) as.numeric(quantile(abs(values(rM)), q, na.rm = TRUE))

.plot_mismatch <- function(rM, title) {
  lim <- .sym_clip(rM, 0.99)
  df  <- as.data.frame(rM, xy = TRUE, na.rm = TRUE)
  names(df) <- c("x","y","M")
  df$M <- pmin(pmax(df$M, -lim), lim)
  x11(width = 7, height = 6)
  gg <- ggplot(df, aes(x, y, fill = M)) +
    geom_raster() + coord_equal() +
    scale_fill_gradient2(low = "#2F8F2F", mid = "white", high = "#FF8A33",
                         limits = c(-lim, lim), oob = squish,
                         breaks = pretty(c(-lim, lim), 5),
                         labels = label_number(accuracy = 0.1)) +
    labs(title = title, fill = "Mismatch\n(Dz − Sz)") +
    theme_void() + theme(legend.position = "right")
  print(gg)
}

# Draw & export (POP>0 masked)
.plot_mismatch(pmr_M, "PMR — Mismatch (Dz − Sz) — POP>0")
.plot_mismatch(par_M, "PAR — Mismatch (Dz − Sz) — POP>0")

writeRaster(pmr_M, "PMR_mismatch_zD_minus_zS_50m_POPgt0.tif",
            overwrite = TRUE, datatype = "FLT4S", gdal = "COMPRESS=LZW")
writeRaster(par_M, "PAR_mismatch_zD_minus_zS_50m_POPgt0.tif",
            overwrite = TRUE, datatype = "FLT4S", gdal = "COMPRESS=LZW")


# =========================================================
# C) Four-category positive mismatch (Dz − Sz > 0),
#    POP==0 forced to "None positive"
#    + exports of binary GeoTIFFs for PMR+, PAR+, BOTH+
# =========================================================

# Positive mismatch only where POP>0
pmr_pos <- (pmr_M_full > 0) & pop_has_people
par_pos <- (par_M_full > 0) & pop_has_people
both_pos <- pmr_pos & par_pos

# 4-class combo: 0=None, 1=PMR>0 only, 2=PAR>0 only, 3=Both>0
combo <- (pmr_pos * 1) + (par_pos * 2)           # logical*integer -> 0/1/2/3
combo <- classify(combo, rcl = matrix(c(0,0, 1,1, 2,2, 3,3), ncol = 2, byrow = TRUE))
combo <- as.factor(combo)
levels(combo) <- data.frame(
  ID    = 0:3,
  label = c("None positive","PMR>0 only","PAR>0 only","Both >0")
)

# Plot (POP==0 shows as "None positive")
x11(width = 7, height = 6)
plot(combo,
     main = "Cells with positive mismatch (Dz − Sz > 0) — POP==0 forced to None",
     col  = c("grey90", "#FFA866", "#66B2FF", "#FF4D4D"),
     plg  = list(title = "Positive mismatch"))

# Exports
writeRaster(combo,    "CES_positive_mismatch_4cats_50m_POPgt0.tif",
            overwrite = TRUE, datatype = "INT1U", gdal = "COMPRESS=LZW")
writeRaster(pmr_pos,  "PMR_positive_mismatch_50m_POPgt0.tif",
            overwrite = TRUE, datatype = "INT1U", gdal = "COMPRESS=LZW")
writeRaster(par_pos,  "PAR_positive_mismatch_50m_POPgt0.tif",
            overwrite = TRUE, datatype = "INT1U", gdal = "COMPRESS=LZW")
writeRaster(both_pos, "BOTH_positive_mismatch_50m_POPgt0.tif",
            overwrite = TRUE, datatype = "INT1U", gdal = "COMPRESS=LZW")

# Robust area shares (no terra::freq useNA arg)
vals <- as.vector(values(combo))
vals <- vals[!is.na(vals)]
tab  <- as.data.frame(table(vals), stringsAsFactors = FALSE)
names(tab) <- c("class_id","count")
tab$class_id <- as.integer(as.character(tab$class_id))

lev  <- levels(combo)[[1]]              # ID, label
area_tbl <- merge(tab, lev, by.x = "class_id", by.y = "ID", all.x = TRUE)
area_tbl$share <- area_tbl$count / sum(area_tbl$count)

cat("\nArea shares by class (grid-cell based; POP>0 logic applied):\n")
print(area_tbl[order(area_tbl$class_id), c("label","count","share")], row.names = FALSE)


# =========================================================
# D) Same-cell association (POP>0 only):
#    Pearson & Spearman + point scatters
# =========================================================
suppressPackageStartupMessages({ library(ggplot2); library(broom) })

# Mask originals to POP>0
mask_pop <- function(r) mask(r, pop_has_people)
PMR <- mask_pop(r_S_PMR); DMR <- mask_pop(r_D_PMR)
PAR <- mask_pop(r_S_PAR); DAR <- mask_pop(r_D_PAR)

.mk_df <- function(Sr, Dr, sx = "Supply", dy = "Demand") {
  df <- as.data.frame(c(Sr, Dr), na.rm = TRUE); names(df) <- c(sx, dy)
  df[is.finite(df[[1]]) & is.finite(df[[2]]), , drop = FALSE]
}
pmr_df <- .mk_df(PMR, DMR, "PMR", "DMR")
par_df <- .mk_df(PAR, DAR, "PAR", "DAR")

.corr_summary <- function(df, x, y, label, digits = 6) {
  pear  <- suppressWarnings(cor.test(df[[x]], df[[y]], method = "pearson"))
  spear <- suppressWarnings(cor.test(df[[x]], df[[y]], method = "spearman", exact = FALSE))
  cat(sprintf("\n[%s] Same-cell association (POP>0)\n", label))
  cat(sprintf("  Pearson:  r = %.3f,  95%% CI [%.3f, %.3f],  p = %.*g,  n = %d\n",
              pear$estimate, pear$conf.int[1], pear$conf.int[2], digits, pear$p.value, nrow(df)))
  cat(sprintf("  Spearman:  ρ = %.3f,                p = %.*g,  n = %d\n",
              spear$estimate, digits, spear$p.value, nrow(df)))
  list(pearson = tidy(pear), spearman = tidy(spear), n = nrow(df))
}

.plot_scatter <- function(df, x, y, title) {
  x11(width = 7, height = 6)
  gg <- ggplot(df, aes_string(x, y)) +
    geom_point(alpha = .2, size = .4) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.9, color = "black") +
    labs(title = title, x = x, y = y) +
    theme_minimal()
  print(gg)
}

pmr_corr <- .corr_summary(pmr_df, "PMR", "DMR", "PMR × DMR (cell-wise, POP>0)")
par_corr <- .corr_summary(par_df, "PAR", "DAR", "PAR × DAR (cell-wise, POP>0)")

.plot_scatter(pmr_df, "PMR", "DMR", "PMR × DMR — same-cell association (POP>0)")
.plot_scatter(par_df, "PAR", "DAR", "PAR × DAR — same-cell association (POP>0)")

###################

# =========================================================
# E) Actionable mismatch (LISA-HL) — POP>0 only
#    - HL = "High demand – Low supply" (significant)
#    - Build binary HL masks for PMR~DMR and PAR~DAR
#    - 4-class combo: 0=None, 1=PMR-HL only, 2=PAR-HL only, 3=Both HL
#    - x11 plot + GeoTIFF exports
#    Requires: pmr_lisa$quad_r, par_lisa$quad_r, pop_has_people
# =========================================================
suppressPackageStartupMessages({ library(terra) })

# ---- 1) Find the level ID that corresponds to "HL_deficit" in each LISA map
lev_pmr <- levels(pmr_lisa$quad_r)[[1]]  # columns: ID, label
lev_par <- levels(par_lisa$quad_r)[[1]]

hl_id_pmr <- lev_pmr$ID[lev_pmr$label == "HL_deficit"]
hl_id_par <- lev_par$ID[lev_par$label == "HL_deficit"]

if (length(hl_id_pmr) != 1 || length(hl_id_par) != 1) {
  stop("Could not identify unique 'HL_deficit' level in one of the LISA rasters.")
}

# ---- 2) Binary HL masks (significant HL cells), then restrict to POP>0
pmr_HL <- (pmr_lisa$quad_r == hl_id_pmr)
par_HL <- (par_lisa$quad_r == hl_id_par)

pmr_HL <- mask(pmr_HL, pop_has_people)
par_HL <- mask(par_HL, pop_has_people)

both_HL <- pmr_HL & par_HL

# ---- 3) 4-class combo: 0=None, 1=PMR-HL only, 2=PAR-HL only, 3=Both HL
combo_HL <- (pmr_HL * 1) + (par_HL * 2)  # logical*int -> 0/1/2/3
combo_HL <- classify(combo_HL, rcl = matrix(c(0,0, 1,1, 2,2, 3,3), ncol = 2, byrow = TRUE))
combo_HL <- as.factor(combo_HL)
levels(combo_HL) <- data.frame(
  ID    = 0:3,
  label = c("None (or POP=0)","PMR HL only","PAR HL only","Both HL")
)

# ---- 4) Plot (x11) — actionable HL hotspots
x11(width = 7, height = 6)
plot(combo_HL,
     main = "Actionable mismatch (LISA: HL = high demand, low supply)\n(POP>0 only)",
     col  = c("grey90", "#FFA866", "#66B2FF", "#FF4D4D"),
     plg  = list(title = "HL classes"))

# ---- 5) Exports — GeoTIFFs (Byte where appropriate)
writeRaster(combo_HL, "CES_LISA_HL_4classes_50m_POPgt0.tif",
            overwrite = TRUE, datatype = "INT1U", gdal = "COMPRESS=LZW")
writeRaster(pmr_HL,   "PMR_LISA_HL_binary_50m_POPgt0.tif",
            overwrite = TRUE, datatype = "INT1U", gdal = "COMPRESS=LZW")
writeRaster(par_HL,   "PAR_LISA_HL_binary_50m_POPgt0.tif",
            overwrite = TRUE, datatype = "INT1U", gdal = "COMPRESS=LZW")
writeRaster(both_HL,  "BOTH_LISA_HL_binary_50m_POPgt0.tif",
            overwrite = TRUE, datatype = "INT1U", gdal = "COMPRESS=LZW")

# (Optional) quick area shares of HL classes (by cell count)
vals <- as.vector(values(combo_HL)); vals <- vals[!is.na(vals)]
tab  <- as.data.frame(table(vals), stringsAsFactors = FALSE)
names(tab) <- c("class_id","count"); tab$class_id <- as.integer(as.character(tab$class_id))
lev  <- levels(combo_HL)[[1]]
share_tbl <- merge(tab, lev, by.x = "class_id", by.y = "ID", all.x = TRUE)
share_tbl$share <- share_tbl$count / sum(share_tbl$count)
cat("\nLISA-HL area shares (POP>0):\n")
print(share_tbl[order(share_tbl$class_id), c("label","count","share")], row.names = FALSE)



# --------------------------------------------
# 4×4 biplot classification on POP>0 cells only
# (no z-scores; quantiles of raw S and D)
# --------------------------------------------
bivar_class_pop <- function(S_raw, D_raw, mask_pop, n = 4) {
  Sx <- terra::mask(S_raw, mask_pop)
  Dx <- terra::mask(D_raw, mask_pop)
  
  # Quantile breaks computed on exposed cells
  qsS <- quantile(terra::values(Sx), probs = seq(0,1,length.out = n+1), na.rm = TRUE)
  qsD <- quantile(terra::values(Dx), probs = seq(0,1,length.out = n+1), na.rm = TRUE)
  
  cutS <- terra::classify(Sx, cbind(qsS[-length(qsS)], qsS[-1], 1:n), include.lowest = TRUE)
  cutD <- terra::classify(Dx, cbind(qsD[-length(qsD)], qsD[-1], 1:n), include.lowest = TRUE)
  
  code <- (cutD - 1L) * n + cutS   # 1..n^2 row-major (D rows × S cols)
  code <- terra::ifel(mask_pop, code, NA)  # NA outside POP>0
  code <- terra::as.factor(code)
  
  # colors (Supply green → x-axis; Demand orange → y-axis)
  low_col <- "#f0f0f0"
  rampS   <- grDevices::colorRampPalette(c(low_col, "#008000"))(n)
  rampD   <- grDevices::colorRampPalette(c(low_col, "#FF5C00"))(n)
  pal_mat <- matrix(NA_character_, n, n)
  for (iy in 1:n) for (ix in 1:n) {
    cx <- grDevices::col2rgb(rampS[ix]); cy <- grDevices::col2rgb(rampD[iy])
    mix <- (cx + cy)/2; pal_mat[iy, ix] <- grDevices::rgb(mix[1]/255, mix[2]/255, mix[3]/255)
  }
  pal_vec <- as.vector(t(pal_mat))  # row-major
  
  # label table
  labs <- as.vector(t(matrix(
    paste0(rep(1:n, each=n), "-", rep(1:n, times=n)),
    nrow = n, byrow = TRUE)))
  terra::levels(code) <- data.frame(ID = 1:(n*n), label = labs)
  
  list(code = code, palette = pal_vec)
}

# Build and plot for PMR~DMR and PAR~DAR (POP>0)
bP <- bivar_class_pop(r_S_PMR, r_D_PMR, pop_has_people, n = 4)
bR <- bivar_class_pop(r_S_PAR, r_D_PAR, pop_has_people, n = 4)

x11(width=7, height=6); plot(bP$code, col=bP$palette, main="Bivariate 4×4 — PMR × DMR (POP>0)", axes=FALSE, plg=FALSE, box=FALSE)
x11(width=7, height=6); plot(bR$code, col=bR$palette, main="Bivariate 4×4 — PAR × DAR (POP>0)", axes=FALSE, plg=FALSE, box=FALSE)

# Optional GeoTIFF export
terra::writeRaster(bP$code, "Bivar4x4_PMRxDMR_POPgt0.tif", overwrite=TRUE, datatype="INT1U", gdal="COMPRESS=LZW")
terra::writeRaster(bR$code, "Bivar4x4_PARxDAR_POPgt0.tif", overwrite=TRUE, datatype="INT1U", gdal="COMPRESS=LZW")

