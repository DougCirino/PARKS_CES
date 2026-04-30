## ============================
## PAR — Recreation model in parks
## ============================

## Packages
library(lme4)
library(lmerTest)
library(ggplot2)
library(dplyr)
library(car)
library(cvTools)

## Working directory and data
setwd('G:/Meu Drive/PESQUISA - Ecologia Urbana e Serviços Ecossistêmicos/Doutorado/Dados/CES_Parques/Statistical Analysis')
data_CES <- read.csv("CES_parks_Final.csv", dec = ',')

## ----------------------------
## 1) Variables
## ----------------------------

# Response
recreation <- as.numeric(data_CES$act_active)

# Landscape / infrastructure (scaled)
lpi1 <- scale(as.numeric(data_CES$X1_lpi_0m))
lpi2 <- scale(as.numeric(data_CES$X2_lpi_0m))
lpi4 <- scale(as.numeric(data_CES$X4_lpi_0m))
lpi5 <- scale(as.numeric(data_CES$X5_lpi_0m))
lpi6 <- scale(as.numeric(data_CES$X6_lpi_0m))

ed1 <- scale(as.numeric(data_CES$X1_ed_0m))
ed2 <- scale(as.numeric(data_CES$X2_ed_0m))
ed4 <- scale(as.numeric(data_CES$X4_ed_0m))
ed5 <- scale(as.numeric(data_CES$X5_ed_0m))
ed6 <- scale(as.numeric(data_CES$X6_ed_0m))

pland1 <- scale(as.numeric(data_CES$X1_pland_0m))
pland2 <- scale(as.numeric(data_CES$X2_pland_0m))
pland4 <- scale(as.numeric(data_CES$X4_pland_0m))
pland5 <- scale(as.numeric(data_CES$X5_pland_0m))
pland6 <- scale(as.numeric(data_CES$X6_pland_0m))
pland_tree <- scale(as.numeric(data_CES$X1_pland_0m)+as.numeric(data_CES$X4_pland_0m))

np1 <- scale(as.numeric(data_CES$X1_np_0m))
np2 <- scale(as.numeric(data_CES$X2_np_0m))
np4 <- scale(as.numeric(data_CES$X4_np_0m))
np5 <- scale(as.numeric(data_CES$X5_np_0m))
np6 <- scale(as.numeric(data_CES$X6_np_0m))
nptree <- scale(as.numeric(data_CES$X1_np_0m)+as.numeric(data_CES$X4_np_0m))

nlsi1 <- scale(as.numeric(data_CES$X1_nlsi_0m))
nlsi2 <- scale(as.numeric(data_CES$X2_nlsi_0m))
nlsi4 <- scale(as.numeric(data_CES$X4_nlsi_0m))
nlsi5 <- scale(as.numeric(data_CES$X5_nlsi_0m))
nlsi6 <- scale(as.numeric(data_CES$X6_nlsi_0m))

contig_mn1 <- scale(as.numeric(data_CES$X1_contig_mn_0m))
contig_mn2 <- scale(as.numeric(data_CES$X2_contig_mn_0m))
contig_mn4 <- scale(as.numeric(data_CES$X4_contig_mn_0m))
contig_mn5 <- scale(as.numeric(data_CES$X5_contig_mn_0m))
contig_mn6 <- scale(as.numeric(data_CES$X6_contig_mn_0m))
contig_mn_tree <- scale(as.numeric(data_CES$contig_mn_tree))

area <- scale(as.numeric(data_CES$pq_area))
area_agua <- scale(as.numeric(data_CES$area_agua)/as.numeric(data_CES$pq_area))
data_CES$presence_agua <- ifelse(data_CES$area_agua == 0, 0, 1)
presence_agua <- data_CES$presence_agua

vol_tree     <- scale(as.numeric(data_CES$vol_asum)/as.numeric(data_CES$pq_area))
vol_tree_tot <- scale(as.numeric(data_CES$vol_asum))
alt_tree     <- scale(as.numeric(data_CES$altmean))

Mobilidade <- scale(as.numeric(data_CES$Mobilidade))
Acessibili <- scale(as.numeric(data_CES$Acessibili))
Equip_cria <- scale(as.numeric(data_CES$Equip_cria))
Equip_espo <- scale(as.numeric(data_CES$Equip_espo))
Sinalizaca <- scale(as.numeric(data_CES$Sinalizaca))
Nota_Infra  <- scale(as.numeric(data_CES$Nota_Infra))
Manejo_Con  <- scale(as.numeric(data_CES$Manejo_Con))
Residuos_s  <- scale(as.numeric(data_CES$Residuos_s))
Ambientes_aqua <- scale(as.numeric(data_CES$Ambientes_))
Nota_manej  <- scale(as.numeric(data_CES$Nota_manej))
Nota_segur  <- scale(as.numeric(data_CES$Nota_segur))
Nota_gesta  <- scale(as.numeric(data_CES$Nota_gesta))
Nota_geral  <- scale(Nota_Infra + Nota_gesta + Nota_manej + Nota_segur)

shannon_land <- scale(as.numeric(data_CES$shannon_index))
total_edge   <- scale(as.numeric(data_CES$total_edge))
contig       <- scale(as.numeric(data_CES$contig))
lsi          <- scale(as.numeric(data_CES$lsi))
circle       <- scale(as.numeric(data_CES$circle))

# Controls / people
park   <- as.factor(data_CES$pq_nome)
GeoSES <- scale(as.numeric(data_CES$Media_GeoSES))
income <- scale(as.numeric(data_CES$Media_renda))
gender <- as.factor(data_CES$female)

data_CES$income_p <- ifelse(data_CES$income == "10 to 20 MW", 15,
                            ifelse(data_CES$income == "2 to 4 MW", 3,
                                   ifelse(data_CES$income == "> 20 MW", 20,
                                          ifelse(data_CES$income == "< 2 MW", 2,
                                                 ifelse(data_CES$income == "4 to 10 MW", 7, NA)))))
income_p <- scale(as.numeric(data_CES$income_p))

age <- scale(as.numeric(data_CES$age))

data_CES$frequency_p <- ifelse(data_CES$frequency == "1 or 2 times a week", 6,
                               ifelse(data_CES$frequency == "3 to 6 times a week", 18,
                                      ifelse(data_CES$frequency == "Every day", 30,
                                             ifelse(data_CES$frequency == "Once a month", 1,
                                                    ifelse(data_CES$frequency == "1 or 2 times a year", 0.125,
                                                           ifelse(data_CES$frequency == "Never/First time", 0.01,
                                                                  ifelse(data_CES$frequency == "Every 2 or 3 months", 0.4,
                                                                         ifelse(data_CES$frequency == "2 or 3 times a month", 2.5, NA))))))))
data_CES$duration_p <- ifelse(data_CES$duration == "1 to 2h", 120,
                              ifelse(data_CES$duration == "30 min to 1h", 60,
                                     ifelse(data_CES$duration == ">3h", 180,
                                            ifelse(data_CES$duration == "2 to 3h", 150,
                                                   ifelse(data_CES$duration == "15 to 30 min", 30,
                                                          ifelse(data_CES$duration == "<15 min", 15, NA))))))
nature_dose <- scale(as.numeric(data_CES$frequency_p) * as.numeric(data_CES$duration_p))

## Modeling frame (row = person)
predictor_df <- data.frame(
  lpi1,lpi4,lpi6,
  ed1,ed4,ed6,
  pland1,pland4,pland6, pland_tree,
  np1,np2,np4,np6,nptree,
  nlsi1,nlsi4,nlsi6,
  contig_mn1,contig_mn4,contig_mn6,
  area,area_agua,presence_agua, vol_tree,vol_tree_tot,alt_tree,
  Mobilidade,Acessibili,Equip_cria,Equip_espo,Sinalizaca,
  Nota_Infra,Manejo_Con,Residuos_s,Ambientes_aqua,
  Nota_manej,Nota_segur,Nota_gesta,Nota_geral,
  shannon_land,total_edge,contig,lsi,
  circle,contig_mn_tree,nptree,
  income,GeoSES,gender,income_p,age,nature_dose
)
model_data <- data.frame(predictor_df, recreation, park)
model_data <- model_data[complete.cases(model_data$recreation), ]

## ----------------------------
## 2) Mixed model and park signal
## ----------------------------

model_base_PAR <- lm(recreation ~ park, data = model_data)

m1_PAR <- glmer(recreation ~ age + income_p + nature_dose + (1|park),
                data = model_data, family = binomial(link="logit"))
summary(m1_PAR)

ri_PAR <- ranef(m1_PAR)$park
ri_PAR_df <- data.frame(park = rownames(ri_PAR),
                        random_intercept = as.numeric(ri_PAR[,1]))

## ----------------------------
## 3) Collapse to one row per park and fit park-level LM
## ----------------------------

selected_columns <- c("park","pland_tree","circle","Nota_manej","shannon_land","contig_mn_tree",
                      "Residuos_s","pland6","lpi1","lpi4","lpi6","ed1","ed4","ed6","pland1",
                      "pland4","np1","np2","np4","np6","nptree","nlsi1","nlsi4","nlsi6",
                      "contig_mn1","contig_mn4","contig_mn6","area","area_agua","presence_agua",
                      "vol_tree","vol_tree_tot","alt_tree","Mobilidade","Acessibili","Equip_cria",
                      "Equip_espo","Sinalizaca","Nota_Infra","Manejo_Con","Ambientes_aqua",
                      "Nota_segur","Nota_gesta","Nota_geral","total_edge","contig","lsi","circle",
                      "contig_mn_tree","nptree","income","GeoSES")

model_data_selected_PAR <- model_data %>%
  dplyr::select(dplyr::all_of(selected_columns)) %>%
  dplyr::distinct(park, .keep_all = TRUE)

merged_data_PAR <- model_data_selected_PAR %>%
  dplyr::left_join(ri_PAR_df, by = "park")

model_adjusted_PAR <- lm(random_intercept ~ circle + Equip_espo + np6 + Ambientes_aqua,
                         data = merged_data_PAR, na.action = na.exclude)
summary(model_adjusted_PAR); vif(model_adjusted_PAR)

## ----------------------------
## 4) Effect-size plot (from the hard-coded table)
## ----------------------------
model_data_es <- data.frame(
  Variable  = c("(Intercept)","Patch elongation (Circle)","Sports infrastructure","NP herbaceous-shrub","Quality of aquatic environments"),
  Estimate  = c(-0.03276, -0.16936, 0.13821, -0.11900, 0.09224),
  Std_Error = c(0.04146, 0.05144, 0.05010, 0.04688, 0.04261),
  t_value   = c(-0.790, -3.292, 2.759, -2.539, 2.165),
  p_value   = c(0.44267, 0.00534, 0.01538, 0.02364, 0.04819)
) %>% mutate(Significance = case_when(
  p_value < 0.001 ~ "***",
  p_value < 0.01  ~ "**",
  p_value < 0.05  ~ "*",
  p_value < 0.1   ~ ".",
  TRUE ~ ""
))
x11()
ggplot(model_data_es, aes(x = Estimate, y = reorder(Variable, Estimate))) +
  geom_point() +
  geom_errorbarh(aes(xmin = Estimate - Std_Error, xmax = Estimate + Std_Error), height = 0.2) +
  geom_text(aes(label = Significance), vjust = -0.5, hjust = -0.2) +
  labs(x = "Effect Size (Estimate)", y = "Variable") +
  theme_minimal()

## ----------------------------
## 5) Effect-curve plots
## ----------------------------
p1 <- ggplot(merged_data_PAR, aes(x = circle, y = random_intercept)) +
  geom_point() + geom_smooth(method="lm", se=FALSE, color="blue") +
  labs(title="Effect of Circle on Random Intercept (PAR)", x="Circle", y="Random Intercept")
p2 <- ggplot(merged_data_PAR, aes(x = Equip_espo, y = random_intercept)) +
  geom_point() + geom_smooth(method="lm", se=FALSE, color="blue") +
  labs(title="Effect of Sports infrastructure on Random Intercept (PAR)", x="Equip_espo", y="Random Intercept")
p3 <- ggplot(merged_data_PAR, aes(x = np6, y = random_intercept)) +
  geom_point() + geom_smooth(method="lm", se=FALSE, color="blue") +
  labs(title="Effect of NP Herbaceous–Shrub on Random Intercept (PAR)", x="NP Herbaceous–Shrub", y="Random Intercept")
p4 <- ggplot(merged_data_PAR, aes(x = Ambientes_aqua, y = random_intercept)) +
  geom_point() + geom_smooth(method="lm", se=FALSE, color="blue") +
  labs(title="Effect of Aquatic environments on Random Intercept (PAR)", x="Aquatic management score", y="Random Intercept")
x11(); gridExtra::grid.arrange(p1,p2,p3,p4,ncol=2)

## ----------------------------
## 6) Extrapolation to all parks (+ CI bars)
## ----------------------------

# Scaling params from the original sample (for 4 predictors)
equip_espo_scaled   <- scale(as.numeric(data_CES$Equip_espo))
circle_scaled       <- scale(as.numeric(data_CES$circle))
np6_scaled          <- scale(as.numeric(data_CES$X6_np_0m))
amb_aqua_scaled     <- scale(as.numeric(data_CES$Ambientes_))
scaling_params_PAR <- list(
  equip_espo_mean = attr(equip_espo_scaled, "scaled:center"),
  equip_espo_sd   = attr(equip_espo_scaled, "scaled:scale"),
  circle_mean     = attr(circle_scaled, "scaled:center"),
  circle_sd       = attr(circle_scaled, "scaled:scale"),
  np6_mean        = attr(np6_scaled, "scaled:center"),
  np6_sd          = attr(np6_scaled, "scaled:scale"),
  amb_mean        = attr(amb_aqua_scaled, "scaled:center"),
  amb_sd          = attr(amb_aqua_scaled, "scaled:scale")
)

all_parks <- read.csv("all_parks_PRS_error.csv", dec = ',')
all_parks$Equip_espo     <- scale(as.numeric(all_parks$Equip_espo), center=scaling_params_PAR$equip_espo_mean, scale=scaling_params_PAR$equip_espo_sd)
all_parks$circle         <- scale(as.numeric(all_parks$circle),     center=scaling_params_PAR$circle_mean,     scale=scaling_params_PAR$circle_sd)
all_parks$np6            <- scale(as.numeric(all_parks$X6_np_0m),   center=scaling_params_PAR$np6_mean,        scale=scaling_params_PAR$np6_sd)
all_parks$Ambientes_aqua <- scale(as.numeric(all_parks$Ambientes_), center=scaling_params_PAR$amb_mean,        scale=scaling_params_PAR$amb_sd)

# Prediction and simple SE for CI bars
fe_PAR <- coef(model_adjusted_PAR)
X <- model.matrix(~ Equip_espo + circle + np6 + Ambientes_aqua, data = all_parks)
vcov_PAR <- vcov(model_adjusted_PAR)
pred_se  <- sqrt(rowSums((X %*% vcov_PAR) * X))
all_parks$pred_se_rec <- pred_se

all_parks$predicted_recreation <- fe_PAR["(Intercept)"] +
  fe_PAR["Equip_espo"]     * all_parks$Equip_espo +
  fe_PAR["circle"]         * all_parks$circle +
  fe_PAR["np6"]            * all_parks$np6 +
  fe_PAR["Ambientes_aqua"] * all_parks$Ambientes_aqua

# Add intercept from m1_PAR (the original step) and clamp at 0
intercept_m1_PAR <- fixef(m1_PAR)["(Intercept)"]
all_parks$predicted_recreation <- pmax(all_parks$predicted_recreation + intercept_m1_PAR, 0)

# 95% CI
z_value <- qnorm(0.975)
all_parks$conf_lower_rec <- all_parks$predicted_recreation - z_value * pred_se
all_parks$conf_upper_rec <- all_parks$predicted_recreation + z_value * pred_se

# Rank plot (ordered by predicted)
all_parks <- all_parks %>% mutate(park_name = factor(park_name, levels = park_name[order(predicted_recreation)]))
x11()
ggplot(all_parks, aes(x = park_name, y = predicted_recreation)) +
  geom_point() +
  geom_errorbar(aes(ymin = conf_lower_rec, ymax = conf_upper_rec), width = 0) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Park", y = "Predicted recreation")

# Ranking table for export/merge later
predicted_recreation_table_PAR <- data.frame(
  park_name = as.character(all_parks$park_name),
  recreation = all_parks$predicted_recreation
)

write.csv(all_parks, file = "all_parks_PRS_RECREATION_error.csv", row.names = FALSE)

## -----------------------------------------------
## STEP 6 — PAR: park-signal only (NO correction)
## Normalize predictions & CIs to [0,1]
## -----------------------------------------------

# ---- Safety: point to the right LM object
if (exists("model_adjusted_PAR")) {
  lm_PAR_final <- model_adjusted_PAR
} else if (exists("model_adjusted")) {
  lm_PAR_final <- model_adjusted
} else {
  stop("Could not find 'model_adjusted_PAR' or 'model_adjusted'.")
}

# ---- Ensure park_name exists
if (!"park_name" %in% names(all_parks)) {
  if ("park" %in% names(all_parks)) {
    all_parks$park_name <- as.character(all_parks$park)
  } else {
    stop("all_parks must have 'park_name' or 'park'.")
  }
}

# ---- Linear predictions (log-odds scale, park signal only)
X_par <- model.matrix(~ circle + Equip_espo + np6 + Ambientes_aqua, data = all_parks)
b_par <- coef(lm_PAR_final)
V_par <- vcov(lm_PAR_final)

par_mean_signal <- as.numeric(X_par %*% b_par)
par_se_signal   <- sqrt(rowSums((X_par %*% V_par) * X_par))
zv <- qnorm(0.975)

par_lo_signal <- par_mean_signal - zv * par_se_signal
par_hi_signal <- par_mean_signal + zv * par_se_signal

# ---- Min–max normalization to [0,1] using CI envelope
min_all <- min(par_lo_signal, na.rm = TRUE)
max_all <- max(par_hi_signal, na.rm = TRUE)
rng     <- max_all - min_all
if (rng <= 0) stop("Non-positive range during normalization.")

norm01 <- function(x) (x - min_all) / rng

PAR_norm_mean <- norm01(par_mean_signal)
PAR_norm_lo   <- pmax(0, pmin(1, norm01(par_lo_signal)))  # clamp to [0,1]
PAR_norm_hi   <- pmax(0, pmin(1, norm01(par_hi_signal)))  # clamp to [0,1]

par_pred_df_signal <- data.frame(
  park_name     = all_parks$park_name,
  PAR_norm_mean = PAR_norm_mean,
  PAR_norm_lo   = PAR_norm_lo,
  PAR_norm_hi   = PAR_norm_hi
)

# ---- Reorder by normalized mean
par_pred_df_signal$park_name <- factor(
  par_pred_df_signal$park_name,
  levels = par_pred_df_signal$park_name[order(par_pred_df_signal$PAR_norm_mean)]
)


  # ---- Plot (same style as PMR; y in 0–1) ----
  x11()
  ggplot(par_pred_df_signal, aes(x = park_name, y = PAR_norm_mean)) +
    geom_point(size = 1.3) +
    geom_errorbar(aes(ymin = PAR_norm_lo, ymax = PAR_norm_hi), width = 0.2) +
    labs(x = "Park",
         y = "Potential Active Recreation (normalized 0–1)") +
    coord_cartesian(ylim = c(0, 1)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank())
  
  ## ----------------------------
  ## 7) Simple validation — Observed vs Predicted (19 parks only)
  ##    - Uses the known per-park random-intercept coefficients (recreation_coef_values)
  ##    - Fits the park-level LM on those 19 parks and compares predictions
  ## ----------------------------
  
  # 19-park observed 'recreation_coef' (the provided values)
  recreation_coef_values  <- c(
    `ALFREDO VOLPI` = -0.31471829,
    `BARRAGEM DE GUARAPIRANGA` = -0.12061018,
    `BENEMERITO JOSE BRAS` = -0.20233742,
    `CIDADE DE TORONTO` =  0.02943055,
    `ECOLOGICO PROFESSORA LYDIA NATALIZIO DIOGO` =  0.09109760,
    `GUARAPIRANGA` =  0.44621809,
    `IBIRAPUERA` = -0.28503448,
    `JARDIM FELICIDADE` =  0.16352497,
    `LAJEADO - IZAURA PEREIRA DE SOUZA FRANZOLIN` = -0.10498509,
    `LINEAR GUARATIBA` =  0.36200885,
    `LINEAR RIO VERDE` =  0.39542680,
    `LUZ` =  0.09893753,
    `PARQUE DO POVO` = -0.10585707,
    `RAPOSO TAVARES` = -0.21140285,
    `SANTA AMELIA` = -0.09479459,
    `SANTO DIAS` =  0.34817191,
    `SAPOPEMBA (ATERRO)` =  0.11250226,
    `SEVERO GOMES` =  0.03133373,
    `TENENTE SIQUEIRA CAMPOS (TRIANON)` = -0.88559829
  )
  
  # Attach observed coefficients to all_parks
  if (!"park_name" %in% names(all_parks)) {
    if ("park" %in% names(all_parks)) {
      all_parks$park_name <- as.character(all_parks$park)
    } else stop("STEP 7: 'all_parks' must contain 'park_name' or 'park'.")
  }
  
  all_parks$recreation_coef <- NA_real_
  all_parks$recreation_coef[match(names(recreation_coef_values), all_parks$park_name)] <- recreation_coef_values
  
  # Keep only the 19 parks with ground-truth recreation_coef
  val_df <- all_parks %>% 
    dplyr::filter(!is.na(recreation_coef)) %>%
    dplyr::select(park_name, circle, Equip_espo, np6, Ambientes_aqua, recreation_coef)
  
  # Fit the *park-level* LM on those 19 (same formula used for extrapolation)
  lm_val_PAR <- lm(recreation_coef ~ circle + Equip_espo + np6 + Ambientes_aqua,
                   data = val_df, na.action = na.exclude)
  
  # Predict and compute R² (in-sample for these 19)
  val_df$predicted_recreation_coef <- predict(lm_val_PAR, newdata = val_df)
  rsq_val <- summary(lm_val_PAR)$r.squared
  
  cat("\n--- STEP 7 (PAR) — Observed vs Predicted on 19 parks ---\n")
  cat("R² (training, 19 parks) =", round(rsq_val, 3), "\n")
  print(summary(lm_val_PAR)$coefficients)
  
  # Plot Observed vs Predicted (19 parks)
  x11()
  ggplot(val_df, aes(x = predicted_recreation_coef, y = recreation_coef, label = park_name)) +
    geom_point(color = "blue", size = 2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    labs(x = "Predicted PAR",
         y = "Observed PAR",
         title = "Observed vs Predicted (PAR)") +
    annotate("text",
             x = min(val_df$predicted_recreation_coef, na.rm = TRUE),
             y = max(val_df$recreation_coef, na.rm = TRUE),
             label = paste0("R² = ", round(rsq_val, 3)),
             hjust = 0, vjust = 1, color = "red", size = 4) +
    theme_minimal()
  
## =========================================================
## STEP 8 — PAR cross-validation of park-level model (as-is)
##  - Uses the cv_function logic:
##    * MSE on TEST folds (lower is better)
##    * R² from TRAINING fit (not CV R²; higher is better in-sample)
##    * TRAINING coefficient CIs and SEs per fold
## =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(cvTools)
})

# ---- 1) Guard: ensure we have recreation_coef and a model formula
if (!"recreation_coef" %in% names(all_parks) || all(all_parks$recreation_coef %in% c(NA, NaN))) {
  # Fill from the provided lookup (kept here for robustness)
  recreation_coef_values  <- c(
    `ALFREDO VOLPI` = -0.31471829,
    `BARRAGEM DE GUARAPIRANGA` = -0.12061018,
    `BENEMERITO JOSE BRAS` = -0.20233742,
    `CIDADE DE TORONTO` =  0.02943055,
    `ECOLOGICO PROFESSORA LYDIA NATALIZIO DIOGO` =  0.09109760,
    `GUARAPIRANGA` =  0.44621809,
    `IBIRAPUERA` = -0.28503448,
    `JARDIM FELICIDADE` =  0.16352497,
    `LAJEADO - IZAURA PEREIRA DE SOUZA FRANZOLIN` = -0.10498509,
    `LINEAR GUARATIBA` =  0.36200885,
    `LINEAR RIO VERDE` =  0.39542680,
    `LUZ` =  0.09893753,
    `PARQUE DO POVO` = -0.10585707,
    `RAPOSO TAVARES` = -0.21140285,
    `SANTA AMELIA` = -0.09479459,
    `SANTO DIAS` =  0.34817191,
    `SAPOPEMBA (ATERRO)` =  0.11250226,
    `SEVERO GOMES` =  0.03133373,
    `TENENTE SIQUEIRA CAMPOS (TRIANON)` = -0.88559829
  )
  if (!"park_name" %in% names(all_parks)) {
    stop("all_parks must contain park_name to map recreation_coef.")
  }
  all_parks$recreation_coef <- NA_real_
  all_parks$recreation_coef[match(names(recreation_coef_values), all_parks$park_name)] <- recreation_coef_values
}

# This is the park-level PAR LM you’ve been using:
#   recreation_coef ~ circle + Equip_espo + np6 + Ambientes_aqua
model_formula <- recreation_coef ~ circle + Equip_espo + np6 + Ambientes_aqua

# ---- 2) CV function (kept consistent with the PMR version)
cv_function <- function(train_indices, test_indices, data, formula) {
  train_data <- data[train_indices, ]
  test_data  <- data[test_indices, ]
  fit <- lm(formula, data = train_data, na.action = na.exclude)
  
  test_data$predicted_recreation <- predict(fit, newdata = test_data)
  mse <- mean((test_data$recreation_coef - test_data$predicted_recreation)^2, na.rm = TRUE)
  
  # Training R² (note: not CV R²)
  rsq_in <- summary(fit)$r.squared
  
  list(
    MSE = mse,
    R_squared = rsq_in,
    Confidence_Intervals = confint(fit),
    SE = summary(fit)$coefficients[, "Std. Error"]
  )
}

# ---- 3) Run K-fold CV
set.seed(123)
n_folds <- 5L
folds <- cvTools::cvFolds(nrow(all_parks), K = n_folds)

cv_results <- lapply(seq_len(n_folds), function(i) {
  train_idx <- folds$subsets[folds$which != i]
  test_idx  <- folds$subsets[folds$which == i]
  cv_function(train_idx, test_idx, all_parks, model_formula)
})

# ---- 4) Summaries: per-fold table + overall means/sd
fold_metrics <- dplyr::tibble(
  Fold      = seq_len(n_folds),
  MSE       = vapply(cv_results, function(x) x$MSE, numeric(1)),
  R2_train  = vapply(cv_results, function(x) x$R_squared, numeric(1))
)

overall_summary <- dplyr::tibble(
  Metric = c("MSE (test)", "R² (training)"),
  Mean   = c(mean(fold_metrics$MSE, na.rm = TRUE),
             mean(fold_metrics$R2_train, na.rm = TRUE)),
  SD     = c(sd(fold_metrics$MSE, na.rm = TRUE),
             sd(fold_metrics$R2_train, na.rm = TRUE))
)

cat("\n--- PAR CV (park-level LM) — per-fold metrics ---\n")
print(fold_metrics, n = Inf)

cat("\n--- PAR CV (park-level LM) — overall summary ---\n")
print(overall_summary)

## =========================================================
## STEP 9 — PAR quick diagnostics
##  - Inputs: all_parks with recreation_coef (19 parks filled) and
##            PAR park-level formula: recreation_coef ~ circle + Equip_espo + np6 + Ambientes_aqua
##  - Outputs:
##    (a) Observed vs Predicted (with 1:1 and calibration fit)
##    (b) Rank agreement (Spearman rho)
##    (c) Residuals vs Fitted
## =========================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
})

# 1) Fit on the 19 measured parks only
par_formula <- recreation_coef ~ circle + Equip_espo + np6 + Ambientes_aqua
par19 <- all_parks %>% filter(!is.na(recreation_coef))
par_fit <- lm(par_formula, data = par19, na.action = na.exclude)

# Predictions & residuals
par19$pred <- predict(par_fit, newdata = par19)
par19$res  <- par19$recreation_coef - par19$pred

# 2) (a) Observed vs Predicted
x11()
ggplot(par19, aes(x = pred, y = recreation_coef, label = park_name)) +
  geom_point(size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = 50) +
  labs(x = "Predicted PAR (park signal)", y = "Observed PAR (recreation_coef)",
       title = "PAR — Observed vs Predicted (19 parks)") +
  theme_bw()

# 2) (b) Rank agreement
par19 <- par19 %>%
  mutate(rank_obs = rank(recreation_coef, ties.method = "average"),
         rank_pred = rank(pred,            ties.method = "average"))

rho_par <- suppressWarnings(cor(par19$rank_obs, par19$rank_pred, method = "spearman", use = "complete.obs"))

x11()
ggplot(par19, aes(x = rank_pred, y = rank_obs, label = park_name)) +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = 50) +
  labs(x = "Predicted rank", y = "Observed rank",
       title = paste0("PAR — Rank agreement (Spearman ρ = ", sprintf("%.2f", rho_par), ")")) +
  theme_bw()

# 2) (c) Residuals vs Fitted
x11()
ggplot(par19, aes(x = pred, y = res, label = park_name)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = 50) +
  labs(x = "Fitted (predicted)", y = "Residual (obs - pred)",
       title = "PAR — Residuals vs Fitted") +
  theme_bw()


# ---- Spearman rank agreement for PAR (observed vs predicted) ----

sp_par <- suppressWarnings(cor.test(par19$recreation_coef, par19$pred,
                                    method = "spearman", exact = FALSE))

rho_par <- unname(sp_par$estimate)
p_par   <- sp_par$p.value

cat(sprintf("\nPAR — Spearman rank correlation: ρ = %.3f (p = %.4f)\n",
            rho_par, p_par))

# Optional: compact row for the validation summary table
par_spearman_row <- data.frame(
  CES = "PAR",
  Spearman_rho = rho_par,
  p_value = p_par,
  stringsAsFactors = FALSE
)
print(par_spearman_row)


