## ============================
## PMR — PRS model in parks
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
data_CES <- read.csv("CES_parks_final.csv", dec = ',')

## ----------------------------
## 1) Variables
## ----------------------------

# Response
prs <- as.numeric(data_CES$prs_score)

# Landscape metrics (scaled)
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
  area,area_agua,presence_agua,
  vol_tree,vol_tree_tot,alt_tree,
  Mobilidade,Acessibili,Equip_cria,Equip_espo,Sinalizaca,
  Nota_Infra,Manejo_Con,Residuos_s,Ambientes_aqua,
  Nota_manej,Nota_segur,Nota_gesta,Nota_geral,
  shannon_land,total_edge,contig,lsi,
  circle,contig_mn_tree,nptree,
  income,GeoSES,gender,income_p,age,nature_dose
)
model_data <- data.frame(predictor_df, prs, park)
model_data <- model_data[complete.cases(model_data$prs), ]

## ----------------------------
## 2) Mixed model and park signal
## ----------------------------

# Baseline check (fixed park as factor)
model_base_PMR <- lm(prs ~ park, data = model_data)

# Mixed model with park as random intercept (people-level controls)
m1_PMR <- lmer(prs ~ age + income_p + nature_dose + (1|park), data = model_data)
summary(m1_PMR)

# Random intercepts by park
ri_PMR <- ranef(m1_PMR)$park
ri_PMR_df <- data.frame(park = rownames(ri_PMR),
                        random_intercept = as.numeric(ri_PMR[,1]))

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

model_data_selected_PMR <- model_data %>%
  dplyr::select(dplyr::all_of(selected_columns)) %>%
  dplyr::distinct(park, .keep_all = TRUE)

merged_data_PMR <- model_data_selected_PMR %>%
  dplyr::left_join(ri_PMR_df, by = "park")

# Park-level LM explaining park signal
model_adjusted_PMR <- lm(random_intercept ~ circle + Nota_manej + shannon_land + contig_mn_tree,
                         data = merged_data_PMR, na.action = na.exclude)
summary(model_adjusted_PMR); vif(model_adjusted_PMR)

## ----------------------------
## 4) Effect-size plot (from the hard-coded table)
## ----------------------------
model_data_es <- data.frame(
  Variable  = c("(Intercept)","Circle","Management score","Shannon landscape","Contiguity of tree coverage"),
  Estimate  = c(0.01420, 0.16592, 0.22872, 0.22989, 0.15324),
  Std_Error = c(0.06588, 0.07785, 0.07316, 0.06961, 0.07081),
  t_value   = c(0.215, 2.131, 3.126, 3.303, 2.164),
  p_value   = c(0.83249, 0.05127, 0.00744, 0.00524, 0.04823)
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
## 5) Effect-curve plots (no grid removal requested)
## ----------------------------
p1 <- ggplot(merged_data_PMR, aes(x = circle, y = random_intercept)) +
  geom_point() + geom_smooth(method="lm", se=FALSE, color="blue") +
  labs(title="Effect of Circle on Random Intercept (PMR)", x="Circle", y="Random Intercept")
p2 <- ggplot(merged_data_PMR, aes(x = Nota_manej, y = random_intercept)) +
  geom_point() + geom_smooth(method="lm", se=FALSE, color="blue") +
  labs(title="Effect of Management score on Random Intercept (PMR)", x="Management score", y="Random Intercept")
p3 <- ggplot(merged_data_PMR, aes(x = shannon_land, y = random_intercept)) +
  geom_point() + geom_smooth(method="lm", se=FALSE, color="blue") +
  labs(title="Effect of Shannon landscape on Random Intercept (PMR)", x="Shannon_land", y="Random Intercept")
p4 <- ggplot(merged_data_PMR, aes(x = contig_mn_tree, y = random_intercept)) +
  geom_point() + geom_smooth(method="lm", se=FALSE, color="blue") +
  labs(title="Effect of Tree contiguity on Random Intercept (PMR)", x="Contig_mn_tree", y="Random Intercept")
x11(); gridExtra::grid.arrange(p1,p2,p3,p4,ncol=2)


## ----------------------------
## 6) Extrapolation to all parks
## ----------------------------

# Store scaling parameters from the original sample (for the 4 predictors in the park LM)
nota_manej_scaled   <- scale(as.numeric(data_CES$Nota_manej))
circle_scaled       <- scale(as.numeric(data_CES$circle))
shannon_land_scaled <- scale(as.numeric(data_CES$shannon_index))
contig_tree_scaled  <- scale(as.numeric(data_CES$contig_mn_tree))
scaling_params_PMR <- list(
  nota_manej_mean = attr(nota_manej_scaled, "scaled:center"),
  nota_manej_sd   = attr(nota_manej_scaled, "scaled:scale"),
  circle_mean     = attr(circle_scaled, "scaled:center"),
  circle_sd       = attr(circle_scaled, "scaled:scale"),
  shannon_mean    = attr(shannon_land_scaled, "scaled:center"),
  shannon_sd      = attr(shannon_land_scaled, "scaled:scale"),
  contig_mean     = attr(contig_tree_scaled, "scaled:center"),
  contig_sd       = attr(contig_tree_scaled, "scaled:scale")
)

# Load all parks and scale them with the same parameters
all_parks <- read.csv("all_parks.csv", dec = ',')
all_parks$Nota_manej     <- scale(as.numeric(all_parks$Nota_manej),     center=scaling_params_PMR$nota_manej_mean, scale=scaling_params_PMR$nota_manej_sd)
all_parks$circle         <- scale(as.numeric(all_parks$circle),         center=scaling_params_PMR$circle_mean,     scale=scaling_params_PMR$circle_sd)
all_parks$shannon_land   <- scale(as.numeric(all_parks$shannon_index),  center=scaling_params_PMR$shannon_mean,    scale=scaling_params_PMR$shannon_sd)
all_parks$contig_mn_tree <- scale(as.numeric(all_parks$contig_mn_tree), center=scaling_params_PMR$contig_mean,     scale=scaling_params_PMR$contig_sd)

# Predict park signal with the park-level LM
fe_PMR <- coef(model_adjusted_PMR)
all_parks$predicted_prs_core <- fe_PMR["(Intercept)"] +
  fe_PMR["Nota_manej"]     * all_parks$Nota_manej +
  fe_PMR["circle"]         * all_parks$circle +
  fe_PMR["shannon_land"]   * all_parks$shannon_land +
  fe_PMR["contig_mn_tree"] * all_parks$contig_mn_tree

# Add the mixed-model intercept from m1_PMR to bring predictions to the PRS scale (the original step)
intercept_value_PMR <- fixef(m1_PMR)["(Intercept)"]
all_parks$final_PRS <- all_parks$predicted_prs_core + intercept_value_PMR

# Ranking table
predicted_prs_table_PMR <- data.frame(
  park_name = all_parks$park_name,
  PRS = all_parks$final_PRS
) %>% arrange(desc(PRS))
print(predicted_prs_table_PMR)


## =========================================================
## STEP 6 — PMR: corrected predictions per park + 95% CIs
## Uses: model_adjusted_PMR, m1_PMR, all_parks
## Predictors in all_parks must be the scaled versions:
##   Nota_manej, circle, shannon_land, contig_mn_tree
## =========================================================

# 1) Build the model matrix that matches the park-level LM
X_PMR <- model.matrix(~ Nota_manej + circle + shannon_land + contig_mn_tree,
                      data = all_parks)

# 2) Get coefficients and covariance from the park-level LM
beta_PMR <- coef(model_adjusted_PMR)
V_PMR    <- vcov(model_adjusted_PMR)

# 3) Linear prediction (without mixed-model intercept) and its SE
pred_core_PMR <- as.numeric(X_PMR %*% beta_PMR)
pred_se_PMR   <- sqrt(rowSums((X_PMR %*% V_PMR) * X_PMR))  # == sqrt(diag(X V X^T))
z97 <- qnorm(0.975)

# 4) Add the mixed-model intercept from m1_PMR to mean and CI bounds
mm_intercept <- fixef(m1_PMR)["(Intercept)"]

all_parks$PMR_mean <- pred_core_PMR + mm_intercept
all_parks$PMR_lo   <- pred_core_PMR - z97 * pred_se_PMR + mm_intercept
all_parks$PMR_hi   <- pred_core_PMR + z97 * pred_se_PMR + mm_intercept

# 5) Produce the table (sorted, descending) with the corrected PRS
predicted_prs_table <- all_parks |>
  dplyr::transmute(
    park_name = park_name,
    PRS       = PMR_mean,          # corrected mean (includes mixed-model intercept)
    PRS_lo    = PMR_lo,
    PRS_hi    = PMR_hi
  ) |>
  dplyr::arrange(dplyr::desc(PRS))

print(predicted_prs_table)

# 6) Plot: Predicted PRS by park with 95% CI, ordered by corrected mean
plot_df <- all_parks |>
  dplyr::mutate(park_name = factor(park_name, levels = park_name[order(PMR_mean)]))

x11()
ggplot(plot_df, aes(x = park_name, y = PMR_mean)) +
  geom_point(size = 1.3) +
  geom_errorbar(aes(ymin = PMR_lo, ymax = PMR_hi), width = 0.2) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),                 # no grid lines
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(x = "Park", y = "Predicted PRS")


## ----------------------------
## 7) Simple validation — Observed vs Predicted (19 parks only)
## ----------------------------

# Observed signal for a subset of parks (the provided coefficients)
prs_coef_values <- c(
  `ALFREDO VOLPI` = 0.84895609, `BARRAGEM DE GUARAPIRANGA` = 0.26749544,
  `BENEMERITO JOSE BRAS` = -0.51621363, `CIDADE DE TORONTO` = 0.14638793,
  `ECOLOGICO PROFESSORA LYDIA NATALIZIO DIOGO` = 0.08426581, `GUARAPIRANGA` = 0.49146719,
  `IBIRAPUERA` = 0.76496225, `JARDIM FELICIDADE` = -0.22130796,
  `LAJEADO - IZAURA PEREIRA DE SOUZA FRANZOLIN` = 0.11296983, `LINEAR GUARATIBA` = -0.67917508,
  `LINEAR RIO VERDE` = -1.01724127, `LUZ` = -0.06621319, `PARQUE DO POVO` = 0.30017353,
  `RAPOSO TAVARES` = -0.35517816, `SANTA AMELIA` = -0.39062124, `SANTO DIAS` = -0.04224277,
  `SAPOPEMBA (ATERRO)` = -0.48456852, `SEVERO GOMES` = 0.50726665,
  `TENENTE SIQUEIRA CAMPOS (TRIANON)` = 0.24881712
)

# Attach observed to all_parks (NA for parks not in the list)
all_parks$prs_coef <- NA_real_
all_parks$prs_coef[match(names(prs_coef_values), all_parks$park_name)] <- prs_coef_values

# Fit simple LM on those parks only and plot Observed vs Predicted
obs_pred_df <- all_parks %>%
  filter(!is.na(prs_coef)) %>%
  dplyr::select(prs_coef, circle, Nota_manej, shannon_land, contig_mn_tree)

lm_val_PMR <- lm(prs_coef ~ circle + Nota_manej + shannon_land + contig_mn_tree,
                 data = obs_pred_df, na.action = na.exclude)
obs_pred_df$predicted <- predict(lm_val_PMR, newdata = obs_pred_df)
R2_PMR <- summary(lm_val_PMR)$r.squared

x11()
ggplot(obs_pred_df, aes(x = predicted, y = prs_coef)) +
  geom_point(color = "blue", size = 2) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(x = "Predicted PRS", y = "Observed PRS", title = "Observed vs Predicted (PMR)") +
  annotate("text", x = max(obs_pred_df$predicted), y = min(obs_pred_df$prs_coef),
           label = paste("R² =", round(R2_PMR, 3)), hjust = 1, vjust = 0, color = "red", size = 5)+
  theme_minimal()

## =========================================================
## STEP 8 — PMR cross-validation of park-level model (as-is)
##  - Uses the cv_function exactly as written:
##    * MSE on TEST folds (lower is better)
##    * R² from TRAINING fit (not CV R²; higher is better in-sample)
##    * TRAINING coefficient CIs and SEs per fold
##  - Adds an automatic interpretation block:
##    * Compares mean MSE to baseline variance of prs_coef
##    * Flags stability via SD across folds
## =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(cvTools)
})

# ---- 0) Guard: ensure required columns exist (prs_coef + predictors) ----
if (!all(c("prs_coef","circle","Nota_manej","shannon_land","contig_mn_tree") %in% names(all_parks))) {
  stop("all_parks must have: prs_coef, circle, Nota_manej, shannon_land, contig_mn_tree.")
}
if (!"park_name" %in% names(all_parks)) {
  warning("park_name not found in all_parks; continuing without it (not strictly required for CV).")
}

# ---- 1) Model formula (park-level PMR LM) ----
model_formula <- prs_coef ~ circle + Nota_manej + shannon_land + contig_mn_tree

# ---- 2) the CV function (unchanged) ----
cv_function <- function(train_indices, test_indices, data, formula) {
  train_data <- data[train_indices, ]
  test_data  <- data[test_indices, ]
  fit <- lm(formula, data = train_data, na.action = na.exclude)
  
  test_data$predicted_prs <- predict(fit, newdata = test_data)
  mse <- mean((test_data$prs_coef - test_data$predicted_prs)^2, na.rm = TRUE)
  
  rsq_in <- summary(fit)$r.squared # training R² (note: not a CV R²)
  
  list(
    MSE = mse,
    R_squared = rsq_in,
    Confidence_Intervals = confint(fit),
    SE = summary(fit)$coefficients[, "Std. Error"]
  )
}

# ---- 3) Run K-fold CV ----
set.seed(123)
n_folds <- 5L
folds <- cvTools::cvFolds(nrow(all_parks), K = n_folds)

cv_results <- lapply(seq_len(n_folds), function(i) {
  train_idx <- folds$subsets[folds$which != i]
  test_idx  <- folds$subsets[folds$which == i]
  cv_function(train_idx, test_idx, all_parks, model_formula)
})

# ---- 4) Summaries: fold table + overall means/sd ----
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

cat("\n--- PMR CV (park-level LM) — per-fold metrics ---\n")
print(fold_metrics, n=Inf)

cat("\n--- PMR CV (park-level LM) — overall summary ---\n")
print(overall_summary)

# ---- 5) Automatic evaluation (transparent, simple heuristics) ----
# Baseline variance of the target: how much error a "predict mean" model would have.
baseline_var <- var(all_parks$prs_coef, na.rm = TRUE)
mean_mse     <- overall_summary$Mean[overall_summary$Metric == "MSE (test)"]
sd_mse       <- overall_summary$SD[overall_summary$Metric == "MSE (test)"]
mean_r2_tr   <- overall_summary$Mean[overall_summary$Metric == "R² (training)"]
sd_r2_tr     <- overall_summary$SD[overall_summary$Metric == "R² (training)"]

ratio_to_baseline <- as.numeric(mean_mse / baseline_var)

cat("\n--- PMR CV — automatic interpretation ---\n")
cat(sprintf("Baseline variance of prs_coef (bigger = harder task): %.3f\n", baseline_var))
cat(sprintf("Mean MSE (test): %.3f  | SD across folds: %.3f  | Relative to baseline (Mean MSE / Var): %.2f\n",
            mean_mse, sd_mse, ratio_to_baseline))
cat(sprintf("Training R²: mean = %.3f  | SD = %.3f\n", mean_r2_tr, sd_r2_tr))

## =========================================================
## STEP 9 — PMR quick diagnostics
##  - Inputs: all_parks with prs_coef (19 parks filled) and
##            PMR park-level formula: prs_coef ~ circle + Nota_manej + shannon_land + contig_mn_tree
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
pmr_formula <- prs_coef ~ circle + Nota_manej + shannon_land + contig_mn_tree
pmr19 <- all_parks %>% filter(!is.na(prs_coef))
pmr_fit <- lm(pmr_formula, data = pmr19, na.action = na.exclude)

# Predictions & residuals
pmr19$pred <- predict(pmr_fit, newdata = pmr19)
pmr19$res  <- pmr19$prs_coef - pmr19$pred

# 2) (a) Observed vs Predicted
x11()
ggplot(pmr19, aes(x = pred, y = prs_coef, label = park_name)) +
  geom_point(size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = 50) +
  labs(x = "Predicted PMR (park signal)", y = "Observed PMR (prs_coef)",
       title = "PMR — Observed vs Predicted (19 parks)") +
  theme_bw()

# 2) (b) Rank agreement
pmr19 <- pmr19 %>%
  mutate(rank_obs = rank(prs_coef, ties.method = "average"),
         rank_pred = rank(pred,     ties.method = "average"))

rho_pmr <- suppressWarnings(cor(pmr19$rank_obs, pmr19$rank_pred, method = "spearman", use = "complete.obs"))

x11()
ggplot(pmr19, aes(x = rank_pred, y = rank_obs, label = park_name)) +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = 50) +
  labs(x = "Predicted rank", y = "Observed rank",
       title = paste0("PMR — Rank agreement (Spearman ρ = ", sprintf("%.2f", rho_pmr), ")")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

# 2) (c) Residuals vs Fitted
x11()
ggplot(pmr19, aes(x = pred, y = res, label = park_name)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = 50) +
  labs(x = "Fitted (predicted)", y = "Residual (obs - pred)",
       title = "PMR — Residuals vs Fitted") +
  theme_bw()

# ---- Spearman rank agreement (observed vs predicted) with p-value ----
sp_pmr <- suppressWarnings(cor.test(pmr19$prs_coef, pmr19$pred,
                                    method = "spearman", exact = FALSE))

rho_pmr <- unname(sp_pmr$estimate)
p_pmr   <- sp_pmr$p.value

cat(sprintf("\nPMR — Spearman rank correlation: ρ = %.3f (p = %.4f)\n",
            rho_pmr, p_pmr))

# PRINT
pmr_spearman_row <- data.frame(
  CES = "PMR",
  Spearman_rho = rho_pmr,
  p_value = p_pmr,
  stringsAsFactors = FALSE
)
print(pmr_spearman_row)

