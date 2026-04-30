############################################################
## AICc-based model selection for PMR and PAR (park-level)
## "Previous approach" replication:
## 
## - Correlation-based drop is applied (protected core vars never dropped)
## - Then all 4-predictor models from the resulting pool are tested
## - Filters: within-model |cor| <= cor_cutoff AND max VIF <= vif_cutoff
## - Exports AICc table + wide coefficients (estimates + p-values)
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(lme4)
  library(lmerTest)
  library(car)      # vif()
})

## =========================
## 0) USER SETTINGS
## =========================

# setwd("G:/Meu Drive/PESQUISA - Ecologia Urbana e Serviços Ecossistêmicos/Doutorado/Dados/CES_Parques/Statistical Analysis")

candidate_files <- c("CES_parks_final.csv", "CES_parks_Final.csv", "CES_parks_FINAL.csv")
infile <- candidate_files[file.exists(candidate_files)][1]
if (is.na(infile)) stop("Could not find CES parks file in the working directory. Edit candidate_files or setwd().")

out_prefix <- "AICc_selection_PREV"

cor_cutoff <- 0.60
vif_cutoff <- 5.00
model_size <- 4L

## =========================
## 1) HELPERS
## =========================

read_csv_br <- function(path) {
  read.csv(path, dec = ",", stringsAsFactors = FALSE)
}

scale_num <- function(x) {
  as.numeric(scale(as.numeric(x)))
}

# Greedy correlation-drop; never drops protected variables
drop_high_cor <- function(df, vars, cutoff = 0.60, protected = character(0)) {
  vars <- unique(vars)
  protected <- unique(protected)
  
  vars <- vars[vars %in% names(df)]
  if (length(vars) < 2) return(vars)
  
  safe_cor <- function(d, v) {
    X <- dplyr::select(d, dplyr::all_of(v))
    X <- as.data.frame(lapply(X, as.numeric))
    suppressWarnings(cor(X, use = "pairwise.complete.obs"))
  }
  
  keep <- vars
  cm <- safe_cor(df, keep)
  diag(cm) <- 0
  
  # Ensure protected are present
  if (length(setdiff(protected, keep)) > 0) {
    keep <- unique(c(keep, protected))
    cm <- safe_cor(df, keep)
    diag(cm) <- 0
  }
  
  while (max(abs(cm), na.rm = TRUE) > cutoff) {
    mean_abs <- colMeans(abs(cm), na.rm = TRUE)
    
    # Drop the "worst" variable that is NOT protected
    drop_candidates <- setdiff(names(sort(mean_abs, decreasing = TRUE)), protected)
    
    if (length(drop_candidates) == 0) {
      warning("Correlation > cutoff among protected vars; keeping them anyway.")
      break
    }
    
    drop_var <- drop_candidates[1]
    keep <- setdiff(keep, drop_var)
    
    if (length(keep) < 2) break
    cm <- safe_cor(df, keep)
    diag(cm) <- 0
  }
  
  keep
}

aicc_lm <- function(fit) {
  n <- stats::nobs(fit)
  k <- length(stats::coef(fit))  # includes intercept
  aic <- stats::AIC(fit)
  if ((n - k - 1) <= 0) return(NA_real_)
  aic + (2 * k * (k + 1)) / (n - k - 1)
}

akaike_weights <- function(aicc_vec) {
  dd <- aicc_vec - min(aicc_vec, na.rm = TRUE)
  rel <- exp(-0.5 * dd)
  rel / sum(rel, na.rm = TRUE)
}

coef_table_wide <- function(fit, model_id, candidate_vars) {
  sm <- summary(fit)
  coefs <- sm$coefficients
  terms <- rownames(coefs)
  
  out <- data.frame(model_id = model_id, stringsAsFactors = FALSE)
  
  out$`(Intercept)_est` <- if ("(Intercept)" %in% terms) coefs["(Intercept)", "Estimate"] else NA_real_
  out$`(Intercept)_p`   <- if ("(Intercept)" %in% terms) coefs["(Intercept)", "Pr(>|t|)"] else NA_real_
  
  for (v in candidate_vars) {
    est_col <- paste0(v, "_est")
    p_col   <- paste0(v, "_p")
    if (v %in% terms) {
      out[[est_col]] <- coefs[v, "Estimate"]
      out[[p_col]]   <- coefs[v, "Pr(>|t|)"]
    } else {
      out[[est_col]] <- NA_real_
      out[[p_col]]   <- NA_real_
    }
  }
  out
}

run_aicc_pool <- function(df, response, pool, model_size = 4L,
                          cor_cutoff = 0.60, vif_cutoff = 5.0, verbose = TRUE) {
  
  pool <- unique(pool)
  pool <- pool[pool %in% names(df)]
  if (length(pool) < model_size) stop("Pool smaller than model_size after checking names(df).")
  
  df[[response]] <- as.numeric(df[[response]])
  
  # Remove near-zero variance predictors
  var_ok <- vapply(pool, function(v) {
    x <- as.numeric(df[[v]])
    sd(x, na.rm = TRUE) > 1e-8
  }, logical(1))
  pool <- pool[var_ok]
  if (length(pool) < model_size) stop("After removing near-zero-variance variables, pool < model_size.")
  
  combos <- combn(pool, model_size, simplify = FALSE)
  
  ok_cor <- function(vars) {
    X <- df[, vars, drop = FALSE]
    X <- as.data.frame(lapply(X, as.numeric))
    cm <- suppressWarnings(cor(X, use = "pairwise.complete.obs"))
    diag(cm) <- 0
    max(abs(cm), na.rm = TRUE) <= cor_cutoff
  }
  
  results <- list()
  n_tested <- 0L
  n_kept <- 0L
  
  for (i in seq_along(combos)) {
    vars <- combos[[i]]
    n_tested <- n_tested + 1L
    
    if (!ok_cor(vars)) next
    
    f <- as.formula(paste(response, "~", paste(vars, collapse = " + ")))
    fit <- lm(f, data = df, na.action = na.exclude)
    
    v <- tryCatch(car::vif(fit), error = function(e) NA)
    if (any(is.na(v))) next
    if (max(v, na.rm = TRUE) > vif_cutoff) next
    
    aic  <- AIC(fit)
    aicc <- aicc_lm(fit)
    if (is.na(aicc)) next
    
    n_kept <- n_kept + 1L
    model_id <- paste0("m", sprintf("%03d", n_kept))
    
    results[[model_id]] <- list(
      model_id   = model_id,
      predictors = paste(vars, collapse = " + "),
      fit        = fit,
      AIC        = aic,
      AICc       = aicc,
      max_vif    = max(v, na.rm = TRUE)
    )
  }
  
  if (verbose) {
    cat("\nPool size:", length(pool),
        "| combinations tested:", n_tested,
        "| admissible models kept:", n_kept, "\n")
  }
  
  if (length(results) == 0) stop("No admissible models survived cor/VIF filters.")
  
  rank_df <- do.call(rbind, lapply(results, function(x) {
    data.frame(
      model_id   = x$model_id,
      predictors = x$predictors,
      AIC        = x$AIC,
      AICc       = x$AICc,
      max_vif    = x$max_vif,
      stringsAsFactors = FALSE
    )
  })) |> arrange(AICc)
  
  rank_df$delta_AICc <- rank_df$AICc - min(rank_df$AICc, na.rm = TRUE)
  rank_df$weight     <- akaike_weights(rank_df$AICc)
  
  list(rank_df = rank_df, models = results, pool = pool)
}

## =========================
## 2) LOAD + TRANSFORM DATA
## =========================

data_CES <- read_csv_br(infile)
data_CES$park <- as.factor(data_CES$pq_nome)

# People-level controls
data_CES$age_num <- as.numeric(data_CES$age)

data_CES$income_p <- ifelse(data_CES$income == "10 to 20 MW", 15,
                            ifelse(data_CES$income == "2 to 4 MW", 3,
                                   ifelse(data_CES$income == "> 20 MW", 20,
                                          ifelse(data_CES$income == "< 2 MW", 2,
                                                 ifelse(data_CES$income == "4 to 10 MW", 7, NA)))))

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

data_CES$age_s       <- scale_num(data_CES$age_num)
data_CES$income_p_s  <- scale_num(data_CES$income_p)
data_CES$nature_dose <- scale_num(as.numeric(data_CES$frequency_p) * as.numeric(data_CES$duration_p))

# Park-level predictors used in the "previous approach"
data_CES$circle_s       <- scale_num(data_CES$circle)
data_CES$Nota_manej_s   <- scale_num(data_CES$Nota_manej)
data_CES$Equip_espo_s   <- scale_num(data_CES$Equip_espo)
data_CES$Amb_aqua_s     <- scale_num(data_CES$Ambientes_)        # aquatic environments quality
data_CES$shannon_s      <- scale_num(data_CES$shannon_index)
data_CES$contig_tree_s  <- scale_num(data_CES$contig_mn_tree)

data_CES$area_s         <- scale_num(data_CES$pq_area)
data_CES$pland_tree_s   <- scale_num(as.numeric(data_CES$X1_pland_0m) + as.numeric(data_CES$X4_pland_0m))

data_CES$pland6_s       <- scale_num(data_CES$X6_pland_0m)
data_CES$np6_s          <- scale_num(data_CES$X6_np_0m)

data_CES$area_agua_s    <- scale_num(as.numeric(data_CES$area_agua) / as.numeric(data_CES$pq_area))
data_CES$presence_agua  <- ifelse(as.numeric(data_CES$area_agua) == 0, 0, 1)

park_predictors <- data_CES |>
  dplyr::select(
    park,
    circle_s, Nota_manej_s, shannon_s, contig_tree_s,
    Equip_espo_s, np6_s, Amb_aqua_s,
    area_s, pland_tree_s, pland6_s, area_agua_s, presence_agua
  ) |>
  dplyr::distinct(park, .keep_all = TRUE)

cat("\nParks in park_predictors:", n_distinct(park_predictors$park), "\n")

## =========================
## 3) PMR: mixed model -> park signal -> AICc
## =========================

data_CES$prs <- as.numeric(data_CES$prs_score)

pmr_people <- data_CES |>
  dplyr::select(park, prs, age_s, income_p_s, nature_dose) |>
  dplyr::filter(!is.na(prs))

m1_PMR <- lmer(prs ~ age_s + income_p_s + nature_dose + (1 | park), data = pmr_people)

cat("\n================ PMR mixed model summary ================\n")
print(summary(m1_PMR))

ri_PMR <- ranef(m1_PMR)$park
ri_PMR_df <- data.frame(
  park = as.factor(rownames(ri_PMR)),
  random_intercept = as.numeric(ri_PMR[, 1]),
  stringsAsFactors = FALSE
)

merged_PMR <- park_predictors |>
  left_join(ri_PMR_df, by = "park") |>
  dplyr::filter(!is.na(random_intercept))

# PMR pool (previous approach)
pmr_pool_raw <- c("circle_s", "Nota_manej_s", "shannon_s", "contig_tree_s", "area_s", "pland_tree_s")
pmr_protected <- c("circle_s", "Nota_manej_s", "shannon_s", "contig_tree_s")

pmr_pool <- drop_high_cor(merged_PMR, pmr_pool_raw, cutoff = cor_cutoff, protected = pmr_protected)

cat("\nPMR pool after correlation drop:\n")
print(pmr_pool)

pmr_sel  <- run_aicc_pool(merged_PMR, "random_intercept", pmr_pool,
                          model_size = model_size, cor_cutoff = cor_cutoff, vif_cutoff = vif_cutoff, verbose = TRUE)
pmr_rank <- pmr_sel$rank_df

cat("\n================ PMR AICc ranking (top 15) ================\n")
print(head(pmr_rank, 15))

pmr_best_id  <- pmr_rank$model_id[1]
pmr_best_fit <- pmr_sel$models[[pmr_best_id]]$fit

cat("\n================ PMR best model summary ================\n")
print(summary(pmr_best_fit))
cat("\nPMR best model formula:\n")
print(formula(pmr_best_fit))
cat("\nPMR best model VIF:\n")
print(car::vif(pmr_best_fit))

pmr_coef_wide <- do.call(rbind, lapply(pmr_rank$model_id, function(mid) {
  coef_table_wide(pmr_sel$models[[mid]]$fit, mid, pmr_sel$pool)
}))
pmr_out <- pmr_rank |> left_join(pmr_coef_wide, by = "model_id")

write.csv(pmr_rank, paste0(out_prefix, "_PMR_AICc_ranking.csv"), row.names = FALSE)
write.csv(pmr_out,  paste0(out_prefix, "_PMR_models_with_coefs.csv"), row.names = FALSE)

## =========================
## 4) PAR: mixed model -> park signal -> AICc
## =========================

data_CES$recreation <- as.numeric(data_CES$act_active)

par_people <- data_CES |>
  dplyr::select(park, recreation, age_s, income_p_s, nature_dose) |>
  dplyr::filter(!is.na(recreation))

m1_PAR <- glmer(recreation ~ age_s + income_p_s + nature_dose + (1 | park),
                data = par_people, family = binomial(link = "logit"))

cat("\n================ PAR mixed model summary ================\n")
print(summary(m1_PAR))

ri_PAR <- ranef(m1_PAR)$park
ri_PAR_df <- data.frame(
  park = as.factor(rownames(ri_PAR)),
  random_intercept = as.numeric(ri_PAR[, 1]),
  stringsAsFactors = FALSE
)

merged_PAR <- park_predictors |>
  left_join(ri_PAR_df, by = "park") |>
  dplyr::filter(!is.na(random_intercept))

# PAR pool (previous approach; IMPORTANT: no area_s, no accessibility/mobility/infra)
par_pool_raw <- c(
  "circle_s", "Equip_espo_s", "np6_s", "Amb_aqua_s",
  "pland6_s", "area_agua_s", "presence_agua"
)
par_protected <- c("circle_s", "Equip_espo_s", "np6_s", "Amb_aqua_s")

par_pool <- drop_high_cor(merged_PAR, par_pool_raw, cutoff = cor_cutoff, protected = par_protected)

cat("\nPAR pool after correlation drop:\n")
print(par_pool)

par_sel  <- run_aicc_pool(merged_PAR, "random_intercept", par_pool,
                          model_size = model_size, cor_cutoff = cor_cutoff, vif_cutoff = vif_cutoff, verbose = TRUE)
par_rank <- par_sel$rank_df

cat("\n================ PAR AICc ranking (top 15) ================\n")
print(head(par_rank, 15))

par_best_id  <- par_rank$model_id[1]
par_best_fit <- par_sel$models[[par_best_id]]$fit

cat("\n================ PAR best model summary ================\n")
print(summary(par_best_fit))
cat("\nPAR best model formula:\n")
print(formula(par_best_fit))
cat("\nPAR best model VIF:\n")
print(car::vif(par_best_fit))

par_coef_wide <- do.call(rbind, lapply(par_rank$model_id, function(mid) {
  coef_table_wide(par_sel$models[[mid]]$fit, mid, par_sel$pool)
}))
par_out <- par_rank |> left_join(par_coef_wide, by = "model_id")

write.csv(par_rank, paste0(out_prefix, "_PAR_AICc_ranking.csv"), row.names = FALSE)
write.csv(par_out,  paste0(out_prefix, "_PAR_models_with_coefs.csv"), row.names = FALSE)

## =========================
## 5) SAVE OBJECTS
## =========================

saveRDS(list(
  data_file = infile,
  settings = list(cor_cutoff = cor_cutoff, vif_cutoff = vif_cutoff, model_size = model_size),
  PMR = list(mixed_model = m1_PMR, park_data = merged_PMR, pool = pmr_sel$pool, aicc_ranking = pmr_rank,
             best_model_id = pmr_best_id, best_model = pmr_best_fit),
  PAR = list(mixed_model = m1_PAR, park_data = merged_PAR, pool = par_sel$pool, aicc_ranking = par_rank,
             best_model_id = par_best_id, best_model = par_best_fit)
), file = paste0(out_prefix, "_objects.rds"))

cat("\nSaved RDS:", paste0(out_prefix, "_objects.rds"), "\n")

getwd()


#############

## =========================
## 8) COVERAGE AUDIT: what exists vs what was tested
## =========================

# Split "a + b + c" predictor strings
split_predictors <- function(x) trimws(unlist(strsplit(x, "\\+")))

# Variables that actually appeared in any admissible model (post filters)
vars_in_rank <- function(rank_df) unique(unlist(lapply(rank_df$predictors, split_predictors)))

pmr_used_any <- vars_in_rank(pmr_rank)
par_used_any <- vars_in_rank(par_rank)

# Pools you *intended* to test (pre-filter)
pmr_pool_full <- pmr_sel$pool
par_pool_full <- par_sel$pool

# All candidate columns you created at park-level
park_cols <- setdiff(names(park_predictors), "park")

# Helper to compute SD / NA counts on the park-level table
park_var_stats <- function(df, vars) {
  data.frame(
    var = vars,
    n_total = nrow(df),
    n_na = vapply(vars, function(v) sum(is.na(df[[v]])), integer(1)),
    sd = vapply(vars, function(v) sd(as.numeric(df[[v]]), na.rm = TRUE), numeric(1)),
    stringsAsFactors = FALSE
  )
}

# 1) Summary of what was in park_predictors vs what was pooled/tested
audit_core <- data.frame(
  var = sort(unique(park_cols)),
  in_park_predictors = TRUE,
  in_pmr_pool = sort(unique(park_cols)) %in% pmr_pool_full,
  in_par_pool = sort(unique(park_cols)) %in% par_pool_full,
  appeared_in_any_PMR_model = sort(unique(park_cols)) %in% pmr_used_any,
  appeared_in_any_PAR_model = sort(unique(park_cols)) %in% par_used_any,
  stringsAsFactors = FALSE
) |>
  dplyr::left_join(park_var_stats(park_predictors, sort(unique(park_cols))), by = "var") |>
  dplyr::arrange(dplyr::desc(in_pmr_pool | in_par_pool), dplyr::desc(appeared_in_any_PMR_model | appeared_in_any_PAR_model), var)

cat("\n================ COVERAGE AUDIT (PARK PREDICTORS) ================\n")
print(audit_core, row.names = FALSE)

write.csv(audit_core, paste0(out_prefix, "_coverage_audit_park_predictors.csv"), row.names = FALSE)

# 2) Search raw data columns for “missing from analysis” candidates by keywords/patterns
#    (this helps you see what exists but you never transformed into *_s / never added to park_predictors)
raw_names <- names(data_CES)

patterns <- c(
  "lpi", "ed_", "ed$", "pland", "np_", "np$", "nl", "nlsi", "contig", "shannon", "edge",
  "altura", "height", "volume", "bambu", "bamboo",
  "crianca", "children",
  "sinal", "sign", "sinaliz",
  "segur", "security",
  "adm", "admin",
  "infra", "manej", "mobil", "acess", "equip", "aqua", "agua", "water"
)

hits <- unique(unlist(lapply(patterns, function(p) grep(p, raw_names, ignore.case = TRUE, value = TRUE))))

raw_inventory <- data.frame(
  raw_col = sort(hits),
  already_in_park_predictors = sort(hits) %in% names(park_predictors),
  stringsAsFactors = FALSE
)

cat("\n================ RAW COLUMNS THAT LOOK LIKE CANDIDATES ================\n")
print(raw_inventory, row.names = FALSE)

write.csv(raw_inventory, paste0(out_prefix, "_raw_candidate_columns_inventory.csv"), row.names = FALSE)

cat("\nSaved:\n",
    paste0(out_prefix, "_coverage_audit_park_predictors.csv\n"),
    paste0(out_prefix, "_raw_candidate_columns_inventory.csv\n"),
    sep = "")
