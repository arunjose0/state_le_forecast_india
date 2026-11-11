# Temporal holdout validation for the life-expectancy model.
#
# Origins are defined by the first test year:
#   Start 2000: train on observations before 2000; score observed years 2000-2019.
#   Start 2010: train on observations before 2010; score observed years 2010-2019.
#
# 2020 onward is excluded from scoring because COVID shocks cannot be forecast
# ex ante from a training sample ending in 1999 or 2009. All forecasts are
# genuinely out of sample: held-out life-expectancy values are not passed to Stan.
# Only the selected main/original prior is refit in holdout validation. Simple
# persistence and linear-trend benchmarks remain for predictive context; the
# tighter and wider priors are assessed separately in prior-sensitivity analysis.

suppressPackageStartupMessages({
  library(here)
  library(rstan)
  library(tidyverse)
})

rstan_options(auto_write = TRUE)
options(mc.cores = min(4L, parallel::detectCores(logical = TRUE)))

# ------------------------------ Settings --------------------------------------
FORECAST_START_YEARS <- c(2000L, 2010L)
HOLDOUT_END_YEAR <- 2019L
PRIOR_MODES <- 2L

ITER <- 2000L
WARMUP <- 1000L
CHAINS <- 4L
ADAPT_DELTA <- 0.92
MAX_TREEDEPTH <- 12L
SEED_BASE <- 52100L
STOP_AFTER_FAILED_SCREEN <- TRUE
REUSE_EXISTING_FITS <- TRUE
FIT_TAG <- "final1"

PRIOR_LABELS <- c(
  `2` = "Main model (regularised baseline)",
  `3` = "Tighter-prior sensitivity",
  `4` = "Wide-prior sensitivity"
)
MODEL_COLOURS <- c(
  "Main model (regularised baseline)" = "#D55E00",
  "Tighter-prior sensitivity" = "#0072B2",
  "Wide-prior sensitivity" = "#009E73",
  # Distinct colour-blind-friendly colours keep both benchmark methods visible.
  "Benchmark: persistence" = "#0072B2",
  "Benchmark: linear trend" = "#009E73"
)

# --------------------------- Project paths -----------------------------------
prep_path <- here::here("fits", "prep_objects.rds")
stan_path <- here::here("stan", "life_expectancy_model.stan")
output_tables <- here::here("tables", "holdout_validation")
output_figures <- here::here("figures", "holdout_validation")

if (!file.exists(prep_path)) stop("Missing input: ", prep_path, call. = FALSE)
if (!file.exists(stan_path)) stop("Missing unified Stan file: ", stan_path, call. = FALSE)
for (path in c(output_tables, output_figures, here::here("fits"))) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

# -------------------------- Load and check data -------------------------------
prep <- readRDS(prep_path)
df_input <- tibble::as_tibble(prep$df)
stan_full <- prep$stan_data_full
state_labels <- as.character(prep$state_levels_full)

stopifnot(nrow(df_input) == stan_full$N)
stopifnot(length(state_labels) == stan_full$N_groups)
stopifnot(all(c("state", "year", "g", "e0_f", "e0_m") %in% names(df_input)))
if (length(stan_full$eag_group) != length(state_labels)) {
  stop("Expected stan_data_full$eag_group to contain one value per state.", call. = FALSE)
}

# This is the same exploratory tier allocation used by the fitting workflow.
tier_levels <- c("Front_runner", "Achiever", "Aspirant")
state_tier_map <- tibble::tribble(
  ~state_name,             ~tier,
  "Andhra Pradesh",        "Front_runner",
  "Telangana",             "Front_runner",
  "Kerala",                "Front_runner",
  "Maharashtra",           "Front_runner",
  "Tamil Nadu",            "Front_runner",
  "West Bengal",           "Front_runner",
  "Chhattisgarh",          "Achiever",
  "Goa",                   "Achiever",
  "Gujarat",               "Achiever",
  "Haryana",               "Achiever",
  "Himachal Pradesh",      "Achiever",
  "Karnataka",             "Achiever",
  "Punjab",                "Achiever",
  "Assam",                 "Aspirant",
  "Bihar",                 "Aspirant",
  "Delhi",                 "Aspirant",
  "Jammu & Kashmir",       "Aspirant",
  "Jharkhand",             "Aspirant",
  "Madhya Pradesh",        "Aspirant",
  "Odisha",                "Aspirant",
  "Rajasthan",             "Aspirant",
  "Uttar Pradesh",         "Aspirant",
  "Uttarakhand",           "Aspirant"
)
normalise_state <- function(x) {
  key <- tolower(gsub("[^[:alnum:]]", "", as.character(x)))
  dplyr::case_when(
    key == "chhatisgarh" ~ "chhattisgarh",
    key == "orissa" ~ "odisha",
    TRUE ~ key
  )
}
state_tier_map <- bind_rows(
  state_tier_map,
  tibble(state_name = c("Chhatisgarh", "Orissa"), tier = c("Achiever", "Aspirant"))
) %>%
  mutate(state_key = normalise_state(state_name)) %>%
  distinct(state_key, .keep_all = TRUE)

tier_match <- match(normalise_state(state_labels), state_tier_map$state_key)
if (anyNA(tier_match)) {
  stop("No Health Index tier mapped for: ",
       paste(state_labels[is.na(tier_match)], collapse = ", "), call. = FALSE)
}
tier_by_state <- match(state_tier_map$tier[tier_match], tier_levels)
if (anyNA(tier_by_state)) stop("Unrecognised tier ID.", call. = FALSE)

# Attach flags before sorting/filtering so every row remains aligned to prep data.
df_full <- df_input %>%
  mutate(
    .old_g = as.integer(g),
    .is_2020 = as.integer(stan_full$is_2020),
    .is_2021 = as.integer(stan_full$is_2021),
    .is_post2021 = as.integer(stan_full$is_post2021),
    .years_since_2021 = as.numeric(stan_full$years_since_2021)
  ) %>%
  filter(is.finite(year), is.finite(e0_f), is.finite(e0_m), year <= HOLDOUT_END_YEAR) %>%
  arrange(.old_g, year)

if (anyDuplicated(df_full[c(".old_g", "year")])) {
  stop("Duplicate state-year rows detected. Resolve before holdout validation.", call. = FALSE)
}

# ------------------------ Build temporal holdout fold -------------------------
make_fold <- function(start_year) {
  train_source <- df_full %>% filter(year < start_year)
  test_source <- df_full %>% filter(year >= start_year, year <= HOLDOUT_END_YEAR)

  if (nrow(train_source) == 0L || nrow(test_source) == 0L) {
    stop("No training or test observations for start year ", start_year, call. = FALSE)
  }

  train_state_ids <- sort(unique(train_source$.old_g))
  group_map <- tibble(old_g = train_state_ids) %>%
    mutate(
      g_new = row_number(),
      state = state_labels[old_g],
      tier_id = as.integer(tier_by_state[old_g]),
      eag_group = as.integer(stan_full$eag_group[old_g])
    )

  # A state without any pre-origin observations cannot be forecast from this fit.
  test_excluded <- test_source %>%
    anti_join(group_map %>% select(old_g), by = c(".old_g" = "old_g")) %>%
    count(state, name = "n_test_rows_excluded") %>%
    mutate(forecast_start_year = start_year, .before = 1)

  train <- train_source %>%
    inner_join(group_map %>% select(old_g, g_new), by = c(".old_g" = "old_g")) %>%
    arrange(g_new, year) %>%
    mutate(local_idx = row_number()) %>%
    group_by(g_new) %>%
    mutate(
      row_in_group = row_number(),
      is_first = as.integer(row_in_group == 1L),
      prev_idx = if_else(row_in_group == 1L, 0L, as.integer(lag(local_idx))),
      delta_t = if_else(row_in_group == 1L, 0, as.numeric(year - lag(year)))
    ) %>%
    ungroup()

  test <- test_source %>%
    inner_join(group_map %>% select(old_g, g_new), by = c(".old_g" = "old_g")) %>%
    arrange(g_new, year) %>%
    mutate(state = as.character(state_labels[.old_g]))

  if (nrow(test) == 0L) stop("No test rows belong to a state with training history for ", start_year)
  if (any(train$delta_t[train$is_first == 0L] <= 0)) stop("Non-positive training time gap.")

  by_group <- train %>%
    group_by(g_new) %>%
    summarise(
      first_year = first(year), last_year = last(year),
      first_idx = first(local_idx), last_idx = last(local_idx),
      .groups = "drop"
    ) %>%
    arrange(g_new)

  forecast_years <- sort(unique(as.numeric(test$year)))
  stan_data_base <- list(
    N = as.integer(nrow(train)),
    g = as.integer(train$g_new),
    is_first = as.integer(train$is_first),
    prev_idx = as.integer(train$prev_idx),
    delta_t = as.numeric(train$delta_t),
    y_f = as.numeric(train$e0_f),
    y_m = as.numeric(train$e0_m),
    N_groups = as.integer(nrow(group_map)),
    N_tiers = as.integer(length(tier_levels)),
    tier_id = as.integer(group_map$tier_id),
    is_2020 = as.integer(train$.is_2020),
    is_2021 = as.integer(train$.is_2021),
    is_post2021 = as.integer(train$.is_post2021),
    years_since_2021 = as.numeric(train$.years_since_2021),
    N_forecast = as.integer(length(forecast_years)),
    forecast_years = forecast_years,
    first_year = as.numeric(by_group$first_year),
    last_year = as.numeric(by_group$last_year),
    last_idx = as.integer(by_group$last_idx),
    prior_mode = 2L,
    prior_only = 0L,
    eag_group = as.integer(group_map$eag_group),
    compute_yrep = 0L
  )

  # Basic index and chronology checks before Stan receives the fold.
  stopifnot(stan_data_base$N == nrow(train))
  stopifnot(sum(stan_data_base$is_first) == stan_data_base$N_groups)
  stopifnot(all(stan_data_base$g >= 1L & stan_data_base$g <= stan_data_base$N_groups))
  stopifnot(all(stan_data_base$last_year < start_year))
  stopifnot(all(test$year >= start_year & test$year <= HOLDOUT_END_YEAR))
  stopifnot(all(stan_data_base$last_idx >= 1L & stan_data_base$last_idx <= stan_data_base$N))

  list(
    start_year = start_year,
    end_year = HOLDOUT_END_YEAR,
    train = train,
    test = test,
    group_map = group_map,
    by_group = by_group,
    forecast_years = forecast_years,
    stan_data = stan_data_base,
    excluded_states = test_excluded
  )
}

# ----------------------------- Metrics ----------------------------------------
crps_ensemble <- function(draws, observed) {
  x <- sort(as.numeric(draws))
  n <- length(x)
  if (n == 0L || !is.finite(observed) || !all(is.finite(x))) return(NA_real_)
  first_term <- mean(abs(x - observed))
  second_term <- sum((2 * seq_len(n) - n - 1) * x) / (n^2)
  max(0, first_term - second_term)
}

summarise_predictive_draws <- function(draws, observed) {
  draws <- as.numeric(draws)
  if (!length(draws) || !all(is.finite(draws))) {
    stop("Non-finite or empty predictive draws; refusing to score this fold.", call. = FALSE)
  }
  q <- quantile(draws, probs = c(.025, .05, .25, .50, .75, .95, .975), names = FALSE, type = 8)
  tibble(
    observed = observed,
    pred_q025 = q[1], pred_q05 = q[2], pred_q25 = q[3], pred_median = q[4],
    pred_q75 = q[5], pred_q95 = q[6], pred_q975 = q[7],
    width50 = q[5] - q[3], width90 = q[6] - q[2], width95 = q[7] - q[1],
    crps = crps_ensemble(draws, observed),
    cover50 = as.integer(observed >= q[3] && observed <= q[5]),
    cover90 = as.integer(observed >= q[2] && observed <= q[6]),
    cover95 = as.integer(observed >= q[1] && observed <= q[7])
  )
}

benchmark_predictions <- function(history, target_year, response) {
  history <- history %>% arrange(year)
  persistence <- dplyr::last(history[[response]])
  linear <- persistence
  if (nrow(history) >= 2L && n_distinct(history$year) >= 2L) {
    lm_fit <- tryCatch(lm(reformulate("year", response = response), data = history), error = function(e) NULL)
    if (!is.null(lm_fit)) {
      pred <- tryCatch(as.numeric(predict(lm_fit, newdata = tibble(year = target_year))), error = function(e) NA_real_)
      if (length(pred) == 1L && is.finite(pred)) linear <- pred
    }
  }
  c(persistence = persistence, linear_trend = linear)
}

# Check sampled parameters only; generated forecast arrays are not convergence parameters.
screen_fit <- function(fit, start_year, prior_mode) {
  summary_df <- as.data.frame(rstan::summary(fit)$summary) %>%
    tibble::rownames_to_column("parameter")
  omit_pattern <- paste0(
    "^(y_rep_f|y_rep_m|mu_curve_f|mu_curve_m|C_tier_f|C_tier_m|tier_effect_f|tier_effect_m|",
    "elapsed_since_first|baseline_f|baseline_m|covid_shock_f|covid_shock_m|",
    "trend_forecast_f|trend_forecast_m|mean_forecast_f|mean_forecast_m|",
    "latent_forecast_f|latent_forecast_m|e0_forecast_f|e0_forecast_m|",
    "u_pred_mean|u_pred_cov|u_filter_mean|u_filter_cov|marginal_loglik|",
    "I2|L_scaled|L_Omega_dev|L_Omega_C|Omega_dev|Omega_C|L_C_scaled)(\\[|$)"
  )
  diag_pars <- summary_df %>%
    filter(parameter != "lp__") %>%
    filter(!grepl(omit_pattern, parameter) | parameter %in% c("L_Omega_dev[2,1]", "L_Omega_C[2,1]"))
  bfmi <- as.numeric(rstan::get_bfmi(fit))
  d <- tibble(
    forecast_start_year = start_year,
    prior_mode = prior_mode,
    model = unname(PRIOR_LABELS[as.character(prior_mode)]),
    divergences = as.integer(rstan::get_num_divergent(fit)),
    max_treedepth_hits = as.integer(rstan::get_num_max_treedepth(fit)),
    min_bfmi = min(bfmi, na.rm = TRUE),
    chains_bfmi_below_0_3 = sum(bfmi < 0.30, na.rm = TRUE),
    max_rhat = max(diag_pars$Rhat[is.finite(diag_pars$Rhat)], na.rm = TRUE),
    min_ess = min(diag_pars$n_eff[is.finite(diag_pars$n_eff)], na.rm = TRUE),
    nonfinite_rhat = sum(!is.finite(diag_pars$Rhat)),
    nonfinite_ess = sum(!is.finite(diag_pars$n_eff))
  ) %>% mutate(
    check_zero_divergences = divergences == 0L,
    check_no_treedepth_hits = max_treedepth_hits == 0L,
    check_bfmi = is.finite(min_bfmi) && min_bfmi >= 0.30 && chains_bfmi_below_0_3 == 0L,
    check_rhat = is.finite(max_rhat) && max_rhat <= 1.01,
    check_ess = is.finite(min_ess) && min_ess >= 100,
    check_finite_diagnostics = nonfinite_rhat == 0L && nonfinite_ess == 0L
  )
  check_cols <- c("check_zero_divergences", "check_no_treedepth_hits", "check_bfmi",
                  "check_rhat", "check_ess", "check_finite_diagnostics")
  passed <- all(unlist(d[1, check_cols], use.names = FALSE))
  d$sampler_status <- if (passed) "PASS" else "FAIL"
  failed <- check_cols[!as.logical(d[1, check_cols])]
  d$failed_checks <- if (!length(failed)) "None" else paste(failed, collapse = "; ")
  cat("Holdout fit ", start_year, " / prior mode ", prior_mode, ": ", d$sampler_status[1],
      if (passed) "" else paste0(" (", d$failed_checks[1], ")"), "\n", sep = "")
  list(summary = d, pass = passed)
}

extract_model_forecasts <- function(fit, fold, prior_mode) {
  ex <- rstan::extract(fit, pars = c(
    "latent_forecast_f", "latent_forecast_m", "e0_forecast_f", "e0_forecast_m"
  ), permuted = TRUE)
  model_label <- unname(PRIOR_LABELS[as.character(prior_mode)])
  result <- list()
  k <- 0L

  # Summarise full forecast grids for plots, plus observed test rows for scores.
  for (type in c("latent", "observed")) {
    for (sex in c("Female", "Male")) {
      suffix <- if (sex == "Female") "f" else "m"
      variable <- if (type == "latent") paste0("latent_forecast_", suffix) else paste0("e0_forecast_", suffix)
      arr <- ex[[variable]]
      if (is.null(arr) || length(dim(arr)) != 3L ||
          dim(arr)[2] != nrow(fold$group_map) || dim(arr)[3] != length(fold$forecast_years)) {
        stop("Unexpected dimensions for ", variable, " in start year ", fold$start_year, call. = FALSE)
      }
      for (gi in seq_len(nrow(fold$group_map))) {
        for (ti in seq_along(fold$forecast_years)) {
          draws <- as.numeric(arr[, gi, ti])
          if (!all(is.finite(draws))) stop("Non-finite forecasts in ", variable, call. = FALSE)
          q <- quantile(draws, probs = c(.025, .05, .25, .50, .75, .95, .975), names = FALSE, type = 8)
          k <- k + 1L
          result[[k]] <- tibble(
            forecast_start_year = fold$start_year,
            prior_mode = prior_mode,
            model = model_label,
            forecast_type = type,
            sex = sex,
            state = fold$group_map$state[gi],
            g_new = gi,
            year = fold$forecast_years[ti],
            q025 = q[1], q05 = q[2], q25 = q[3], median = q[4],
            q75 = q[5], q95 = q[6], q975 = q[7],
            width50 = q[5] - q[3], width90 = q[6] - q[2], width95 = q[7] - q[1]
          )
        }
      }
    }
  }
  grid_summary <- bind_rows(result)

  score_rows <- list()
  s <- 0L
  for (i in seq_len(nrow(fold$test))) {
    gi <- as.integer(fold$test$g_new[i])
    ti <- match(as.numeric(fold$test$year[i]), fold$forecast_years)
    if (is.na(ti)) stop("Held-out year absent from forecast grid.", call. = FALSE)
    history <- fold$train %>% filter(g_new == gi) %>% arrange(year)

    for (sex in c("Female", "Male")) {
      suffix <- if (sex == "Female") "f" else "m"
      draws <- as.numeric(ex[[paste0("e0_forecast_", suffix)]][, gi, ti])
      observed <- if (sex == "Female") fold$test$e0_f[i] else fold$test$e0_m[i]
      sm <- summarise_predictive_draws(draws, observed)
      response <- if (sex == "Female") "e0_f" else "e0_m"
      bench <- benchmark_predictions(history, fold$test$year[i], response)
      last_train <- fold$by_group$last_year[match(gi, fold$by_group$g_new)]

      # Benchmarks do not depend on priors; save them once per test row from mode 2.
      if (prior_mode == 2L) {
        for (bench_name in names(bench)) {
          s <- s + 1L
          score_rows[[s]] <- tibble(
            forecast_start_year = fold$start_year,
            forecast_end_year = fold$end_year,
            model = paste0("Benchmark: ", if (bench_name == "linear_trend") "linear trend" else "persistence"),
            prior_mode = NA_integer_,
            state = as.character(fold$test$state[i]),
            g_new = gi,
            year = as.numeric(fold$test$year[i]),
            last_train_year = last_train,
            lead_years = as.numeric(fold$test$year[i] - last_train),
            sex = sex,
            observed = observed,
            pred_median = unname(bench[[bench_name]]),
            pred_q025 = NA_real_, pred_q05 = NA_real_, pred_q25 = NA_real_,
            pred_q75 = NA_real_, pred_q95 = NA_real_, pred_q975 = NA_real_,
            width50 = NA_real_, width90 = NA_real_, width95 = NA_real_,
            crps = NA_real_, cover50 = NA_real_, cover90 = NA_real_, cover95 = NA_real_
          ) %>% mutate(error = pred_median - observed, abs_error = abs(error), sq_error = error^2)
        }
      }

      s <- s + 1L
      score_rows[[s]] <- tibble(
        forecast_start_year = fold$start_year,
        forecast_end_year = fold$end_year,
        model = model_label,
        prior_mode = prior_mode,
        state = as.character(fold$test$state[i]),
        g_new = gi,
        year = as.numeric(fold$test$year[i]),
        last_train_year = last_train,
        lead_years = as.numeric(fold$test$year[i] - last_train),
        sex = sex
      ) %>% bind_cols(sm) %>%
        mutate(error = pred_median - observed, abs_error = abs(error), sq_error = error^2)
    }
  }

  list(grid = grid_summary, scores = bind_rows(score_rows))
}

# ----------------------------- Fit and score ----------------------------------
cat("Compiling the single unified Stan model...\n")
stan_model <- rstan::stan_model(file = stan_path)
all_grids <- list()
all_scores <- list()
all_diagnostics <- list()
all_exclusions <- list()
grid_i <- score_i <- diagnostic_i <- exclusion_i <- 0L

for (start_year in FORECAST_START_YEARS) {
  fold <- make_fold(start_year)
  exclusion_i <- exclusion_i + 1L
  all_exclusions[[exclusion_i]] <- fold$excluded_states

  cat("\nForecast start: ", start_year,
      " | training rows: ", nrow(fold$train),
      " | training states: ", nrow(fold$group_map),
      " | held-out rows: ", nrow(fold$test),
      " | scoring through: ", HOLDOUT_END_YEAR, "\n", sep = "")

  for (prior_mode in PRIOR_MODES) {
    model_label <- unname(PRIOR_LABELS[as.character(prior_mode)])
    fit_path <- here::here("fits", sprintf(
      "holdout_start_%d_mode%d_%s_%s.rds",
      start_year, prior_mode,
      c(`2` = "main_regularised", `3` = "tighter_prior", `4` = "wide_prior")[as.character(prior_mode)],
      FIT_TAG
    ))

    if (REUSE_EXISTING_FITS && file.exists(fit_path)) {
      cat("Loading cached fit: ", fit_path, "\n", sep = "")
      fit <- readRDS(fit_path)
    } else {
      fold_data <- fold$stan_data
      fold_data$prior_mode <- as.integer(prior_mode)
      fit <- rstan::sampling(
        stan_model,
        data = fold_data,
        iter = ITER, warmup = WARMUP, chains = CHAINS,
        control = list(adapt_delta = ADAPT_DELTA, max_treedepth = MAX_TREEDEPTH),
        seed = SEED_BASE + start_year + prior_mode * 100L,
        refresh = 200
      )
      saveRDS(fit, fit_path)
    }

    diag <- screen_fit(fit, start_year, prior_mode)
    diagnostic_i <- diagnostic_i + 1L
    all_diagnostics[[diagnostic_i]] <- diag$summary
    write_csv(diag$summary, here::here("tables", "holdout_validation",
      sprintf("diagnostics_start_%d_mode_%d.csv", start_year, prior_mode)))
    write_csv(diag$summary %>% select(forecast_start_year, prior_mode, model, sampler_status,
      failed_checks, divergences, max_treedepth_hits, min_bfmi, chains_bfmi_below_0_3,
      max_rhat, min_ess), here::here("tables", "holdout_validation",
      sprintf("diagnostic_status_start_%d_mode_%d.csv", start_year, prior_mode)))

    if (!diag$pass) {
      message("Sampler screen failed for start year ", start_year, ", mode ", prior_mode,
              ". Fit and diagnostics have been saved.")
      if (STOP_AFTER_FAILED_SCREEN) stop("Do not score a failed fit; inspect its diagnostics.", call. = FALSE)
    }

    extracted <- extract_model_forecasts(fit, fold, prior_mode)
    grid_i <- grid_i + 1L
    all_grids[[grid_i]] <- extracted$grid
    score_i <- score_i + 1L
    all_scores[[score_i]] <- extracted$scores
  }
}

forecast_grid <- bind_rows(all_grids) %>%
  mutate(
    forecast_start_year = as.integer(forecast_start_year),
    sex = factor(sex, levels = c("Female", "Male")),
    forecast_type = factor(forecast_type, levels = c("latent", "observed"),
                           labels = c("Latent trajectory", "Future observation")),
    model = factor(model, levels = unname(PRIOR_LABELS[as.character(PRIOR_MODES)]))
  )
all_scores <- bind_rows(all_scores) %>%
  mutate(
    lead_band = case_when(
      lead_years <= 5 ~ "1-5 years",
      lead_years <= 10 ~ "6-10 years",
      TRUE ~ "11+ years"
    ),
    model = factor(model, levels = c(unname(PRIOR_LABELS[as.character(PRIOR_MODES)]),
      "Benchmark: persistence", "Benchmark: linear trend"))
  )

write_csv(forecast_grid, file.path(output_tables, "forecast_grid_main_model.csv"))
write_csv(all_scores, file.path(output_tables, "holdout_predictions_main_model_and_benchmarks.csv"))
write_csv(bind_rows(all_diagnostics), file.path(output_tables, "sampler_diagnostics_main_model_holdout_fits.csv"))
write_csv(bind_rows(all_exclusions), file.path(output_tables, "states_without_pre_origin_training_data.csv"))

# Summaries for probabilistic forecasts and point-forecast benchmarks.
metrics <- all_scores %>%
  group_by(forecast_start_year, model, sex, lead_band) %>%
  summarise(
    n_predictions = sum(is.finite(observed) & is.finite(pred_median)),
    MAE = mean(abs_error, na.rm = TRUE),
    MSE = mean(sq_error, na.rm = TRUE),
    RMSE = sqrt(mean(sq_error, na.rm = TRUE)),
    mean_bias = mean(error, na.rm = TRUE),
    CRPS = if (all(is.na(crps))) NA_real_ else mean(crps, na.rm = TRUE),
    coverage50 = if (all(is.na(cover50))) NA_real_ else mean(cover50, na.rm = TRUE),
    coverage90 = if (all(is.na(cover90))) NA_real_ else mean(cover90, na.rm = TRUE),
    coverage95 = if (all(is.na(cover95))) NA_real_ else mean(cover95, na.rm = TRUE),
    mean_width50 = if (all(is.na(width50))) NA_real_ else mean(width50, na.rm = TRUE),
    mean_width90 = if (all(is.na(width90))) NA_real_ else mean(width90, na.rm = TRUE),
    mean_width95 = if (all(is.na(width95))) NA_real_ else mean(width95, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(metrics, file.path(output_tables, "holdout_metrics_by_origin_sex_lead.csv"))

state_metrics <- all_scores %>%
  group_by(forecast_start_year, model, sex, state) %>%
  summarise(
    n_predictions = n(),
    MAE = mean(abs_error, na.rm = TRUE),
    MSE = mean(sq_error, na.rm = TRUE),
    RMSE = sqrt(mean(sq_error, na.rm = TRUE)),
    mean_bias = mean(error, na.rm = TRUE),
    CRPS = if (all(is.na(crps))) NA_real_ else mean(crps, na.rm = TRUE),
    coverage90 = if (all(is.na(cover90))) NA_real_ else mean(cover90, na.rm = TRUE),
    mean_width90 = if (all(is.na(width90))) NA_real_ else mean(width90, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(state_metrics, file.path(output_tables, "holdout_metrics_by_state.csv"))

# Main comparison table without lead bands, to support concise reporting.
overall_metrics <- all_scores %>%
  group_by(forecast_start_year, model, sex) %>%
  summarise(
    n_predictions = n(),
    MAE = mean(abs_error, na.rm = TRUE),
    MSE = mean(sq_error, na.rm = TRUE),
    RMSE = sqrt(mean(sq_error, na.rm = TRUE)),
    mean_bias = mean(error, na.rm = TRUE),
    CRPS = if (all(is.na(crps))) NA_real_ else mean(crps, na.rm = TRUE),
    coverage50 = if (all(is.na(cover50))) NA_real_ else mean(cover50, na.rm = TRUE),
    coverage90 = if (all(is.na(cover90))) NA_real_ else mean(cover90, na.rm = TRUE),
    coverage95 = if (all(is.na(cover95))) NA_real_ else mean(cover95, na.rm = TRUE),
    mean_width50 = if (all(is.na(width50))) NA_real_ else mean(width50, na.rm = TRUE),
    mean_width90 = if (all(is.na(width90))) NA_real_ else mean(width90, na.rm = TRUE),
    mean_width95 = if (all(is.na(width95))) NA_real_ else mean(width95, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(overall_metrics, file.path(output_tables, "holdout_metrics_overall.csv"))

# --------------------------- Publication figures ------------------------------
# Only two holdout figures are produced: observed versus predicted for the
# 2000-origin and 2010-origin folds. Metrics and diagnostic tables above are
# retained; interval, trajectory, and coverage figures are intentionally omitted.
base_theme <- theme_classic(base_size = 10) +
  theme(
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    plot.title = element_blank(),
    plot.subtitle = element_blank()
  )

plot_models <- c(
  "Main model (regularised baseline)",
  "Benchmark: persistence",
  "Benchmark: linear trend"
)
plot_colours <- MODEL_COLOURS[plot_models]
plot_shapes <- c(
  "Main model (regularised baseline)" = 16,
  "Benchmark: persistence" = 17,
  "Benchmark: linear trend" = 15
)

# Clean obsolete figures from earlier script versions. Keep the output directory
# limited to the two requested observed-versus-predicted plots.
expected_holdout_figures <- sprintf(
  "holdout_observed_vs_predicted_start_%d.png", FORECAST_START_YEARS
)
existing_pngs <- list.files(output_figures, pattern = "\\.png$", full.names = TRUE)
stale_pngs <- existing_pngs[!basename(existing_pngs) %in% expected_holdout_figures]
if (length(stale_pngs)) unlink(stale_pngs)

for (start_year in FORECAST_START_YEARS) {
  scatter <- all_scores %>%
    filter(
      forecast_start_year == start_year,
      as.character(model) %in% plot_models,
      is.finite(observed), is.finite(pred_median)
    ) %>%
    mutate(
      model = factor(as.character(model), levels = plot_models),
      sex = factor(sex, levels = c("Female", "Male"))
    )

  if (!nrow(scatter)) {
    warning("No observed-versus-predicted data for holdout origin ", start_year,
            "; no figure was written.", call. = FALSE)
    next
  }

  # Each point is a state-year-sex prediction. The diagonal denotes perfect
  # prediction; colour and shape distinguish the primary model and benchmarks.
  p_scatter <- ggplot(scatter, aes(x = observed, y = pred_median,
                                    colour = model, shape = model)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                linewidth = 0.45, colour = "grey45") +
    geom_point(alpha = 0.78, size = 1.65, stroke = 0.25) +
    facet_wrap(~sex, nrow = 1) +
    scale_colour_manual(values = plot_colours, breaks = plot_models, drop = FALSE) +
    scale_shape_manual(values = plot_shapes, breaks = plot_models, drop = FALSE) +
    coord_equal() +
    labs(
      x = "Observed life expectancy at birth (years)",
      y = "Predicted life expectancy at birth (years)",
      colour = "Model", shape = "Model"
    ) +
    base_theme

  ggsave(
    file.path(output_figures,
              sprintf("holdout_observed_vs_predicted_start_%d.png", start_year)),
    p_scatter, width = 9.5, height = 5.2, dpi = 600, bg = "white"
  )
}

cat("\nHoldout validation completed.\n")
cat("Prediction and metrics tables: ", output_tables, "\n", sep = "")
cat("Only observed-versus-predicted figures saved under: ", output_figures, "\n", sep = "")
cat("  holdout_observed_vs_predicted_start_2000.png\n")
cat("  holdout_observed_vs_predicted_start_2010.png\n")
cat("Start 2000 means training strictly before 2000, scoring 2000-2019.\n")
cat("Start 2010 means training strictly before 2010, scoring 2010-2019.\n")
