# Shared engine for prior-predictive checks and model fitting.
# Run codes/01_run_models.R from the project root.

suppressPackageStartupMessages({
  library(here)
  library(rstan)
  library(tidyverse)
  library(bayesplot)
})

if (!exists("RUN_MODE", envir = .GlobalEnv, inherits = FALSE)) {
  stop("RUN_MODE is not set. Run one of codes/01_run_models.R.", call. = FALSE)
}
RUN_MODE <- as.integer(get("RUN_MODE", envir = .GlobalEnv))
if (!RUN_MODE %in% 1:4) stop("RUN_MODE must be 1, 2, 3, or 4.", call. = FALSE)

rstan_options(auto_write = TRUE)
options(mc.cores = min(4L, parallel::detectCores(logical = TRUE)))
set.seed(7317)

# Reuse saved fits by default. Set TRUE after changing data, priors, or Stan code.
FORCE_REFIT <- FALSE

# Fitting defaults are shared across posterior runs; prior predictive uses fewer chains.
CHAINS <- 4L
ITER <- 2000L
WARMUP <- 1000L
PRIOR_CHAINS <- 2L
PRIOR_ITER <- 1200L
PRIOR_WARMUP <- 600L
ADAPT_DELTA <- 0.95
MAX_TREEDEPTH <- 10L

project_dir <- here::here()
prep_file <- here::here("fits", "prep_objects.rds")
stan_file <- here::here("stan", "life_expectancy_model.stan")
if (!file.exists(prep_file)) {
  stop("Preparation object not found: ", prep_file,
       "\nRun codes/00_prepare_stan_data.R first to build prep_objects.rds from data/model_data.csv and data/groups.csv.", call. = FALSE)
}
if (!file.exists(stan_file)) stop("Unified Stan file not found: ", stan_file, call. = FALSE)

for (d in c("fits", "figures", "tables", "codes", "stan")) {
  dir.create(here::here(d), recursive = TRUE, showWarnings = FALSE)
}
for (d in c(file.path("tables", "prior_predictive"), file.path("figures", "prior_predictive"))) {
  dir.create(here::here(d), recursive = TRUE, showWarnings = FALSE)
}

# Prior names and outputs. Mode 1 is an operation (prior-only checks), not a fourth prior.
prior_labels <- c(
  `2` = "Main model (regularised baseline)",
  `3` = "Tighter-prior sensitivity",
  `4` = "Wide-prior sensitivity"
)
mode_folder <- c(`2` = "main_regularised", `3` = "tighter_prior_sensitivity", `4` = "wide_prior_sensitivity")
fit_path_for_mode <- function(prior_mode) {
  file.path(here::here("fits"), switch(as.character(prior_mode),
    `2` = "fit_main_regularised.rds",
    `3` = "fit_tighter_priors.rds",
    `4` = "fit_wide_priors.rds"))
}

# Keep the user's exploratory Health Index allocation, including West Bengal as Front-runner.
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
normalise_state_key <- function(x) {
  z <- tolower(gsub("[^[:alnum:]]", "", as.character(x)))
  dplyr::case_when(z == "chhatisgarh" ~ "chhattisgarh",
                   z == "orissa" ~ "odisha",
                   TRUE ~ z)
}
state_tier_map <- bind_rows(
  state_tier_map,
  tibble(state_name = c("Chhatisgarh", "Orissa"), tier = c("Achiever", "Aspirant"))
) %>%
  mutate(state_key = normalise_state_key(state_name)) %>%
  distinct(state_key, .keep_all = TRUE)

prep <- readRDS(prep_file)
df <- as_tibble(prep$df)
stan_base <- prep$stan_data_full
state_labels <- as.character(prep$state_levels_full)
N_groups <- as.integer(stan_base$N_groups)
forecast_years <- as.numeric(prep$forecast_years_full)
N_tiers <- length(tier_levels)

if (nrow(df) != stan_base$N || length(state_labels) != N_groups) {
  stop("prep_objects.rds dimension checks failed.", call. = FALSE)
}
if (!all(c("state", "year", "g", "e0_f", "e0_m") %in% names(df))) {
  stop("prep$df must contain state, year, g, e0_f and e0_m.", call. = FALSE)
}
label_keys <- normalise_state_key(state_labels)
matched_tiers <- match(label_keys, state_tier_map$state_key)
if (anyNA(matched_tiers)) {
  stop("No Health Index tier assigned to: ", paste(state_labels[is.na(matched_tiers)], collapse = ", "),
       ". Edit state_tier_map explicitly.", call. = FALSE)
}
tier_for_state <- state_tier_map$tier[matched_tiers]
tier_id_for_state <- as.integer(match(tier_for_state, tier_levels))
if (anyNA(tier_id_for_state)) stop("Unrecognised tier label.", call. = FALSE)

first_idx <- vapply(seq_len(N_groups), function(s) {
  ii <- which(df$g == s)
  if (!length(ii)) stop("No observations for group index ", s)
  ii[which.min(df$year[ii])]
}, integer(1))
last_idx <- vapply(seq_len(N_groups), function(s) {
  ii <- which(df$g == s)
  if (!length(ii)) stop("No observations for group index ", s)
  ii[which.max(df$year[ii])]
}, integer(1))
first_year <- as.numeric(df$year[first_idx])
last_year <- as.numeric(df$year[last_idx])
stopifnot(sum(stan_base$is_first) == N_groups,
          all(stan_base$is_first[first_idx] == 1L),
          all(stan_base$prev_idx[stan_base$is_first == 0L] >= 1L),
          all(stan_base$delta_t[stan_base$is_first == 0L] > 0))

# Explicit prior register is saved with the analyses.
prior_register <- tibble::tribble(
  ~prior_mode, ~parameter, ~prior, ~demographic_or_statistical_rationale,
  2L, "C_f", "Normal(87, 2)", "Main/original regularised baseline: plausible target centre, not a ceiling",
  2L, "C_m", "Normal(82, 2)", "Main/original regularised baseline: plausible target centre, not a ceiling",
  2L, "tier contrasts", "Raw Normal(0, 4)", "Main baseline with moderate cross-tier variation; less shrinkage than mode 3",
  2L, "tau_C", "Half-Normal(3)", "Main baseline between-state equilibrium heterogeneity; less shrinkage than mode 3",
  2L, "L_Omega_C and L_Omega_dev", "LKJ-Cholesky(2)", "Regularises cross-sex dependence of state targets and OU deviations",
  2L, "z_C", "Standard Normal", "Standardised state-specific equilibrium deviations",
  2L, "phi_covid_f/m", "Beta(2, 2)", "Main baseline prior for geometric post-2021 shock decay fractions",
  2L, "log kappa_f", "Normal(-1, 2.5)", "Main regularised baseline retains this broad rate prior; mode 3 tests stronger demographic regularisation",
  2L, "delta_log_kappa", "Normal(0, 0.5)", "Main baseline uncertainty in the sex-specific convergence-rate ratio; narrower in mode 3",
  2L, "log lambda gap", "Normal(log(0.20), 0.7)", "Main baseline reversion-rate prior; mode 3 regularises temporary deviations to revert faster",
  2L, "sigma_dev", "Half-t(df=3, scale=1.5)", "Main baseline half-t process prior; mode 3 uses a tighter half-normal scale prior",
  2L, "sigma_obs", "Half-t(df=3, scale=1.5)", "Main baseline observation SD prior; the SD is estimated rather than fixed",
  2L, "COVID shocks", "Normal(0, 5)", "Main baseline COVID shock prior in years; mode 3 downweights extreme shocks",
  3L, "C_f", "Normal(87, 2)", "Same centre as the main baseline; no hard ceiling",
  3L, "C_m conditional on C_f", "Normal(C_f - 4.5, 1.5)", "Tighter sensitivity prior on the long-run sex difference; uncertainty remains",
  3L, "tier contrasts", "Raw Normal(0, 2)", "Stronger shrinkage of tier differences; no ordering is imposed",
  3L, "tau_C", "Half-Normal(2)", "Stronger shrinkage of state equilibrium heterogeneity while allowing departures",
  3L, "log kappa_f", "Normal(log(log(2)/35), 0.65)", "Centres smooth-trend half-life at 35 years; approximate 95% prior range 10–126 years",
  3L, "delta_log_kappa", "Normal(0, 0.35)", "Allows sex-specific rates; rate ratio has approx. 95% prior range 0.50–1.99",
  3L, "log lambda gap", "Normal(log(0.05), 0.55)", "Temporary deviations revert faster; broadly allows multi-year to multi-decade decay",
  3L, "sigma_dev", "Half-Normal(1.0)", "Stronger shrinkage of process diffusion; stochastic trajectory variation remains estimated",
  3L, "sigma_obs", "Half-Normal(0.75)", "Stronger shrinkage of observation error; the SD remains estimated, not fixed",
  3L, "COVID shocks", "Normal(0, 3)", "Allows multi-year shocks but downweights extremes more than the main baseline",
  4L, "C_f", "Normal(87, 5)", "Wide equilibrium hypermean sensitivity prior; not a ceiling",
  4L, "C_m", "Normal(82, 5)", "Wide equilibrium hypermean sensitivity prior; not a ceiling",
  4L, "tier contrasts", "Raw Normal(0, 8)", "Very wide cross-tier variation sensitivity",
  4L, "tau_C", "Half-Normal(5)", "Wide between-state target heterogeneity sensitivity",
  4L, "L_Omega_C and L_Omega_dev", "LKJ-Cholesky(1)", "Uniform correlation-matrix prior for wide sensitivity",
  4L, "z_C", "Standard Normal", "Standardised state-specific equilibrium deviations",
  4L, "phi_covid_f/m", "Beta(1, 1)", "Uniform post-2021 shock-decay fractions",
  4L, "log kappa_f", "Normal(-1, 3)", "Deliberately diffuse rate prior for stress-testing prior dependence",
  4L, "delta_log_kappa", "Normal(0, 1)", "Wide uncertainty in sex-specific convergence-rate ratio",
  4L, "log lambda gap", "Normal(log(0.20), 1)", "Wide short-term reversion-rate sensitivity",
  4L, "sigma_dev", "Half-t(df=3, scale=3)", "Wide latent-process scale sensitivity",
  4L, "sigma_obs", "Half-t(df=3, scale=3)", "Wide observation-error scale sensitivity",
  4L, "COVID shocks", "Normal(0, 8)", "Wide shock sensitivity in years"
)
prior_register <- prior_register %>%
  mutate(model_role = recode(as.character(prior_mode),
    `2` = "Main/original regularised baseline",
    `3` = "Tighter-prior sensitivity",
    `4` = "Wide-prior sensitivity"), .before = prior_mode)
write_csv(prior_register, here::here("tables", "prior_specification.csv"))
write_csv(tibble(state = state_labels, group_index = seq_len(N_groups), tier = tier_for_state,
                 tier_id = tier_id_for_state),
          here::here("tables", "state_health_index_tiers_used.csv"))

build_stan_data <- function(prior_mode, prior_only = 0L, compute_yrep = 1L) {
  stopifnot(prior_mode %in% 2:4)
  list(
    N = as.integer(stan_base$N),
    g = as.integer(stan_base$g),
    is_first = as.integer(stan_base$is_first),
    prev_idx = as.integer(stan_base$prev_idx),
    delta_t = as.numeric(stan_base$delta_t),
    y_f = as.numeric(stan_base$y_f),
    y_m = as.numeric(stan_base$y_m),
    N_groups = N_groups,
    N_tiers = as.integer(N_tiers),
    tier_id = tier_id_for_state,
    is_2020 = as.integer(stan_base$is_2020),
    is_2021 = as.integer(stan_base$is_2021),
    is_post2021 = as.integer(stan_base$is_post2021),
    years_since_2021 = as.numeric(stan_base$years_since_2021),
    N_forecast = as.integer(length(forecast_years)),
    forecast_years = forecast_years,
    first_year = first_year,
    last_year = last_year,
    last_idx = as.integer(last_idx),
    prior_mode = as.integer(prior_mode),
    prior_only = as.integer(prior_only),
    eag_group = as.integer(stan_base$eag_group),
    compute_yrep = as.integer(compute_yrep)
  )
}

summarise_draws <- function(x) {
  if (!all(is.finite(x))) {
    stop("Found non-finite forecast draws. Do not silently drop these; inspect the generated quantities and refit.", call. = FALSE)
  }
  # Central and tail quantiles make prior sensitivity easy to inspect.
  probs <- c(.01, .025, .05, .25, .50, .75, .95, .975, .99)
  qs <- as.numeric(quantile(x, probs = probs, names = FALSE, type = 8))
  tibble(
    q01 = qs[1], q025 = qs[2], q05 = qs[3], q25 = qs[4], median = qs[5],
    q75 = qs[6], q95 = qs[7], q975 = qs[8], q99 = qs[9],
    width50 = qs[6] - qs[4], width90 = qs[7] - qs[3], width95 = qs[8] - qs[2]
  )
}

extract_forecast_summary <- function(fit, model_label, prior_mode) {
  ex <- rstan::extract(fit, pars = c("latent_forecast_f", "latent_forecast_m",
                                    "e0_forecast_f", "e0_forecast_m"), permuted = TRUE)
  rows <- list(); k <- 0L
  for (forecast_type in c("Latent trajectory", "Future observation")) {
    for (sex in c("Female", "Male")) {
      suffix <- if (sex == "Female") "f" else "m"
      arr_name <- if (forecast_type == "Latent trajectory") paste0("latent_forecast_", suffix) else paste0("e0_forecast_", suffix)
      arr <- ex[[arr_name]]
      if (is.null(arr) || length(dim(arr)) != 3L || dim(arr)[2] != N_groups || dim(arr)[3] != length(forecast_years)) {
        stop("Unexpected generated quantity dimensions for ", arr_name, call. = FALSE)
      }
      for (gi in seq_len(N_groups)) for (ti in seq_along(forecast_years)) {
        k <- k + 1L
        rows[[k]] <- tibble(model_label = model_label, prior_mode = prior_mode,
                            forecast_type = forecast_type, sex = sex,
                            state = state_labels[gi], g = gi, year = forecast_years[ti]) %>%
          bind_cols(summarise_draws(arr[, gi, ti]))
      }
    }
  }
  bind_rows(rows)
}

summary_parameter_diagnostics <- function(fit) {
  smry <- as.data.frame(rstan::summary(fit)$summary) %>% rownames_to_column("parameter")
  omit_pattern <- paste0(
    "^(y_rep_f|y_rep_m|mu_curve_f|mu_curve_m|C_tier_f|C_tier_m|tier_effect_f|tier_effect_m|",
    "elapsed_since_first|baseline_f|baseline_m|covid_shock_f|covid_shock_m|",
    "trend_forecast_f|trend_forecast_m|mean_forecast_f|mean_forecast_m|",
    "latent_forecast_f|latent_forecast_m|e0_forecast_f|e0_forecast_m|",
    "u_pred_mean|u_pred_cov|u_filter_mean|u_filter_cov|marginal_loglik|",
    "I2|L_scaled|L_Omega_dev|L_Omega_C|Omega_dev|Omega_C|L_C_scaled)(\\[|$)"
  )
  d <- smry %>% filter(parameter != "lp__") %>%
    filter(!grepl(omit_pattern, parameter) |
             parameter %in% c("L_Omega_dev[2,1]", "L_Omega_C[2,1]"))
  bfmi <- as.numeric(rstan::get_bfmi(fit))
  tibble(
    divergences = as.integer(rstan::get_num_divergent(fit)),
    max_treedepth_hit = as.integer(rstan::get_num_max_treedepth(fit)),
    min_bfmi = min(bfmi, na.rm = TRUE),
    low_bfmi_chains = sum(bfmi < .3, na.rm = TRUE),
    max_rhat = max(d$Rhat[is.finite(d$Rhat)], na.rm = TRUE),
    min_n_eff = min(d$n_eff[is.finite(d$n_eff)], na.rm = TRUE),
    nonfinite_rhat = sum(!is.finite(d$Rhat)),
    nonfinite_n_eff = sum(!is.finite(d$n_eff))
  )

}

# A single explicit PASS/FAIL summary is written for every sampled fit.
# Thresholds: zero divergences; zero max-treedepth hits; BFMI >= 0.30 in every
# chain; max R-hat <= 1.01; minimum effective sample size >= 100; finite diagnostics.
add_diagnostic_status <- function(diagnostics) {
  d <- diagnostics %>% mutate(
    check_zero_divergences = divergences == 0L,
    check_no_treedepth_hits = max_treedepth_hit == 0L,
    check_bfmi = is.finite(min_bfmi) && min_bfmi >= 0.30 && low_bfmi_chains == 0L,
    check_rhat = is.finite(max_rhat) && max_rhat <= 1.01,
    check_ess = is.finite(min_n_eff) && min_n_eff >= 100,
    check_finite_rhat = nonfinite_rhat == 0L,
    check_finite_ess = nonfinite_n_eff == 0L
  )
  check_cols <- c("check_zero_divergences", "check_no_treedepth_hits", "check_bfmi",
                  "check_rhat", "check_ess", "check_finite_rhat", "check_finite_ess")
  passed <- all(unlist(d[1, check_cols], use.names = FALSE))
  d$sampler_status <- if (passed) "PASS" else "FAIL"
  failed <- check_cols[!as.logical(d[1, check_cols])]
  d$failed_checks <- if (!length(failed)) "None" else paste(failed, collapse = "; ")
  d
}

print_check_line <- function(label, value_text, passed) {
  cat(sprintf("  %-35s %s  [%s]\n", label, value_text,
              if (isTRUE(passed)) "PASS" else "FAIL"))
}

save_sampler_diagnostics <- function(fit, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  diagnostics <- add_diagnostic_status(summary_parameter_diagnostics(fit))
  write_csv(diagnostics, file.path(output_dir, "sampler_diagnostics.csv"))
  write_csv(diagnostics %>% select(sampler_status, failed_checks, divergences,
    max_treedepth_hit, min_bfmi, low_bfmi_chains, max_rhat, min_n_eff,
    nonfinite_rhat, nonfinite_n_eff,
    check_zero_divergences, check_no_treedepth_hits, check_bfmi, check_rhat,
    check_ess, check_finite_rhat, check_finite_ess),
    file.path(output_dir, "diagnostic_status.csv"))

  bfmi_values <- as.numeric(rstan::get_bfmi(fit))
  bfmi_table <- tibble(chain = seq_along(bfmi_values), BFMI = bfmi_values,
                       meets_threshold = is.finite(bfmi_values) & bfmi_values >= 0.30)
  write_csv(bfmi_table, file.path(output_dir, "bfmi_by_chain.csv"))
  full_summary <- as.data.frame(rstan::summary(fit)$summary) %>% rownames_to_column("parameter")
  write_csv(full_summary %>% filter(parameter == "lp__") %>% select(parameter, mean, sd, n_eff, Rhat),
            file.path(output_dir, "lp_diagnostic_separate.csv"))

  d <- diagnostics[1, ]
  cat("\nSAMPLER DIAGNOSTICS\n")
  print_check_line("Divergences:", sprintf("%d (required = 0)", d$divergences), d$check_zero_divergences)
  print_check_line("Maximum tree-depth hits:", sprintf("%d (required = 0)", d$max_treedepth_hit), d$check_no_treedepth_hits)
  print_check_line("Minimum BFMI:", sprintf("%.3f (required >= 0.300; chains below = %d)", d$min_bfmi, d$low_bfmi_chains), d$check_bfmi)
  print_check_line("Maximum R-hat:", sprintf("%.4f (required <= 1.0100)", d$max_rhat), d$check_rhat)
  print_check_line("Minimum effective sample size:", sprintf("%.1f (required >= 100)", d$min_n_eff), d$check_ess)
  print_check_line("Non-finite R-hat values:", sprintf("%d (required = 0)", d$nonfinite_rhat), d$check_finite_rhat)
  print_check_line("Non-finite ESS values:", sprintf("%d (required = 0)", d$nonfinite_n_eff), d$check_finite_ess)
  for (i in seq_len(nrow(bfmi_table))) {
    print_check_line(sprintf("BFMI chain %d:", bfmi_table$chain[i]),
                     sprintf("%.3f (required >= 0.300)", bfmi_table$BFMI[i]),
                     bfmi_table$meets_threshold[i])
  }
  cat("OVERALL SAMPLER STATUS: ", d$sampler_status, "\n", sep = "")
  if (d$sampler_status == "FAIL") cat("Failed checks: ", d$failed_checks, "\n", sep = "")
  diagnostics
}

validate_generated_quantities <- function(fit, output_dir) {
  # Check important forecast and predictive arrays before any summaries are written.
  variables <- c("latent_forecast_f", "latent_forecast_m", "e0_forecast_f", "e0_forecast_m",
                 "y_rep_f", "y_rep_m")
  rows <- lapply(variables, function(variable) {
    one <- rstan::extract(fit, pars = variable, permuted = TRUE)[[variable]]
    total <- length(one)
    n_nonfinite <- sum(!is.finite(one))
    data.frame(variable = variable, n_values = total, n_nonfinite = n_nonfinite,
               status = if (n_nonfinite == 0L) "PASS" else "FAIL")
  })
  status <- dplyr::bind_rows(rows)
  readr::write_csv(status, file.path(output_dir, "generated_quantities_status.csv"))
  cat("\nGENERATED-QUANTITIES CHECKS\n")
  for (i in seq_len(nrow(status))) {
    cat(sprintf("  %-25s non-finite = %d / %d  [%s]\n",
                status$variable[i], status$n_nonfinite[i], status$n_values[i], status$status[i]))
  }
  overall <- if (all(status$status == "PASS")) "PASS" else "FAIL"
  cat("OVERALL GENERATED-QUANTITIES STATUS: ", overall, "\n", sep = "")
  if (overall != "PASS") {
    failed <- paste(status$variable[status$status != "PASS"], collapse = ", ")
    stop("Undefined/non-finite generated quantities: ", failed,
         ". See generated_quantities_status.csv; forecasts were not summarised.", call. = FALSE)
  }
  invisible(status)
}

save_energy_plot <- function(fit, path) {
  p <- bayesplot::mcmc_nuts_energy(bayesplot::nuts_params(fit), merge_chains = FALSE) +
    ggplot2::labs(title = NULL, subtitle = NULL)
  ggsave(path, p, width = 9, height = 6, dpi = 600, bg = "white")
}

save_forecast_figures <- function(fit, model_label, output_dir, output_table_dir,
                                 max_plot_draws = 1500L) {
  # The figure combines (1) historical fitted latent values conditional on the
  # observed series and (2) latent forecasts started immediately after each
  # state's own final observed year. Both sexes use solid historical lines and
  # dashed forecast lines; ribbons show pointwise 95% credible intervals.
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(output_table_dir, recursive = TRUE, showWarnings = FALSE)

  needed <- c(
    "mu_curve_f", "mu_curve_m", "u_filter_mean", "u_filter_cov",
    "baseline_f", "baseline_m", "C_state_f", "C_state_m",
    "kappa_f", "kappa_m", "lambda_dev", "sigma_dev", "L_Omega_dev",
    "delta2020_f", "delta2021_f", "delta2020_m", "delta2021_m",
    "phi_covid_f", "phi_covid_m"
  )
  # Extract only the quantities needed below; extracting every saved variable
  # first would duplicate large posterior arrays in memory.
  ex <- rstan::extract(fit, pars = needed, permuted = TRUE, inc_warmup = FALSE)
  n_total <- dim(ex$mu_curve_f)[1]
  n_draws <- min(as.integer(max_plot_draws), n_total)
  if (n_draws < 100L) stop("Too few posterior draws available for the trajectory figure.", call. = FALSE)
  draw_ids <- unique(as.integer(round(seq(1, n_total, length.out = n_draws))))
  n_draws <- length(draw_ids)
  n_obs <- nrow(df)

  # Make plot simulations reproducible without changing the caller's RNG state.
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(33107)

  if (length(dim(ex$u_filter_mean)) != 3L ||
      length(dim(ex$u_filter_cov)) != 4L ||
      dim(ex$u_filter_mean)[2] != n_obs || dim(ex$u_filter_mean)[3] != 2L ||
      dim(ex$u_filter_cov)[2] != n_obs || dim(ex$u_filter_cov)[3] != 2L ||
      dim(ex$u_filter_cov)[4] != 2L) {
    stop("Unexpected dimensions for u_filter_mean or u_filter_cov.", call. = FALSE)
  }

  # Posterior draws of historical fitted latent values. The filtering state is
  # conditioned on observations up to and including the current observation;
  # observation-level residual noise is deliberately excluded.
  fitted_draws_f <- matrix(NA_real_, nrow = n_draws, ncol = n_obs)
  fitted_draws_m <- matrix(NA_real_, nrow = n_draws, ncol = n_obs)
  for (i in seq_len(n_obs)) {
    mu_f <- as.numeric(ex$mu_curve_f[draw_ids, i])
    mu_m <- as.numeric(ex$mu_curve_m[draw_ids, i])
    umean_f <- as.numeric(ex$u_filter_mean[draw_ids, i, 1])
    umean_m <- as.numeric(ex$u_filter_mean[draw_ids, i, 2])
    var_f <- pmax(as.numeric(ex$u_filter_cov[draw_ids, i, 1, 1]), 0)
    var_m <- pmax(as.numeric(ex$u_filter_cov[draw_ids, i, 2, 2]), 0)
    fitted_draws_f[, i] <- mu_f + umean_f + sqrt(var_f) * rnorm(n_draws)
    fitted_draws_m[, i] <- mu_m + umean_m + sqrt(var_m) * rnorm(n_draws)
  }

  # A paired draw from the bivariate filtered state at a state's last observed
  # year starts the projected paths, ensuring the forecast is anchored to its
  # own last observation rather than being forced to start in 2025 for every state.
  draw_filtered_state <- function(draw_number, obs_index) {
    mu <- as.numeric(ex$u_filter_mean[draw_number, obs_index, ])
    S <- matrix(as.numeric(ex$u_filter_cov[draw_number, obs_index, , ]), nrow = 2L, ncol = 2L)
    S <- 0.5 * (S + t(S))
    if (max(abs(S)) < 1e-14) return(mu)
    eig_min <- min(eigen(S, symmetric = TRUE, only.values = TRUE)$values)
    if (!is.finite(eig_min)) stop("Non-finite filtered-state covariance in trajectory plot.", call. = FALSE)
    if (eig_min < 1e-10) S <- S + diag(1e-10 - eig_min, 2L)
    as.numeric(mu + t(chol(S)) %*% rnorm(2L))
  }

  covid_addition <- function(year, sex) {
    if (sex == "Female") {
      if (year == 2020) return(as.numeric(ex$delta2020_f[draw_ids]))
      if (year == 2021) return(as.numeric(ex$delta2021_f[draw_ids]))
      if (year > 2021) return(as.numeric(ex$delta2021_f[draw_ids] * ex$phi_covid_f[draw_ids]^(year - 2021)))
    } else {
      if (year == 2020) return(as.numeric(ex$delta2020_m[draw_ids]))
      if (year == 2021) return(as.numeric(ex$delta2021_m[draw_ids]))
      if (year > 2021) return(as.numeric(ex$delta2021_m[draw_ids] * ex$phi_covid_m[draw_ids]^(year - 2021)))
    }
    rep(0, n_draws)
  }

  quantile_row <- function(x, state, year, sex, segment) {
    q <- as.numeric(quantile(x, probs = c(.025, .50, .975), names = FALSE, type = 8))
    tibble(model_label = model_label, state = as.character(state), year = as.numeric(year),
           sex = sex, segment = segment, q025 = q[1], median = q[2], q975 = q[3])
  }

  # Retain each state's last latent-state draw so historical and forecast lines
  # meet at the state-specific forecast origin for both sexes.
  last_u_f <- matrix(NA_real_, nrow = n_draws, ncol = N_groups)
  last_u_m <- matrix(NA_real_, nrow = n_draws, ncol = N_groups)
  for (gi in seq_len(N_groups)) {
    last_i <- last_idx[gi]
    for (j in seq_len(n_draws)) {
      di <- draw_ids[j]
      u_last <- draw_filtered_state(di, last_i)
      last_u_f[j, gi] <- u_last[1]
      last_u_m[j, gi] <- u_last[2]
      fitted_draws_f[j, last_i] <- as.numeric(ex$mu_curve_f[di, last_i]) + u_last[1]
      fitted_draws_m[j, last_i] <- as.numeric(ex$mu_curve_m[di, last_i]) + u_last[2]
    }
  }

  history_summary_rows <- vector("list", n_obs * 2L)
  history_k <- 0L
  for (i in seq_len(n_obs)) {
    history_k <- history_k + 1L
    history_summary_rows[[history_k]] <- quantile_row(
      fitted_draws_f[, i], df$state[i], df$year[i], "Female", "Historical fit")
    history_k <- history_k + 1L
    history_summary_rows[[history_k]] <- quantile_row(
      fitted_draws_m[, i], df$state[i], df$year[i], "Male", "Historical fit")
  }
  history_summary <- bind_rows(history_summary_rows)

  forecast_rows <- list()
  forecast_k <- 0L
  forecast_end <- max(forecast_years, na.rm = TRUE)

  for (gi in seq_len(N_groups)) {
    last_i <- last_idx[gi]
    this_last_year <- as.numeric(df$year[last_i])
    state <- as.character(state_labels[gi])

    # Include the forecast origin as an anchor, then simulate one-year steps
    # from each state's final observed year through the projection endpoint.
    forecast_k <- forecast_k + 1L
    forecast_rows[[forecast_k]] <- quantile_row(fitted_draws_f[, last_i], state, this_last_year,
                                                "Female", "Forecast")
    forecast_k <- forecast_k + 1L
    forecast_rows[[forecast_k]] <- quantile_row(fitted_draws_m[, last_i], state, this_last_year,
                                                "Male", "Forecast")

    if (this_last_year >= forecast_end) next

    prev_u_f <- last_u_f[, gi]
    prev_u_m <- last_u_m[, gi]
    prev_year <- this_last_year
    first_year_g <- min(df$year[df$g == gi], na.rm = TRUE)

    for (yr in seq.int(as.integer(this_last_year) + 1L, as.integer(forecast_end))) {
      dtt <- as.numeric(yr - prev_year)
      lambda <- as.numeric(ex$lambda_dev[draw_ids])
      rho_dev <- exp(-lambda * dtt)
      mult_dev <- ifelse(lambda * dtt < 1e-8, sqrt(dtt),
                         sqrt(pmax(0, -expm1(-2 * lambda * dtt) / (2 * lambda))))

      z1 <- rnorm(n_draws)
      z2 <- rnorm(n_draws)
      l11 <- as.numeric(ex$L_Omega_dev[draw_ids, 1, 1])
      l21 <- as.numeric(ex$L_Omega_dev[draw_ids, 2, 1])
      l22 <- as.numeric(ex$L_Omega_dev[draw_ids, 2, 2])
      zcorr_f <- l11 * z1
      zcorr_m <- l21 * z1 + l22 * z2
      sig_f <- as.numeric(ex$sigma_dev[draw_ids, 1])
      sig_m <- as.numeric(ex$sigma_dev[draw_ids, 2])

      next_u_f <- rho_dev * prev_u_f + sig_f * mult_dev * zcorr_f
      next_u_m <- rho_dev * prev_u_m + sig_m * mult_dev * zcorr_m

      trend_f <- as.numeric(ex$C_state_f[draw_ids, gi]) +
        (as.numeric(ex$baseline_f[draw_ids, gi]) - as.numeric(ex$C_state_f[draw_ids, gi])) *
        exp(-as.numeric(ex$kappa_f[draw_ids]) * (yr - first_year_g)) + covid_addition(yr, "Female")
      trend_m <- as.numeric(ex$C_state_m[draw_ids, gi]) +
        (as.numeric(ex$baseline_m[draw_ids, gi]) - as.numeric(ex$C_state_m[draw_ids, gi])) *
        exp(-as.numeric(ex$kappa_m[draw_ids]) * (yr - first_year_g)) + covid_addition(yr, "Male")

      latent_f <- trend_f + next_u_f
      latent_m <- trend_m + next_u_m
      forecast_k <- forecast_k + 1L
      forecast_rows[[forecast_k]] <- quantile_row(latent_f, state, yr, "Female", "Forecast")
      forecast_k <- forecast_k + 1L
      forecast_rows[[forecast_k]] <- quantile_row(latent_m, state, yr, "Male", "Forecast")

      prev_u_f <- next_u_f
      prev_u_m <- next_u_m
      prev_year <- yr
    }
  }

  history_summary <- history_summary %>% arrange(state, sex, year)
  forecast_summary_plot <- bind_rows(forecast_rows) %>% arrange(state, sex, year)
  write_csv(history_summary,
            file.path(output_table_dir, "historical_fitted_latent_trajectory_summary.csv"))
  write_csv(forecast_summary_plot,
            file.path(output_table_dir, "state_specific_latent_forecast_summary.csv"))

  history <- df %>%
    transmute(state = as.character(state), year = as.numeric(year),
              Female = as.numeric(e0_f), Male = as.numeric(e0_m)) %>%
    pivot_longer(cols = c(Female, Male), names_to = "sex", values_to = "observed")

  p <- ggplot() +
    geom_ribbon(data = history_summary,
                aes(x = year, ymin = q025, ymax = q975, fill = sex,
                    group = interaction(state, sex)),
                alpha = .10, colour = NA) +
    geom_ribbon(data = forecast_summary_plot,
                aes(x = year, ymin = q025, ymax = q975, fill = sex,
                    group = interaction(state, sex)),
                alpha = .15, colour = NA) +
    geom_line(data = history_summary,
              aes(x = year, y = median, colour = sex, linetype = "Historical fit",
                  group = interaction(state, sex)), linewidth = .55) +
    geom_line(data = forecast_summary_plot,
              aes(x = year, y = median, colour = sex, linetype = "Forecast",
                  group = interaction(state, sex)), linewidth = .65) +
    geom_point(data = history,
               aes(x = year, y = observed, shape = sex),
               colour = "grey35", size = .65, alpha = .70) +
    facet_wrap(~state, ncol = 5) +
    scale_colour_manual(values = c("Female" = "#B2182B", "Male" = "#2166AC"), name = "Sex") +
    scale_fill_manual(values = c("Female" = "#B2182B", "Male" = "#2166AC"), name = "Sex") +
    scale_linetype_manual(values = c("Historical fit" = "solid", "Forecast" = "dashed"), name = "Period") +
    scale_shape_manual(values = c("Female" = 16, "Male" = 17), name = "Observed sex") +
    scale_x_continuous(limits = c(min(df$year, na.rm = TRUE), forecast_end),
                       breaks = scales::pretty_breaks(n = 5), expand = expansion(mult = c(.01, .01))) +
    labs(title = NULL, subtitle = NULL,
         x = "Year", y = "Life expectancy at birth (years)") +
    theme_minimal(base_size = 9) +
    theme(legend.position = "bottom", strip.background = element_blank(),
          strip.text = element_text(size = 8), plot.title = element_blank(),
          plot.subtitle = element_blank())
  ggsave(file.path(output_dir, "latent_trajectories_both_sexes.png"), p,
         width = 16, height = 13, dpi = 600, bg = "white")

  cat("Saved historical fitted latent trajectories and state-specific forecasts for ",
      model_label, " to ", output_dir, ".\n", sep = "")
}

save_ppc <- function(fit, model_label, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  include <- stan_base$is_first == 0L
  idx <- which(include)
  if (!length(idx)) stop("No non-initial observations available for PPC.", call. = FALSE)

  index_key <- df %>% mutate(.original_index = row_number()) %>% filter(include) %>%
    transmute(ppc_index = row_number(), state = as.character(state), year = as.numeric(year), g)
  write_csv(index_key, file.path(output_dir, "ppc_observation_index_key.csv"))

  ex <- rstan::extract(fit, pars = c("y_rep_f", "y_rep_m"), permuted = TRUE)
  plot_rows <- list()
  coverage_rows <- list()
  k <- 0L
  for (sx in c("Female", "Male")) {
    obs_all <- if (sx == "Female") stan_base$y_f else stan_base$y_m
    obs <- as.numeric(obs_all[idx])
    rep_all <- if (sx == "Female") ex$y_rep_f[, idx, drop = FALSE] else ex$y_rep_m[, idx, drop = FALSE]
    if (!all(is.finite(rep_all))) {
      stop("Undefined y_rep values for ", model_label, " / ", sx,
           ". PPC must not silently drop missing draws.", call. = FALSE)
    }
    q25 <- apply(rep_all, 2, quantile, probs = .25, names = FALSE, na.rm = FALSE)
    q75 <- apply(rep_all, 2, quantile, probs = .75, names = FALSE, na.rm = FALSE)
    q05 <- apply(rep_all, 2, quantile, probs = .05, names = FALSE, na.rm = FALSE)
    q95 <- apply(rep_all, 2, quantile, probs = .95, names = FALSE, na.rm = FALSE)

    k <- k + 1L
    plot_rows[[k]] <- tibble(
      state = as.character(df$state[idx]), year = as.numeric(df$year[idx]), sex = sx,
      observed = obs, q05 = q05, q25 = q25, q75 = q75, q95 = q95
    )
    coverage_rows[[k]] <- tibble(
      model_label = model_label, sex = sx,
      coverage50 = mean(obs >= q25 & obs <= q75),
      coverage90 = mean(obs >= q05 & obs <= q95)
    )
  }
  ppc_data <- bind_rows(plot_rows)
  observed_data <- ppc_data %>% distinct(state, year, sex, observed)
  p <- ggplot(ppc_data, aes(x = year, colour = sex, group = sex)) +
    geom_linerange(aes(ymin = q05, ymax = q95), linewidth = .35, alpha = .65) +
    geom_linerange(aes(ymin = q25, ymax = q75), linewidth = .8) +
    geom_point(data = observed_data, aes(x = year, y = observed, shape = sex),
               inherit.aes = FALSE, colour = "black", size = .8) +
    facet_wrap(~state, ncol = 5) +
    scale_colour_manual(values = c("Female" = "#B2182B", "Male" = "#2166AC")) +
    scale_shape_manual(values = c("Female" = 16, "Male" = 17)) +
    labs(title = NULL, subtitle = NULL, x = "Year",
         y = "Life expectancy at birth (years)", colour = "Sex", shape = "Observed sex") +
    theme_minimal(base_size = 9) +
    theme(legend.position = "bottom", strip.background = element_blank(),
          strip.text = element_text(size = 8), plot.title = element_blank(),
          plot.subtitle = element_blank())
  ggsave(file.path(output_dir, "ppc_intervals_by_state.png"), p,
         width = 16, height = 13, dpi = 600, bg = "white")

  bind_rows(coverage_rows)
}

# Fit the unified Stan program once. All three posterior runs and prior checks use this exact file.
cat("Project root: ", project_dir, "\nStan file: ", stan_file, "\n", sep = "")
mod <- rstan::stan_model(file = stan_file)

if (RUN_MODE == 1L) {
  # Prior-predictive check all three distributions with no likelihood contribution.
  # The model conditions each state's first level on its initial observed value;
  # the check is therefore conditional on the starting level, not unconditional at birth.
  all_prior_forecasts <- list(); all_prior_outside <- list(); all_prior_ppc <- list(); all_prior_quantiles <- list()
  for (pm in 2:4) {
    model_label <- unname(prior_labels[as.character(pm)])
    fit_file <- file.path(here::here("fits"), paste0("fit_prior_predictive_", c("main_regularised", "tighter_priors", "wide_priors")[pm - 1L], ".rds"))
    cat("\nPrior-predictive check for ", model_label, " (prior_mode=", pm, ")\n", sep = "")
    if (file.exists(fit_file) && !isTRUE(FORCE_REFIT)) {
      cat("LOADING EXISTING PRIOR-PREDICTIVE FIT: ", fit_file, "\n", sep = "")
      fit <- readRDS(fit_file)
    } else {
      cat("FITTING FRESH PRIOR-PREDICTIVE RUN: ", fit_file, "\n", sep = "")
      fit <- rstan::sampling(
        mod, data = build_stan_data(prior_mode = pm, prior_only = 1L, compute_yrep = 1L),
        iter = PRIOR_ITER, warmup = PRIOR_WARMUP, chains = PRIOR_CHAINS,
        control = list(adapt_delta = ADAPT_DELTA, max_treedepth = MAX_TREEDEPTH),
        seed = 9100 + pm, refresh = 100
      )
      saveRDS(fit, fit_file)
      cat("SAVED NEW PRIOR-PREDICTIVE FIT: ", fit_file, "\n", sep = "")
    }
    prior_diag_dir <- here::here("tables", "prior_predictive", paste0("mode_", pm))
    prior_diag <- save_sampler_diagnostics(fit, prior_diag_dir)
    if (!identical(prior_diag$sampler_status[1], "PASS")) {
      stop("Prior-only sampling diagnostics FAILED for ", model_label,
           ". See tables/prior_predictive/mode_", pm, "/diagnostic_status.csv.", call. = FALSE)
    }
    validate_generated_quantities(fit, prior_diag_dir)
    fs <- extract_forecast_summary(fit, model_label, pm)
    all_prior_forecasts[[as.character(pm)]] <- fs
    ex <- rstan::extract(fit, pars = c("e0_forecast_f", "e0_forecast_m", "y_rep_f", "y_rep_m"), permuted = TRUE)
    for (sx in c("Female", "Male")) {
      suffix <- if (sx == "Female") "f" else "m"
      arr <- ex[[paste0("e0_forecast_", suffix)]]
      for (ti in seq_along(forecast_years)) {
        vals <- as.vector(arr[, , ti])
        if (!all(is.finite(vals))) stop("Non-finite prior-predictive forecast draws.", call. = FALSE)
        quantile_probs <- c(.01, .025, .05, .25, .50, .75, .95, .975, .99)
        q <- as.numeric(quantile(vals, probs = quantile_probs, names = FALSE, type = 8))
        all_prior_quantiles[[length(all_prior_quantiles) + 1L]] <- tibble(
          prior_mode = pm, model_label = model_label, sex = sx, year = forecast_years[ti],
          q01 = q[1], q025 = q[2], q05 = q[3], q25 = q[4], median = q[5],
          q75 = q[6], q95 = q[7], q975 = q[8], q99 = q[9],
          n_state_draws = length(vals),
          distribution_scope = "Pooled across states and prior-predictive draws"
        )
        all_prior_outside[[length(all_prior_outside) + 1L]] <- tibble(
          prior_mode = pm, model_label = model_label, sex = sx, year = forecast_years[ti],
          share_outside_20_110 = mean(vals < 20 | vals > 110),
          share_outside_50_100 = mean(vals < 50 | vals > 100)
        )
      }
    }
    # Prior-predictive historical check; exclude first observations used as anchors.
    ppc_dir <- here::here("figures", "prior_predictive")
    dir.create(ppc_dir, recursive = TRUE, showWarnings = FALSE)
    include <- stan_base$is_first == 0L
    idx <- which(include)
    for (sx in c("Female", "Male")) {
      obs <- if (sx == "Female") stan_base$y_f[idx] else stan_base$y_m[idx]
      rep_all <- if (sx == "Female") ex$y_rep_f[, idx, drop = FALSE] else ex$y_rep_m[, idx, drop = FALSE]
      if (!all(is.finite(rep_all))) stop("Undefined prior predictive y_rep values for ", model_label, " / ", sx, call. = FALSE)
      z <- if (sx == "Female") 1L else 2L
      qlo <- apply(rep_all, 2, quantile, probs = .05, na.rm = FALSE)
      qhi <- apply(rep_all, 2, quantile, probs = .95, na.rm = FALSE)
      qlo50 <- apply(rep_all, 2, quantile, probs = .25, na.rm = FALSE)
      qhi50 <- apply(rep_all, 2, quantile, probs = .75, na.rm = FALSE)
      all_prior_ppc[[length(all_prior_ppc) + 1L]] <- tibble(
        prior_mode = pm, model_label = model_label, sex = sx,
        coverage50 = mean(obs >= qlo50 & obs <= qhi50),
        coverage90 = mean(obs >= qlo & obs <= qhi)
      )
      set.seed(9200 + pm + z)
      draw_id <- unique(round(seq(1, nrow(rep_all), length.out = min(300, nrow(rep_all)))))
      p <- bayesplot::ppc_intervals(y = obs, yrep = rep_all[draw_id, , drop = FALSE],
                                    prob = .50, prob_outer = .90) +
        labs(title = NULL, subtitle = NULL, x = "State-year observation index", y = "Life expectancy at birth (years)") +
        theme_minimal(base_size = 10)
      ggsave(file.path(ppc_dir, paste0("prior_ppc_mode", pm, "_", tolower(sx), ".png")), p,
             width = 13, height = 7, dpi = 600, bg = "white")
    }
  }
  prior_forecasts <- bind_rows(all_prior_forecasts)
  outside_summary <- bind_rows(all_prior_outside)
  prior_ppc_summary <- bind_rows(all_prior_ppc)
  prior_quantile_ranges <- bind_rows(all_prior_quantiles)
  write_csv(prior_forecasts, here::here("tables", "prior_predictive", "prior_predictive_forecast_quantiles_all_states_years.csv"))
  write_csv(prior_quantile_ranges, here::here("tables", "prior_predictive", "prior_predictive_quantile_ranges_001_099.csv"))
  write_csv(outside_summary, here::here("tables", "prior_predictive", "prior_predictive_plausibility_shares.csv"))
  write_csv(prior_ppc_summary, here::here("tables", "prior_predictive", "prior_predictive_historical_interval_coverage.csv"))

  horizons <- intersect(c(2030, 2040, 2050), forecast_years)
  pp <- prior_forecasts %>% filter(forecast_type == "Latent trajectory", year %in% horizons) %>%
    mutate(model_label = factor(model_label, levels = unname(prior_labels)),
           state = factor(state, levels = rev(state_labels)), year = factor(year, levels = horizons))
  p <- ggplot(pp, aes(x = state, y = median, colour = model_label)) +
    geom_linerange(aes(ymin = q025, ymax = q975), position = position_dodge(width = .72), linewidth = .35, alpha = .65) +
    geom_linerange(aes(ymin = q05, ymax = q95), position = position_dodge(width = .72), linewidth = .8) +
    geom_point(position = position_dodge(width = .72), size = 1.15) +
    facet_grid(sex ~ year) + coord_flip() +
    labs(title = NULL, subtitle = NULL, x = NULL, y = "Life expectancy at birth (years)", colour = "Prior setting") +
    theme_minimal(base_size = 10) +
    theme(axis.text.y = element_text(size = 6.5), legend.position = "bottom",
          panel.grid.major.y = element_blank(), strip.text = element_text(face = "bold"))
  ggsave(here::here("figures", "prior_predictive", "prior_predictive_forecast_intervals_2030_2040_2050.png"),
         p, width = 15, height = 11, dpi = 600, bg = "white")
  p2 <- outside_summary %>% pivot_longer(c(share_outside_20_110, share_outside_50_100),
                                         names_to = "band", values_to = "share") %>%
    mutate(band = recode(band, share_outside_20_110 = "Outside 20–110 years",
                         share_outside_50_100 = "Outside 50–100 years"),
           model_label = factor(model_label, levels = unname(prior_labels))) %>%
    ggplot(aes(x = year, y = share, colour = model_label)) +
    geom_line(linewidth = .65) + geom_point(size = 1.1) +
    facet_grid(sex ~ band) + scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(title = NULL, subtitle = NULL, x = "Forecast year",
         y = "Share of prior-predictive draws", colour = "Prior setting") +
    theme_minimal(base_size = 10) + theme(legend.position = "bottom")
  ggsave(here::here("figures", "prior_predictive", "prior_predictive_out_of_range_shares.png"),
         p2, width = 12, height = 7, dpi = 600, bg = "white")
  cat("\nMODE 1 COMPLETE: prior-predictive outputs are in figures/prior_predictive and tables/prior_predictive.\n")
  cat("Quantile ranges include 0.01, 0.025, 0.05, 0.25, 0.50, 0.75, 0.95, 0.975 and 0.99.\n")
  cat("Review prior-predictive plausibility before interpreting the main, tighter, or wide specification.\n")
} else {
  pm <- RUN_MODE
  model_label <- unname(prior_labels[as.character(pm)])
  out_stub <- unname(mode_folder[as.character(pm)])
  out_tables <- here::here("tables", out_stub)
  out_figures <- here::here("figures", out_stub)
  dir.create(out_tables, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_figures, recursive = TRUE, showWarnings = FALSE)

  fit_file <- fit_path_for_mode(pm)
  if (file.exists(fit_file) && !isTRUE(FORCE_REFIT)) {
    cat("\nLOADING EXISTING FIT for ", model_label, "\nFile: ", fit_file, "\n", sep = "")
    fit <- readRDS(fit_file)
  } else {
    cat("\nFITTING MODEL FRESHLY: ", model_label, " (prior_mode=", pm, ")\n", sep = "")
    fit <- rstan::sampling(
      mod, data = build_stan_data(prior_mode = pm, prior_only = 0L, compute_yrep = 1L),
      iter = ITER, warmup = WARMUP, chains = CHAINS,
      control = list(adapt_delta = ADAPT_DELTA, max_treedepth = MAX_TREEDEPTH),
      seed = 7300 + pm, refresh = 200
    )
    saveRDS(fit, fit_file)
    cat("SAVED NEW FIT: ", fit_file, "\n", sep = "")
  }

  full_summary <- as.data.frame(rstan::summary(fit)$summary) %>% rownames_to_column("parameter")
  write_csv(full_summary, file.path(out_tables, "mcmc_summary_all_quantities.csv"))
  diagnostics <- save_sampler_diagnostics(fit, out_tables)
  save_energy_plot(fit, file.path(out_figures, "nuts_energy.png"))
  pass <- identical(diagnostics$sampler_status[1], "PASS")
  if (!pass) stop("Sampler diagnostic screen FAILED for ", model_label,
                  ". Fit and diagnostics were saved. Forecast summaries were not produced; inspect diagnostic_status.csv.", call. = FALSE)

  validate_generated_quantities(fit, out_tables)
  forecast_summary <- extract_forecast_summary(fit, model_label, pm)
  write_csv(forecast_summary, file.path(out_tables, "forecast_summary.csv"))
  write_csv(forecast_summary %>% filter(year == max(forecast_years)),
            file.path(out_tables, "forecast_summary_last_year.csv"))
  ppc <- save_ppc(fit, model_label, out_figures)
  write_csv(ppc, file.path(out_tables, "historical_one_step_ppc_coverage.csv"))
  save_forecast_figures(fit, model_label, out_figures, out_tables)

  # Compact equilibrium and timescale summaries for the report.
  pars <- rstan::extract(fit, pars = c("C_f", "C_m", "delta_C", "C_tier_f", "C_tier_m",
    "C_state_f", "C_state_m", "tau_C", "kappa_f", "kappa_m", "lambda_dev",
    "deviation_half_life", "convergence_half_life_f", "convergence_half_life_m",
    "years_to_90pct_f", "years_to_90pct_m", "sigma_dev", "sigma_obs",
    "long_run_sd_dev_f", "long_run_sd_dev_m", "delta2020_f", "delta2021_f",
    "delta2020_m", "delta2021_m", "phi_covid_f", "phi_covid_m"), permuted = TRUE)
  qsummary <- function(x) tibble(median = median(x), q025 = quantile(x, .025), q975 = quantile(x, .975))
  eq <- bind_rows(
    qsummary(pars$C_f) %>% mutate(quantity = "Female overall equilibrium hypermean", .before = 1),
    qsummary(pars$C_m) %>% mutate(quantity = "Male overall equilibrium hypermean", .before = 1),
    qsummary(pars$delta_C) %>% mutate(quantity = "Male minus female hypermean", .before = 1),
    qsummary(pars$kappa_f) %>% mutate(quantity = "Female convergence rate", .before = 1),
    qsummary(pars$kappa_m) %>% mutate(quantity = "Male convergence rate", .before = 1),
    qsummary(pars$deviation_half_life) %>% mutate(quantity = "Temporary-deviation half-life", .before = 1),
    qsummary(pars$sigma_obs) %>% mutate(quantity = "Observation SD", .before = 1)
  )
  write_csv(eq, file.path(out_tables, "key_parameter_summaries.csv"))
  # State-specific equilibrium targets are important to separate from 2050 forecasts.
  state_eq <- map_dfr(seq_len(N_groups), function(gi) {
    tibble(state = state_labels[gi], group_index = gi, tier = tier_for_state[gi],
           sex = c("Female", "Male"),
           median = c(median(pars$C_state_f[, gi]), median(pars$C_state_m[, gi])),
           q025 = c(quantile(pars$C_state_f[, gi], .025), quantile(pars$C_state_m[, gi], .025)),
           q975 = c(quantile(pars$C_state_f[, gi], .975), quantile(pars$C_state_m[, gi], .975)))
  })
  write_csv(state_eq, file.path(out_tables, "state_equilibrium_targets.csv"))
  cat("Outputs saved under: ", out_tables, " and ", out_figures, "\n", sep = "")
}
