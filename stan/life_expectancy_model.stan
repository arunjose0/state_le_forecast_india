// life_expectancy_model_v2.stan
//
// State-level e0 model: asymptotic trend, bivariate OU deviations, COVID shock.
// v2 changes only the log_kappa_f prior (median half-life 50 years) and the
// tighter-set C_f and C_m priors, to match Supplementary Tables ST1 and ST2.
//
// prior_mode: 2 main, 3 tighter, 4 wider. Set prior_only = 1 for a prior run.

functions {
  matrix stabilize_covariance(matrix cov) {
    matrix[2, 2] sym_cov;
    real covariance_scale;
    sym_cov = 0.5 * (cov + cov');
    covariance_scale = fmax(1.0, max(diagonal(sym_cov)));
    return sym_cov + diag_matrix(rep_vector(1e-8 * covariance_scale, 2));
  }

  real ou_diffusion_multiplier(real rate, real delta_t) {
    if (rate * delta_t < 1e-8) {
      return sqrt(delta_t);
    }
    return sqrt(fmax(0, -expm1(-2 * rate * delta_t) / (2 * rate)));
  }

  real covid_adjustment(
    real year_value,
    real shock2020,
    real shock2021,
    real phi
  ) {
    if (abs(year_value - 2020) < 1e-8) return shock2020;
    if (abs(year_value - 2021) < 1e-8) return shock2021;
    if (year_value > 2021) return shock2021 * pow(phi, year_value - 2021);
    return 0;
  }
}

data {
  int<lower=1> N;
  array[N] int<lower=1> g;
  array[N] int<lower=0, upper=1> is_first;
  array[N] int<lower=0> prev_idx;
  vector<lower=0>[N] delta_t;
  vector[N] y_f;
  vector[N] y_m;
  int<lower=1> N_groups;
  int<lower=1> N_tiers;
  array[N_groups] int<lower=1, upper=N_tiers> tier_id;

  array[N] int<lower=0, upper=1> is_2020;
  array[N] int<lower=0, upper=1> is_2021;
  array[N] int<lower=0, upper=1> is_post2021;
  vector<lower=0>[N] years_since_2021;

  int<lower=0> N_forecast;
  vector[N_forecast] forecast_years;
  vector[N_groups] first_year;
  vector[N_groups] last_year;
  array[N_groups] int<lower=1, upper=N> last_idx;

  int<lower=2, upper=4> prior_mode;
  int<lower=0, upper=1> prior_only;
  array[N_groups] int<lower=1, upper=2> eag_group; // unused; kept for data compatibility
  int<lower=0, upper=1> compute_yrep;
}

parameters {
  real C_f;
  real C_m;
  vector[N_tiers] tier_effect_raw_f;
  vector[N_tiers] tier_effect_raw_m;

  vector<lower=0>[2] tau_C;
  cholesky_factor_corr[2] L_Omega_C;
  matrix[2, N_groups] z_C;

  real log_kappa_f;
  real delta_log_kappa;

  real log_lambda_gap;
  vector<lower=0>[2] sigma_dev;
  cholesky_factor_corr[2] L_Omega_dev;

  real<lower=0> sigma_obs;

  real delta2020_f;
  real delta2021_f;
  real delta2020_m;
  real delta2021_m;
  real<lower=0, upper=1> phi_covid_f;
  real<lower=0, upper=1> phi_covid_m;
}

transformed parameters {
  real delta_C = C_m - C_f;
  vector[N_tiers] tier_effect_f = tier_effect_raw_f - mean(tier_effect_raw_f);
  vector[N_tiers] tier_effect_m = tier_effect_raw_m - mean(tier_effect_raw_m);
  vector[N_tiers] C_tier_f = C_f + tier_effect_f;
  vector[N_tiers] C_tier_m = C_m + tier_effect_m;
  vector[N_groups] C_state_f;
  vector[N_groups] C_state_m;
  matrix[2, 2] L_C_scaled = diag_pre_multiply(tau_C, L_Omega_C);
  real kappa_f = exp(log_kappa_f);
  real kappa_m = kappa_f * exp(delta_log_kappa);
  real lambda_dev = fmax(kappa_f, kappa_m) + exp(log_lambda_gap);
  real deviation_half_life = log(2) / lambda_dev;
  real long_run_sd_dev_f = sigma_dev[1] / sqrt(2 * lambda_dev);
  real long_run_sd_dev_m = sigma_dev[2] / sqrt(2 * lambda_dev);

  vector[N_groups] baseline_f;
  vector[N_groups] baseline_m;
  vector[N] elapsed_since_first;
  vector[N] covid_shock_f;
  vector[N] covid_shock_m;
  vector[N] mu_curve_f;
  vector[N] mu_curve_m;

  array[N] vector[2] u_pred_mean;
  array[N] matrix[2, 2] u_pred_cov;
  array[N] vector[2] u_filter_mean;
  array[N] matrix[2, 2] u_filter_cov;
  vector[N] marginal_loglik;

  matrix[2, 2] I2 = diag_matrix(rep_vector(1.0, 2));
  matrix[2, 2] L_scaled = diag_pre_multiply(sigma_dev, L_Omega_dev);

  for (gi in 1:N_groups) {
    vector[2] z_target = col(z_C, gi);
    vector[2] target_deviation = L_C_scaled * z_target;
    C_state_f[gi] = C_tier_f[tier_id[gi]] + target_deviation[1];
    C_state_m[gi] = C_tier_m[tier_id[gi]] + target_deviation[2];
  }

  for (i in 1:N) {
    covid_shock_f[i] = 0;
    covid_shock_m[i] = 0;
    if (is_2020[i] == 1) {
      covid_shock_f[i] = delta2020_f;
      covid_shock_m[i] = delta2020_m;
    } else if (is_2021[i] == 1) {
      covid_shock_f[i] = delta2021_f;
      covid_shock_m[i] = delta2021_m;
    } else if (is_post2021[i] == 1) {
      covid_shock_f[i] = delta2021_f * pow(phi_covid_f, years_since_2021[i]);
      covid_shock_m[i] = delta2021_m * pow(phi_covid_m, years_since_2021[i]);
    }

    if (is_first[i] == 1) {
      baseline_f[g[i]] = y_f[i] - covid_shock_f[i];
      baseline_m[g[i]] = y_m[i] - covid_shock_m[i];
      elapsed_since_first[i] = 0;
    } else {
      elapsed_since_first[i] = elapsed_since_first[prev_idx[i]] + delta_t[i];
    }
  }

  for (i in 1:N) {
    real rho_trend_f = exp(-kappa_f * elapsed_since_first[i]);
    real rho_trend_m = exp(-kappa_m * elapsed_since_first[i]);
    mu_curve_f[i] = C_state_f[g[i]] + (baseline_f[g[i]] - C_state_f[g[i]]) * rho_trend_f + covid_shock_f[i];
    mu_curve_m[i] = C_state_m[g[i]] + (baseline_m[g[i]] - C_state_m[g[i]]) * rho_trend_m + covid_shock_m[i];
  }

  // Kalman filter; u = 0 at each group's first observation
  for (i in 1:N) {
    if (is_first[i] == 1) {
      u_pred_mean[i] = rep_vector(0.0, 2);
      u_filter_mean[i] = rep_vector(0.0, 2);
      u_pred_cov[i] = rep_matrix(0.0, 2, 2);
      u_filter_cov[i] = rep_matrix(0.0, 2, 2);
      marginal_loglik[i] = 0;
    } else {
      int pi = prev_idx[i];
      real dtt = fmax(delta_t[i], 1e-6);
      real rho_dev = exp(-lambda_dev * dtt);
      real mult_dev = ou_diffusion_multiplier(lambda_dev, dtt);
      matrix[2, 2] Q = square(mult_dev) * (L_scaled * L_scaled');
      matrix[2, 2] S;
      matrix[2, 2] K;
      matrix[2, 2] IK;
      vector[2] residual;
      vector[2] innovation;

      u_pred_mean[i] = rho_dev * u_filter_mean[pi];
      u_pred_cov[i] = square(rho_dev) * u_filter_cov[pi] + Q;
      S = u_pred_cov[i] + square(sigma_obs) * I2;
      residual[1] = y_f[i] - mu_curve_f[i];
      residual[2] = y_m[i] - mu_curve_m[i];
      innovation = residual - u_pred_mean[i];
      marginal_loglik[i] = multi_normal_lpdf(residual | u_pred_mean[i], S);

      if (prior_only == 1) {
        u_filter_mean[i] = u_pred_mean[i];
        u_filter_cov[i] = u_pred_cov[i];
      } else {
        K = u_pred_cov[i] * mdivide_left_spd(S, I2);
        u_filter_mean[i] = u_pred_mean[i] + K * innovation;
        IK = I2 - K;
        u_filter_cov[i] = IK * u_pred_cov[i] * IK' + square(sigma_obs) * K * K';
        u_filter_cov[i] = 0.5 * (u_filter_cov[i] + u_filter_cov[i]');
      }
    }
  }
}

model {
  if (prior_mode == 2) {
    // main (ST1)
    C_f ~ normal(87, 2);
    C_m ~ normal(82, 2);
    tier_effect_raw_f ~ normal(0, 4);
    tier_effect_raw_m ~ normal(0, 4);
    tau_C ~ normal(0, 3);
    L_Omega_C ~ lkj_corr_cholesky(2);
    log_kappa_f ~ normal(log(log(2) / 50), log(4) / 1.96);
    delta_log_kappa ~ normal(0, 0.5);
    log_lambda_gap ~ normal(log(0.20), 0.7);
    sigma_dev ~ student_t(3, 0, 1.5);
    L_Omega_dev ~ lkj_corr_cholesky(2);
    sigma_obs ~ student_t(3, 0, 1.5);
    delta2020_f ~ normal(0, 5);
    delta2021_f ~ normal(0, 5);
    delta2020_m ~ normal(0, 5);
    delta2021_m ~ normal(0, 5);
    phi_covid_f ~ beta(2, 2);
    phi_covid_m ~ beta(2, 2);
  } else if (prior_mode == 3) {
    // tighter (ST2)
    C_f ~ normal(87, 1.5);
    C_m ~ normal(82, 1.5);
    tier_effect_raw_f ~ normal(0, 2);
    tier_effect_raw_m ~ normal(0, 2);
    tau_C ~ normal(0, 2);
    L_Omega_C ~ lkj_corr_cholesky(2);
    log_kappa_f ~ normal(log(log(2) / 50), log(2) / 1.96);
    delta_log_kappa ~ normal(0, 0.35);
    log_lambda_gap ~ normal(log(0.05), 0.55);
    sigma_dev ~ normal(0, 1.0);
    L_Omega_dev ~ lkj_corr_cholesky(2);
    sigma_obs ~ normal(0, 0.75);
    delta2020_f ~ normal(0, 3);
    delta2021_f ~ normal(0, 3);
    delta2020_m ~ normal(0, 3);
    delta2021_m ~ normal(0, 3);
    phi_covid_f ~ beta(2, 2);
    phi_covid_m ~ beta(2, 2);
  } else {
    // wider (ST2)
    C_f ~ normal(87, 5);
    C_m ~ normal(82, 5);
    tier_effect_raw_f ~ normal(0, 8);
    tier_effect_raw_m ~ normal(0, 8);
    tau_C ~ normal(0, 5);
    L_Omega_C ~ lkj_corr_cholesky(1);
    log_kappa_f ~ normal(log(log(2) / 50), log(8) / 1.96);
    delta_log_kappa ~ normal(0, 1.0);
    log_lambda_gap ~ normal(log(0.20), 1.0);
    sigma_dev ~ student_t(3, 0, 3.0);
    L_Omega_dev ~ lkj_corr_cholesky(1);
    sigma_obs ~ student_t(3, 0, 3.0);
    delta2020_f ~ normal(0, 8);
    delta2021_f ~ normal(0, 8);
    delta2020_m ~ normal(0, 8);
    delta2021_m ~ normal(0, 8);
    phi_covid_f ~ beta(1, 1);
    phi_covid_m ~ beta(1, 1);
  }

  to_vector(z_C) ~ std_normal();

  if (prior_only == 0) {
    for (i in 1:N) {
      if (is_first[i] == 0) target += marginal_loglik[i];
    }
  }
}

generated quantities {
  array[N] real y_rep_f;
  array[N] real y_rep_m;
  array[N_groups, N_forecast] real trend_forecast_f;
  array[N_groups, N_forecast] real trend_forecast_m;
  array[N_groups, N_forecast] real mean_forecast_f;
  array[N_groups, N_forecast] real mean_forecast_m;
  array[N_groups, N_forecast] real latent_forecast_f;
  array[N_groups, N_forecast] real latent_forecast_m;
  array[N_groups, N_forecast] real e0_forecast_f;
  array[N_groups, N_forecast] real e0_forecast_m;

  matrix[2, 2] Omega_dev = multiply_lower_tri_self_transpose(L_Omega_dev);
  matrix[2, 2] Omega_C = multiply_lower_tri_self_transpose(L_Omega_C);
  real rho_dev_corr = Omega_dev[1, 2];
  real rho_C_corr = Omega_C[1, 2];
  real kappa_ratio_m_to_f = exp(delta_log_kappa);
  real convergence_half_life_f = log(2) / kappa_f;
  real convergence_half_life_m = log(2) / kappa_m;
  real years_to_90pct_f = log(10) / kappa_f;
  real years_to_90pct_m = log(10) / kappa_m;

  // one-step-ahead predictive replicates
  for (i in 1:N) {
    if (is_first[i] == 1) {
      y_rep_f[i] = y_f[i];
      y_rep_m[i] = y_m[i];
    } else if (compute_yrep == 1) {
      vector[2] u_draw = multi_normal_rng(
        u_pred_mean[i], stabilize_covariance(u_pred_cov[i])
      );
      y_rep_f[i] = mu_curve_f[i] + u_draw[1] + normal_rng(0, sigma_obs);
      y_rep_m[i] = mu_curve_m[i] + u_draw[2] + normal_rng(0, sigma_obs);
    } else {
      y_rep_f[i] = mu_curve_f[i] + u_pred_mean[i][1];
      y_rep_m[i] = mu_curve_m[i] + u_pred_mean[i][2];
    }
  }

  // forecasts from the filtered last latent state
  for (gi in 1:N_groups) {
    real first_de_shocked_f = baseline_f[gi];
    real first_de_shocked_m = baseline_m[gi];
    real prev_year = last_year[gi];
    vector[2] u_last;
    real prev_u_f;
    real prev_u_m;
    if (is_first[last_idx[gi]] == 1) {
      u_last = rep_vector(0.0, 2);
    } else {
      u_last = multi_normal_rng(
        u_filter_mean[last_idx[gi]],
        stabilize_covariance(u_filter_cov[last_idx[gi]])
      );
    }
    prev_u_f = u_last[1];
    prev_u_m = u_last[2];

    for (t in 1:N_forecast) {
      real target_year = forecast_years[t];
      real dtt = fmax(target_year - prev_year, 1e-6);
      real elapsed = fmax(target_year - first_year[gi], 0);
      real rho_trend_f = exp(-kappa_f * elapsed);
      real rho_trend_m = exp(-kappa_m * elapsed);
      real rho_dev = exp(-lambda_dev * dtt);
      real mult_dev = ou_diffusion_multiplier(lambda_dev, dtt);
      real covid_f = covid_adjustment(target_year, delta2020_f, delta2021_f, phi_covid_f);
      real covid_m = covid_adjustment(target_year, delta2020_m, delta2021_m, phi_covid_m);
      real trend_f = C_state_f[gi] + (first_de_shocked_f - C_state_f[gi]) * rho_trend_f + covid_f;
      real trend_m = C_state_m[gi] + (first_de_shocked_m - C_state_m[gi]) * rho_trend_m + covid_m;
      real expected_u_f = rho_dev * prev_u_f;
      real expected_u_m = rho_dev * prev_u_m;
      vector[2] z_new;
      vector[2] z_corr;
      real next_u_f;
      real next_u_m;

      z_new[1] = normal_rng(0, 1);
      z_new[2] = normal_rng(0, 1);
      z_corr = L_Omega_dev * z_new;
      next_u_f = expected_u_f + sigma_dev[1] * mult_dev * z_corr[1];
      next_u_m = expected_u_m + sigma_dev[2] * mult_dev * z_corr[2];

      trend_forecast_f[gi, t] = trend_f;
      trend_forecast_m[gi, t] = trend_m;
      mean_forecast_f[gi, t] = trend_f + expected_u_f;
      mean_forecast_m[gi, t] = trend_m + expected_u_m;
      latent_forecast_f[gi, t] = trend_f + next_u_f;
      latent_forecast_m[gi, t] = trend_m + next_u_m;
      e0_forecast_f[gi, t] = latent_forecast_f[gi, t] + normal_rng(0, sigma_obs);
      e0_forecast_m[gi, t] = latent_forecast_m[gi, t] + normal_rng(0, sigma_obs);

      prev_u_f = next_u_f;
      prev_u_m = next_u_m;
      prev_year = target_year;
    }
  }
}
