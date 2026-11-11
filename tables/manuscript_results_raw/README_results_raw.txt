MANUSCRIPT RESULTS RAW-NUMBER EXTRACTION
Project root: /home/arunjose/subnational_india_life_expectancy_forecast
Output directory: /home/arunjose/subnational_india_life_expectancy_forecast/tables/manuscript_results_raw
No posterior model was fitted or refitted by this script.

Main output files:
- Historical: historical_extremes_and_ranges_1971_2024.csv; historical_health_index_group_summary_1971_2024.csv; historical_state_values_1971_2024.csv; historical_first_latest_values_and_sex_gaps_by_state.csv; historical_sex_gap_counts.csv
- Projections: latent_projection_estimates_by_state_2030_2040_2050.csv; latent_projection_estimates_2050_by_state.csv; projection_2050_summary_by_health_index_group.csv; interstate_dispersion_2024_observed_vs_2050_projected.csv; projected_female_minus_male_gap_2050_by_state.csv
- Equilibrium: state_equilibrium_targets_and_gaps.csv; equilibrium_target_summary_by_health_index_group.csv; health_index_tier_equilibrium_centres.csv when available
- Threshold probability: threshold_probability_below_80_by_state_2024_2050.csv; threshold_probability_change_2024_to_2050_by_state.csv; threshold_probability_below_80_by_health_index_group.csv; threshold_probability_change_by_health_index_group.csv
- COVID: covid_shock_and_decay_posterior_summary.csv
- Validation: holdout_metrics_overall_all_models.csv; holdout_metrics_by_origin_sex_lead.csv; holdout_metrics_by_state.csv when available; holdout predictions and diagnostics when available
- Diagnostics/PPC: posterior_sampler_diagnostic_status_by_model.csv; generated_quantities_status_by_model.csv; posterior_parameter_rhat_ess_all_models.csv; posterior_predictive_coverage_overall_by_sex.csv; posterior_predictive_coverage_by_state_and_sex.csv when available
- Prior checks: prior predictive files and diagnostic status, where available
- Sensitivity: prior_sensitivity_comparison_2050.csv and summary tables, where available
- key_numbers_for_results_review.csv is a compact roll-up. Consult the fuller tables for complete raw estimates.

Definitions and cautions:
1. Forecast projections and threshold probabilities refer to the latent trajectory, not future observation draws.
2. 2024 threshold probabilities for states lacking a 2024 observation are model-based latent-state estimates propagated from the latest observation.
3. Health Index group means are state-level/equal-weighted summaries, not population-weighted pooled life expectancy.
4. 2024 observed interstate ranges use only states with actual 2024 observations; 2050 ranges use states with 2050 projections. Check n_states when comparing ranges.
5. Target-minus-latest-observation intervals translate uncertainty in the target while treating the latest observed SRS value as fixed; input SRS standard errors are not modelled.
6. Point-only holdout benchmarks have no CRPS or predictive interval coverage; these cells should remain NA.
7. Convergence diagnostics should be interpreted separately from predictive validation and posterior predictive checks.
8. Missing optional files are listed as exists=FALSE in input_file_inventory.csv; do not interpret missing outputs as passed checks.
