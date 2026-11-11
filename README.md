# Subnational life expectancy at birth in India: Bayesian projections to 2050

R and Stan code for the paper "Subnational Life Expectancy at Birth in India: Historical Trends and Bayesian Projections to 2050".

The code fits a Bayesian hierarchical model to female and male life expectancy at birth for 23 Indian states, using Sample Registration System (SRS) estimates for 1971–2024, and forecasts it to 2050. The model is described in the paper and its supplementary material.

## Folders

```
data/      input panel (model_data.csv) and state key (groups.csv)
stan/      life_expectancy_model.stan
codes/     analysis scripts, run in numbered order
fits/      prepared data and Stan fits (not uploaded; created by the scripts)
tables/    output tables
figures/   output figures
```

## Running the analysis

The fitted models are not included in the repository because of their size. Running the scripts from the project root, in this order, generates them along with all tables and figures:

| Script | Purpose |
|---|---|
| `00_prepare_stan_data.R` | Prepares the model data |
| `01_run_models.R` | Prior predictive checks, then the main, tighter-prior and wider-prior fits |
| `02_compare_prior_sensitivity.R` | Compares the three prior settings |
| `03_temporal_holdout_validation_2000_2010.R` | Holdout validation with forecast origins in 2000 and 2010 |
| `04_threshold_probabilities_2050.R` | Probability that life expectancy is below 80 years in 2024 and 2050 |
| `05_extract_manuscript_results.R` | Collects the numbers reported in the paper into `tables/manuscript_results_raw/` |
| `06_holdout_combined_figure.R` | Combined holdout figure |

Scripts 01 and 03 reuse saved fits. After changing the data, priors or Stan code, delete the old fits in `fits/` (keep `prep_objects.rds`) and rerun.

## Requirements

R 4.6.0 with `rstan`, `tidyverse`, `here`, `bayesplot` and `scales`.

Figures use Times New Roman. Where that font is not installed, change `base_family` to `"serif"`.

## Data

Life expectancy estimates come from SRS publications of the Office of the Registrar General, India: the SRS compendium for 1970–75 to 2009–13, and abridged life tables calculated from SRS age-specific mortality rates for later years.

## Licence

PolyForm Noncommercial License 1.0.0. The code may be used, changed and shared for noncommercial purposes only. See [LICENSE](LICENSE) for the full terms.
