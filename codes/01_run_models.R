# Run the four analysis steps in order: prior check, main fit, tighter sensitivity, wide sensitivity.
# Run from the project root after codes/00_prepare_stan_data.R has completed.

for (run_mode in 1:4) {
  RUN_MODE <- as.integer(run_mode)
  message("\n==================== RUN MODE ", RUN_MODE, " ====================")
  source(here::here("codes", "_run_model_engine.R"), local = FALSE)
}
