# Prepare model-ready state-year data for the life-expectancy model.
# Inputs are the existing panel and group key; no life-expectancy values are fabricated.
# Run from the project root: Rscript codes/00_prepare_stan_data.R

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
})

# ------------------------------ Settings --------------------------------------
MODEL_DATA_FILE <- here::here("data", "model_data.csv")
GROUPS_FILE <- here::here("data", "groups.csv")
PREP_FILE <- here::here("fits", "prep_objects.rds")
TABLE_DIR <- here::here("tables", "preparation")
FIGURE_DIR <- here::here("figures", "preparation")
MIN_FORECAST_START_YEAR <- 2025L
FORECAST_END_YEAR <- 2050L

# Existing preparation objects are backed up before rebuilding.
OVERWRITE_PREP_OBJECTS <- TRUE

EAG_PLUS_ASSAM <- c(
  "Assam", "Bihar", "Chhattisgarh", "Jharkhand", "Madhya Pradesh",
  "Odisha", "Rajasthan", "Uttar Pradesh", "Uttarakhand"
)

for (path in c(here::here("fits"), TABLE_DIR, FIGURE_DIR)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}
if (!file.exists(MODEL_DATA_FILE)) {
  stop("Model panel not found: ", MODEL_DATA_FILE,
       "\nExpected the existing file data/model_data.csv.", call. = FALSE)
}
if (!file.exists(GROUPS_FILE)) {
  stop("Group key not found: ", GROUPS_FILE,
       "\nExpected the existing file data/groups.csv.", call. = FALSE)
}
if (file.exists(PREP_FILE) && !OVERWRITE_PREP_OBJECTS) {
  stop("Preparation object already exists: ", PREP_FILE,
       "\nSet OVERWRITE_PREP_OBJECTS <- TRUE only when intentionally rebuilding it.", call. = FALSE)
}

# Known spelling variants only; other labels are left visible for validation.
canonical_state <- function(x) {
  original <- trimws(as.character(x))
  key <- tolower(gsub("[^[:alnum:]]", "", original))
  dplyr::case_when(
    key == "andhrapradesh" ~ "Andhra Pradesh",
    key == "assam" ~ "Assam",
    key == "bihar" ~ "Bihar",
    key %in% c("chhattisgarh", "chhatisgarh") ~ "Chhattisgarh",
    key %in% c("delhi", "nctofdelhi", "nationalcapitalterritoryofdelhi") ~ "Delhi",
    key == "goa" ~ "Goa",
    key == "gujarat" ~ "Gujarat",
    key == "haryana" ~ "Haryana",
    key == "himachalpradesh" ~ "Himachal Pradesh",
    key %in% c("jammukashmir", "jammuandkashmir") ~ "Jammu & Kashmir",
    key == "jharkhand" ~ "Jharkhand",
    key == "karnataka" ~ "Karnataka",
    key == "kerala" ~ "Kerala",
    key == "madhyapradesh" ~ "Madhya Pradesh",
    key == "maharashtra" ~ "Maharashtra",
    key %in% c("odisha", "orissa") ~ "Odisha",
    key == "punjab" ~ "Punjab",
    key == "rajasthan" ~ "Rajasthan",
    key == "tamilnadu" ~ "Tamil Nadu",
    key == "telangana" ~ "Telangana",
    key == "uttarpradesh" ~ "Uttar Pradesh",
    key %in% c("uttarakhand", "uttaranchal") ~ "Uttarakhand",
    key == "westbengal" ~ "West Bengal",
    TRUE ~ original
  )
}
normalise_key <- function(x) tolower(gsub("[^[:alnum:]]", "", trimws(as.character(x))))

# ------------------------------ Read inputs -----------------------------------
model_raw <- readr::read_csv(MODEL_DATA_FILE, show_col_types = FALSE, progress = FALSE)
groups_raw <- readr::read_csv(GROUPS_FILE, show_col_types = FALSE, progress = FALSE)

required_model <- c("g", "year", "e0_f", "e0_m")
required_groups <- c("g", "state")
missing_model <- setdiff(required_model, names(model_raw))
missing_groups <- setdiff(required_groups, names(groups_raw))
if (length(missing_model)) stop("Missing column(s) in model_data.csv: ", paste(missing_model, collapse = ", "), call. = FALSE)
if (length(missing_groups)) stop("Missing column(s) in groups.csv: ", paste(missing_groups, collapse = ", "), call. = FALSE)
if (nrow(model_raw) == 0L || nrow(groups_raw) == 0L) stop("Input file(s) are empty.", call. = FALSE)

# The group file defines the canonical mapping from original group IDs to states.
groups <- groups_raw %>%
  transmute(g_source = suppressWarnings(as.integer(g)), state = canonical_state(state)) %>%
  arrange(g_source)
if (anyNA(groups$g_source) || anyDuplicated(groups$g_source) || anyDuplicated(groups$state)) {
  stop("groups.csv must have unique integer g values and unique state names.", call. = FALSE)
}
if (!identical(groups$g_source, seq_len(nrow(groups)))) {
  stop("groups.csv g must be consecutive integers 1, ..., number of states.", call. = FALSE)
}
groups <- groups %>% mutate(g = row_number())
state_levels_full <- as.character(groups$state)

panel0 <- model_raw %>%
  transmute(
    g_source = suppressWarnings(as.integer(g)),
    year = suppressWarnings(as.numeric(year)),
    e0_f = suppressWarnings(as.numeric(e0_f)),
    e0_m = suppressWarnings(as.numeric(e0_m)),
    state_in_panel = if ("state" %in% names(model_raw)) as.character(model_raw$state) else NA_character_
  ) %>%
  left_join(groups %>% select(g_source, g, state), by = "g_source")

if (anyNA(panel0$g) || anyNA(panel0$state)) {
  bad_groups <- sort(unique(panel0$g_source[is.na(panel0$g)]))
  stop("model_data.csv contains g values absent from groups.csv: ", paste(bad_groups, collapse = ", "), call. = FALSE)
}
if ("state" %in% names(model_raw)) {
  panel_state <- canonical_state(panel0$state_in_panel)
  mismatch <- !is.na(panel_state) & normalise_key(panel_state) != normalise_key(panel0$state)
  if (any(mismatch)) {
    bad <- panel0 %>% filter(mismatch) %>% distinct(g_source, state_in_panel, state)
    readr::write_csv(bad, file.path(TABLE_DIR, "state_group_mismatches.csv"))
    stop("State labels in model_data.csv do not agree with groups.csv for some rows. See tables/preparation/state_group_mismatches.csv.", call. = FALSE)
  }
}

df <- panel0 %>%
  transmute(state, year, e0_f, e0_m, g) %>%
  arrange(g, year)

invalid_rows <- df %>%
  filter(!is.finite(year) | !is.finite(e0_f) | !is.finite(e0_m) | e0_f <= 0 | e0_m <= 0)
if (nrow(invalid_rows)) {
  readr::write_csv(invalid_rows, file.path(TABLE_DIR, "invalid_model_rows.csv"))
  stop(nrow(invalid_rows), " row(s) have missing/non-finite years or non-positive/non-finite life expectancy. See tables/preparation/invalid_model_rows.csv.", call. = FALSE)
}
if (anyDuplicated(df[c("g", "year")])) {
  duplicated_rows <- df %>% count(state, g, year) %>% filter(n > 1L)
  readr::write_csv(duplicated_rows, file.path(TABLE_DIR, "duplicate_state_year_rows.csv"))
  stop("Duplicate state-year rows detected. See tables/preparation/duplicate_state_year_rows.csv.", call. = FALSE)
}
if (!setequal(unique(df$g), seq_along(state_levels_full))) {
  stop("At least one state in groups.csv has no panel observations, or an invalid group ID exists.", call. = FALSE)
}

# Rebuild the time index after sorting. This avoids relying on stale index columns
# in model_data.csv and guarantees prev_idx points to the previous global row.
df <- df %>%
  arrange(g, year) %>%
  mutate(global_idx = row_number()) %>%
  group_by(g) %>%
  mutate(
    row_in_group = row_number(),
    first_year_g = min(year),
    last_year_g = max(year),
    is_first = as.integer(row_in_group == 1L),
    prev_idx = if_else(row_in_group == 1L, 0L, as.integer(lag(global_idx)), missing = 0L),
    prev_global_idx = prev_idx,
    delta_t = if_else(row_in_group == 1L, 0, as.numeric(year - lag(year)), missing = 0),
    t_since_start = as.numeric(year - first_year_g),
    is_2020 = as.integer(abs(year - 2020) < 1e-8),
    is_2021 = as.integer(abs(year - 2021) < 1e-8),
    is_post2021 = as.integer(year > 2021),
    years_since_2021 = as.numeric(pmax(year - 2021, 0))
  ) %>%
  ungroup() %>%
  select(-global_idx, -row_in_group)

if (any(df$delta_t[df$is_first == 0L] <= 0)) stop("Non-positive time gap within a state.", call. = FALSE)
if (sum(df$is_first) != nrow(groups)) stop("Expected exactly one first observation per state.", call. = FALSE)

first_idx <- which(df$is_first == 1L)
last_idx <- vapply(seq_len(nrow(groups)), function(gi) max(which(df$g == gi)), integer(1))
first_year <- as.numeric(df$year[first_idx])
last_year <- as.numeric(df$year[last_idx])

forecast_start_year <- max(as.integer(MIN_FORECAST_START_YEAR), floor(max(df$year)) + 1L)
if (forecast_start_year > FORECAST_END_YEAR) {
  stop("Forecast start year ", forecast_start_year, " exceeds forecast end year ", FORECAST_END_YEAR, ".", call. = FALSE)
}
forecast_years_full <- seq.int(forecast_start_year, as.integer(FORECAST_END_YEAR))
eag_key <- normalise_key(EAG_PLUS_ASSAM)
# Model data convention: 1 = EAG + Assam, 2 = non-EAG.
eag_group <- if_else(normalise_key(state_levels_full) %in% eag_key, 1L, 2L)

stan_data_full <- list(
  N = as.integer(nrow(df)),
  g = as.integer(df$g),
  is_first = as.integer(df$is_first),
  prev_idx = as.integer(df$prev_idx),
  delta_t = as.numeric(df$delta_t),
  t_since_start = as.numeric(df$t_since_start),
  y_f = as.numeric(df$e0_f),
  y_m = as.numeric(df$e0_m),
  is_2020 = as.integer(df$is_2020),
  is_2021 = as.integer(df$is_2021),
  is_post2021 = as.integer(df$is_post2021),
  years_since_2021 = as.numeric(df$years_since_2021),
  N_groups = as.integer(nrow(groups)),
  N_forecast = as.integer(length(forecast_years_full)),
  forecast_years = as.numeric(forecast_years_full),
  first_year = first_year,
  last_year = last_year,
  last_idx = as.integer(last_idx),
  eag_group = as.integer(eag_group)
)

prep <- list(
  config = list(
    eag_states = EAG_PLUS_ASSAM,
    holdout_start_years = c(2000L, 2010L),
    holdout_end_year = 2019L,
    forecast_start_year = forecast_start_year,
    forecast_end_year = as.integer(FORECAST_END_YEAR)
  ),
  df = df,
  df_groups = groups %>% select(g, state),
  state_levels_full = state_levels_full,
  eag_group = eag_group,
  forecast_years_full = as.numeric(forecast_years_full),
  stan_data_full = stan_data_full,
  source_files = c(
    model_data = normalizePath(MODEL_DATA_FILE, winslash = "/", mustWork = TRUE),
    groups = normalizePath(GROUPS_FILE, winslash = "/", mustWork = TRUE)
  ),
  prepared_at = as.character(Sys.time())
)

# Audit table and state-specific time coverage.
state_coverage <- df %>%
  group_by(state, g) %>%
  summarise(n_observations = n(), first_year = min(year), last_year = max(year),
            n_gaps_over_1_year = sum(diff(year) > 1), .groups = "drop") %>%
  arrange(g)
readr::write_csv(state_coverage, file.path(TABLE_DIR, "state_year_coverage.csv"))
readr::write_csv(df, file.path(TABLE_DIR, "model_panel_used.csv"))
readr::write_csv(tibble(
  item = c("model_data.csv rows read", "groups.csv states", "Rows retained", "First observation year",
           "Last observation year", "Forecast start year", "Forecast end year"),
  value = c(nrow(model_raw), nrow(groups), nrow(df), min(df$year), max(df$year),
            min(forecast_years_full), max(forecast_years_full))
), file.path(TABLE_DIR, "preparation_audit.csv"))

# Save a title-free overview plot for inspection or later use in a manuscript.
plot_data <- df %>%
  select(state, year, e0_f, e0_m) %>%
  pivot_longer(c(e0_f, e0_m), names_to = "sex", values_to = "life_expectancy") %>%
  mutate(sex = recode(sex, e0_f = "Female", e0_m = "Male"),
         state = factor(state, levels = rev(state_levels_full)))
p_history <- ggplot(plot_data, aes(year, life_expectancy, colour = sex, group = sex)) +
  geom_line(linewidth = 0.4, alpha = 0.8) +
  geom_point(size = 0.9, alpha = 0.9) +
  facet_wrap(~state, ncol = 5) +
  scale_colour_manual(values = c(Female = "#D55E00", Male = "#0072B2")) +
  labs(title = NULL, subtitle = NULL, x = "Reference year", y = "Life expectancy at birth (years)", colour = NULL) +
  theme_classic(base_size = 9) +
  theme(plot.title = element_blank(), plot.subtitle = element_blank(),
        strip.background = element_blank(), strip.text = element_text(size = 8),
        panel.grid = element_blank(), panel.spacing = grid::unit(0.35, "lines"),
        axis.text = element_text(colour = "black", size = 7), legend.position = "bottom")
ggsave(file.path(FIGURE_DIR, "observed_life_expectancy_by_state.png"), p_history,
       width = 16, height = 13, dpi = 600, bg = "white")
ggsave(file.path(FIGURE_DIR, "observed_life_expectancy_by_state.pdf"), p_history,
       width = 16, height = 13, bg = "white")

# If replacing a previous preparation object, retain a timestamped backup.
if (file.exists(PREP_FILE)) {
  backup_file <- file.path(dirname(PREP_FILE), paste0("prep_objects_backup_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds"))
  file.copy(PREP_FILE, backup_file, overwrite = FALSE)
  message("Existing prep object backed up to: ", backup_file)
}
saveRDS(prep, PREP_FILE)
cat("\nPREPARATION COMPLETE\n")
cat("Input panel: ", normalizePath(MODEL_DATA_FILE, winslash = "/", mustWork = TRUE), "\n", sep = "")
cat("Group key: ", normalizePath(GROUPS_FILE, winslash = "/", mustWork = TRUE), "\n", sep = "")
cat("Saved preparation object: ", PREP_FILE, "\n", sep = "")
cat("Rows/states: ", nrow(df), " / ", nrow(groups), "\n", sep = "")
cat("Observed years: ", min(df$year), " to ", max(df$year), "\n", sep = "")
cat("Forecast years: ", min(forecast_years_full), " to ", max(forecast_years_full), "\n", sep = "")
cat("Tables: ", TABLE_DIR, " | Figures: ", FIGURE_DIR, "\n", sep = "")
