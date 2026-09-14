#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(openxlsx)
  library(dplyr)
  library(purrr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[[1]] else "rice_lineage_support_tables.xlsx"
output_file <- if (length(args) >= 2) args[[2]] else "rice_founder_lineage_metrics.xlsx"
rho_max <- if (length(args) >= 3) as.numeric(args[[3]]) else 0.15
min_clade_support <- if (length(args) >= 4) as.numeric(args[[4]]) else 2

stopifnot(file.exists(input_file))

mapping_df <- read.xlsx(input_file, sheet = "rice_mapping") %>%
  mutate(across(everything(), as.character))

clade_df <- read.xlsx(input_file, sheet = "rice_selected_clades") %>%
  mutate(across(everything(), as.character)) %>%
  mutate(
    size = as.integer(size),
    support_n = as.numeric(support_n),
    descendant_list = strsplit(descendant_samples, ";", fixed = TRUE)
  )

get_top_level_lineages <- function(seed_group, mapping_df, clade_df) {
  tips <- mapping_df %>%
    filter(seed_group == !!seed_group) %>%
    pull(new_id)

  seed_clades <- clade_df %>%
    mutate(
      within_seed = map_lgl(descendant_list, ~ all(.x %in% tips)),
      proper_subset = map_int(descendant_list, length) < length(tips),
      enough_support = support_n >= min_clade_support
    ) %>%
    filter(within_seed, proper_subset, enough_support)

  keep <- rep(TRUE, nrow(seed_clades))

  if (nrow(seed_clades) > 0) {
    for (i in seq_len(nrow(seed_clades))) {
      for (j in seq_len(nrow(seed_clades))) {
        if (i == j) next
        a <- seed_clades$descendant_list[[i]]
        b <- seed_clades$descendant_list[[j]]
        if (length(a) < length(b) && all(a %in% b)) {
          keep[i] <- FALSE
        }
      }
    }
  }

  maximal_clades <- seed_clades[keep, , drop = FALSE]
  covered_tips <- sort(unique(unlist(maximal_clades$descendant_list, use.names = FALSE)))
  singleton_tips <- setdiff(tips, covered_tips)

  lineage_sizes <- c(lengths(maximal_clades$descendant_list), rep(1L, length(singleton_tips)))
  lineage_sizes <- lineage_sizes[lineage_sizes > 0]

  list(
    tips = tips,
    maximal_clades = maximal_clades,
    singleton_tips = singleton_tips,
    lineage_sizes = lineage_sizes
  )
}

estimate_seed_metrics <- function(seed_group, mapping_df, clade_df, rho_max = 0.15) {
  lineage_info <- get_top_level_lineages(seed_group, mapping_df, clade_df)
  lineage_sizes <- lineage_info$lineage_sizes

  sample_n <- length(lineage_info$tips)
  L_obs <- length(lineage_sizes)
  f1 <- sum(lineage_sizes == 1)
  f2 <- sum(lineage_sizes == 2)
  coverage <- if (sample_n > 0) 1 - f1 / sample_n else NA_real_
  Chao1 <- if (f2 > 0) {
    L_obs + f1 * f1 / (2 * f2)
  } else {
    L_obs + f1 * (f1 - 1) / 2
  }

  over_split_lower <- Chao1 * (1 - rho_max)
  over_split_upper <- Chao1
  unresolved_founder_note <- ifelse(f1 > 0, "possible_under_resolution", "no_singleton_signal")

  tibble(
    seed_group = seed_group,
    sample_n = sample_n,
    L_obs = L_obs,
    clade_children = nrow(lineage_info$maximal_clades),
    singleton_children = length(lineage_info$singleton_tips),
    lineage_sizes = paste(sort(lineage_sizes, decreasing = TRUE), collapse = ","),
    f1 = f1,
    f2 = f2,
    coverage = round(coverage, 3),
    Chao1 = round(Chao1, 2),
    min_clade_support = min_clade_support,
    sensitivity_rho_max = rho_max,
    over_split_lower = round(over_split_lower, 2),
    over_split_upper = round(over_split_upper, 2),
    unresolved_founder_note = unresolved_founder_note,
    reliability = case_when(
      sample_n >= 25 & coverage >= 0.85 ~ "high",
      sample_n >= 15 & coverage >= 0.60 ~ "moderate",
      TRUE ~ "low"
    )
  )
}

seed_metrics <- bind_rows(lapply(unique(mapping_df$seed_group), estimate_seed_metrics,
                                 mapping_df = mapping_df, clade_df = clade_df, rho_max = rho_max)) %>%
  arrange(desc(sample_n), seed_group)

anchor_seed_metrics <- seed_metrics %>%
  filter(reliability %in% c("high", "moderate"))

sensitivity_scenarios <- tibble(
  scenario = c(
    "equal_contribution_occupancy_baseline",
    "dirichlet_multinomial_unequal_contribution",
    "over_splitting_sensitivity",
    "under_resolution_merging_sensitivity",
    "include_vs_exclude_low_sampled_seeds"
  ),
  role = c(
    "baseline_interpretive_framework",
    "publication_grade_extension_not_fitted",
    "supporting_sensitivity_only",
    "supporting_sensitivity_only",
    "robustness_diagnostic"
  ),
  implementation_in_current_output = c(
    "empirical_L_obs_plus_completeness_diagnostics",
    "conceptual_only_no_fitted_alpha",
    paste0("reported as Chao1*(1-rho), rho_max=", rho_max),
    "reported qualitatively via unresolved_founder_note",
    paste0("anchor seeds retained as reliability high/moderate: ", paste(anchor_seed_metrics$seed_group, collapse = ", "))
  ),
  interpretation_note = c(
    "L_obs is the primary empirical readout under finite founder occupancy",
    "unequal founder contribution is biologically plausible but not fitted here",
    "scenario-based and not a fitted final corrected value",
    "true founder diversity may be underestimated if distinguishable shared mutations are absent",
    "low-sampled seeds support multicellular origin but are not used alone for strict quantification"
  )
)

model_notes <- tibble(
  item = c(
    "model_status",
    "target_quantity",
    "L_obs_role",
    "coverage_role",
    "Chao1_role",
    "observation_layer",
    "recommended_working_range"
  ),
  value = c(
    "publication_grade_probabilistic_framework_not_fully_fitted",
    "effective_founder_cell_diversity_under_current_regeneration_system",
    "primary_empirical_readout",
    "sampling_completeness_diagnostic",
    "auxiliary_richness_extrapolation",
    "finite_sampling_plus_split_merge_plus_genotype_calling_uncertainty",
    "approximately_10_to_15_cells_under_current_calling_and_lineage_definition_framework"
  )
)

summary_df <- seed_metrics %>%
  summarise(
    seed_n = n(),
    sample_total = sum(sample_n),
    median_L_obs = median(L_obs),
    max_L_obs = max(L_obs),
    best_supported_seed = seed_group[which.max(sample_n)],
    best_supported_L_obs = L_obs[which.max(sample_n)],
    best_supported_Chao1 = Chao1[which.max(sample_n)],
    anchor_seed_set = paste(anchor_seed_metrics$seed_group, collapse = ", "),
    anchor_L_obs_range = paste(range(anchor_seed_metrics$L_obs), collapse = "-"),
    high_reliability_seed_n = sum(reliability == "high"),
    moderate_reliability_seed_n = sum(reliability == "moderate")
  )

write.xlsx(
  list(
    per_seed_metrics = seed_metrics,
    anchor_seed_metrics = anchor_seed_metrics,
    sensitivity_scenarios = sensitivity_scenarios,
    model_notes = model_notes,
    summary = summary_df
  ),
  file = output_file,
  overwrite = TRUE
)

write.table(
  seed_metrics,
  file = sub("\\.xlsx$", ".tsv", output_file),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

message("Wrote founder-lineage metrics to: ", normalizePath(output_file, winslash = "/", mustWork = FALSE))
