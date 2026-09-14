#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
  library(openxlsx)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ape)
  library(ggplot2)
  library(ggtree)
  library(ggnewscale)
  library(scatterpie)
})

input_mutation_file <- file.path("..", "mutation_id.xlsx")
input_mapping_file <- file.path("..", "sample_id_mapping.xlsx")
output_dir <- "."

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

species_cfg <- list(
  rice = list(sheet = "rice", exclude_ids = character(0), exclude_pattern = NULL),
  Liriodendron = list(sheet = "Liriodendron", exclude_ids = character(0), exclude_pattern = "[Cc]"),
  Arabidopsis = list(sheet = "Arabidopsis", exclude_ids = c("col-H-D12-1"), exclude_pattern = NULL)
)

read_mapping_sheet <- function(sheet_name) {
  read_xlsx(input_mapping_file, sheet = sheet_name) %>%
    mutate(across(everything(), as.character))
}

read_mutation_sheet <- function(sheet_name, exclude_ids = character(0), exclude_pattern = NULL) {
  raw_df <- read_xlsx(input_mutation_file, sheet = sheet_name)
  stopifnot(ncol(raw_df) >= 3)

  cleaned <- raw_df %>%
    transmute(
      chrom = as.character(.data[[names(raw_df)[1]]]),
      pos = as.character(.data[[names(raw_df)[2]]]),
      mutation_key = paste(chrom, pos, sep = ":"),
      id_raw = as.character(.data[[names(raw_df)[3]]])
    ) %>%
    separate_longer_delim(id_raw, delim = ";") %>%
    mutate(sample_id = str_trim(id_raw)) %>%
    filter(!is.na(sample_id), sample_id != "")

  if (length(exclude_ids) > 0) {
    cleaned <- cleaned %>% filter(!sample_id %in% exclude_ids)
  }
  if (!is.null(exclude_pattern)) {
    cleaned <- cleaned %>% filter(!str_detect(sample_id, exclude_pattern))
  }

  cleaned %>%
    distinct(mutation_key, chrom, pos, sample_id)
}

join_with_mapping <- function(mutation_df, mapping_df) {
  mapping_df <- mapping_df %>%
    distinct(original_id, .keep_all = TRUE)

  mutation_df %>%
    inner_join(mapping_df, by = c("sample_id" = "original_id"))
}

is_laminar_pair <- function(a, b) {
  overlap <- intersect(a, b)
  length(overlap) == 0 || length(overlap) == length(a) || length(overlap) == length(b)
}

collapse_set_support <- function(filtered_df, tip_labels) {
  n_tip <- length(tip_labels)

  filtered_df %>%
    distinct(mutation_key, new_id) %>%
    reframe(
      descendants = list(sort(unique(new_id))),
      size = dplyr::n_distinct(new_id),
      .by = mutation_key
    ) %>%
    filter(size >= 2, size < n_tip) %>%
    mutate(sample_signature = map_chr(descendants, ~ paste(.x, collapse = ";"))) %>%
    reframe(
      descendants = list(descendants[[1]]),
      size = first(size),
      support_n = dplyr::n(),
      mutation_keys = paste(sort(mutation_key), collapse = ";"),
      .by = sample_signature
    ) %>%
    arrange(desc(support_n), desc(size), sample_signature)
}

select_laminar_family <- function(candidate_sets, tip_labels) {
  root_desc <- sort(unique(tip_labels))
  kept_sets <- list(root_desc)
  keep_idx <- logical(nrow(candidate_sets))
  reject_reason <- rep(NA_character_, nrow(candidate_sets))

  for (i in seq_len(nrow(candidate_sets))) {
    current_set <- candidate_sets$descendants[[i]]
    ok <- all(vapply(kept_sets, function(existing) is_laminar_pair(current_set, existing), logical(1)))
    if (ok) {
      kept_sets[[length(kept_sets) + 1]] <- current_set
      keep_idx[i] <- TRUE
    } else {
      reject_reason[i] <- "non_laminar_overlap"
    }
  }

  list(
    selected = candidate_sets %>% filter(keep_idx),
    rejected = candidate_sets %>% mutate(reject_reason = reject_reason) %>% filter(!keep_idx)
  )
}

build_tree_structure <- function(selected_clades, tip_order) {
  root_tbl <- tibble(
    node_name = "Root",
    descendants = list(tip_order),
    node_size = length(tip_order),
    support_n = NA_integer_,
    node_type = "root"
  )

  clade_tbl <- selected_clades %>%
    transmute(
      node_name = paste0("Clade", row_number()),
      descendants = map(descendants, identity),
      node_size = size,
      support_n = support_n,
      sample_signature = sample_signature,
      mutation_keys = mutation_keys,
      node_type = "internal"
    )

  node_tbl <- bind_rows(root_tbl, clade_tbl) %>%
    mutate(node_idx = row_number())

  desc_list <- node_tbl$descendants
  names(desc_list) <- node_tbl$node_name

  is_subset <- function(x, y) {
    all(x %in% y)
  }

  get_internal_children <- function(idx) {
    current <- desc_list[[idx]]
    candidates <- which(
      seq_along(desc_list) != idx &
        vapply(desc_list, function(x) length(x) < length(current) && is_subset(x, current), logical(1))
    )

    if (length(candidates) == 0) {
      return(integer(0))
    }

    candidates[!vapply(
      candidates,
      function(candidate_idx) {
        any(vapply(
          setdiff(candidates, candidate_idx),
          function(other_idx) {
            length(desc_list[[candidate_idx]]) < length(desc_list[[other_idx]]) &&
              is_subset(desc_list[[candidate_idx]], desc_list[[other_idx]])
          },
          logical(1)
        ))
      },
      logical(1)
    )]
  }

  edge_rows <- list()
  child_order_key <- function(descendants) {
    min(match(descendants, tip_order))
  }

  for (idx in seq_len(nrow(node_tbl))) {
    internal_children <- get_internal_children(idx)
    covered_tips <- character(0)

    if (length(internal_children) > 0) {
      internal_children <- internal_children[order(vapply(internal_children, function(i) child_order_key(desc_list[[i]]), numeric(1)))]
      covered_tips <- sort(unique(unlist(desc_list[internal_children], use.names = FALSE)))
      edge_rows[[length(edge_rows) + 1]] <- tibble(
        parent_name = node_tbl$node_name[idx],
        child_name = node_tbl$node_name[internal_children],
        child_type = "internal"
      )
    }

    leaf_children <- setdiff(desc_list[[idx]], covered_tips)
    if (length(leaf_children) > 0) {
      leaf_children <- leaf_children[order(match(leaf_children, tip_order))]
      edge_rows[[length(edge_rows) + 1]] <- tibble(
        parent_name = node_tbl$node_name[idx],
        child_name = leaf_children,
        child_type = "tip"
      )
    }
  }

  edge_tbl <- bind_rows(edge_rows)
  internal_tbl <- node_tbl %>%
    filter(node_type != "tip") %>%
    mutate(phylo_id = length(tip_order) + row_number())

  tip_tbl <- tibble(
    child_name = tip_order,
    phylo_id = seq_along(tip_order)
  )

  edge_phylo <- edge_tbl %>%
    left_join(internal_tbl %>% select(parent_name = node_name, parent_id = phylo_id), by = "parent_name") %>%
    left_join(
      bind_rows(
        internal_tbl %>% transmute(child_name = node_name, child_id = phylo_id),
        tip_tbl %>% transmute(child_name, child_id = phylo_id)
      ),
      by = "child_name"
    ) %>%
    select(parent_id, child_id)

  tree <- list(
    edge = as.matrix(edge_phylo),
    tip.label = tip_order,
    Nnode = nrow(internal_tbl)
  )
  class(tree) <- "phylo"
  tree <- reorder.phylo(tree, order = "cladewise")

  node_export <- internal_tbl %>%
    select(node_name, phylo_id) %>%
    left_join(node_tbl %>% select(node_name, descendants, node_size, support_n, node_type), by = "node_name") %>%
    mutate(descendant_samples = map_chr(descendants, ~ paste(.x, collapse = ";")))

  list(tree = tree, node_table = node_export)
}

compute_tip_counts <- function(filtered_df) {
  mutation_sizes <- filtered_df %>%
    distinct(mutation_key, new_id) %>%
    count(mutation_key, name = "sample_n")

  filtered_df %>%
    distinct(mutation_key, new_id) %>%
    left_join(mutation_sizes, by = "mutation_key") %>%
    mutate(mutation_type = if_else(sample_n == 1, "unique", "shared")) %>%
    count(new_id, mutation_type, name = "mutation_n") %>%
    tidyr::pivot_wider(names_from = mutation_type, values_from = mutation_n, values_fill = 0L) %>%
    mutate(total = shared + unique)
}

make_rect_df <- function(tip_plot_df, x_start, x_scale) {
  long_df <- tip_plot_df %>%
    select(label, y, seed_group, shared, unique) %>%
    pivot_longer(cols = c(shared, unique), names_to = "mutation_type", values_to = "count") %>%
    mutate(
      mutation_type = factor(mutation_type, levels = c("shared", "unique")),
      cumulative_end = ave(count, label, FUN = cumsum),
      cumulative_start = cumulative_end - count,
      xmin = x_start + cumulative_start * x_scale,
      xmax = x_start + cumulative_end * x_scale,
      ymin = y - 0.32,
      ymax = y + 0.32
    )

  long_df %>% filter(count > 0)
}

get_descendant_tips <- function(tree) {
  n_tip <- Ntip(tree)
  n_node <- tree$Nnode
  children_map <- split(tree$edge[, 2], tree$edge[, 1])
  descendant_cache <- vector("list", n_tip + n_node)

  get_node_tips <- function(node_id) {
    cached <- descendant_cache[[node_id]]
    if (!is.null(cached)) {
      return(cached)
    }

    if (node_id <= n_tip) {
      descendant_cache[[node_id]] <<- node_id
      return(node_id)
    }

    child_ids <- children_map[[as.character(node_id)]]
    desc_ids <- sort(unique(unlist(lapply(child_ids, get_node_tips), use.names = FALSE)))
    descendant_cache[[node_id]] <<- desc_ids
    desc_ids
  }

  all_nodes <- seq_len(n_tip + n_node)
  descendant_ids <- lapply(all_nodes, get_node_tips)
  names(descendant_ids) <- as.character(all_nodes)
  descendant_ids
}

make_branch_color_df <- function(tree, rice_meta, seed_palette) {
  descendant_ids <- get_descendant_tips(tree)
  seed_lookup <- setNames(rice_meta$seed_group, rice_meta$new_id)

  tibble(
    node = seq_along(descendant_ids),
    descendant_tips = lapply(descendant_ids, function(idx) tree$tip.label[idx])
  ) %>%
    mutate(
      seed_groups = lapply(descendant_tips, function(x) unique(seed_lookup[x])),
      branch_seed_group = vapply(
        seed_groups,
        function(x) if (length(x) == 1) x else NA_character_,
        character(1)
      ),
      branch_color = ifelse(
        is.na(branch_seed_group),
        "#8d8d8d",
        unname(seed_palette[branch_seed_group])
      )
    ) %>%
    select(node, branch_seed_group, branch_color)
}

plot_rice_tree_numeric_nodes <- function(tree, rice_meta, counts_df, node_table, output_base) {
  seed_levels <- unique(rice_meta$seed_group)
  seed_palette <- c(
    "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
    "#e6ab02", "#a6761d", "#1f78b4", "#b2df8a", "#fb9a99",
    "#cab2d6", "#ff7f00", "#6a3d9a"
  )
  names(seed_palette) <- seed_levels

  tree_plot <- ggtree(tree, branch.length = "none", linewidth = 0.45, color = "#666666")
  branch_color_df <- make_branch_color_df(tree, rice_meta, seed_palette)
  tree_plot$data <- tree_plot$data %>%
    left_join(branch_color_df, by = "node") %>%
    mutate(branch_color = replace_na(branch_color, "#8d8d8d"))
  tip_df <- tree_plot$data %>%
    filter(isTip) %>%
    select(label, x, y) %>%
    left_join(rice_meta, by = c("label" = "new_id")) %>%
    left_join(counts_df, by = c("label" = "new_id")) %>%
    mutate(
      shared = replace_na(shared, 0L),
      unique = replace_na(unique, 0L),
      total = replace_na(total, 0L)
    )

  node_text_df <- tree_plot$data %>%
    filter(!isTip) %>%
    select(node, x, y) %>%
    left_join(
      node_table %>%
        transmute(node = phylo_id, shared_mutation_n = replace_na(as.numeric(support_n), 0)),
      by = "node"
    ) %>%
    mutate(shared_mutation_n = replace_na(shared_mutation_n, 0))

  root_node <- node_table %>%
    filter(node_name == "Root") %>%
    pull(phylo_id)
  seed_group_nodes <- tree$edge[tree$edge[, 1] == root_node, 2]
  seed_group_x <- min(node_text_df$x[node_text_df$node %in% seed_group_nodes], na.rm = TRUE)
  node_text_df <- node_text_df %>%
    mutate(
      display_x = if_else(node %in% seed_group_nodes, seed_group_x, x)
    ) %>%
    filter(shared_mutation_n > 0)

  max_tree_x <- max(tree_plot$data$x, na.rm = TRUE)
  max_total <- max(tip_df$total, na.rm = TRUE)
  bar_start <- max_tree_x + 0.7
  bar_span <- 4.0
  x_scale <- if (max_total > 0) bar_span / max_total else 0
  rect_df <- make_rect_df(tip_df, x_start = bar_start, x_scale = x_scale) %>%
    mutate(mutation_type = factor(mutation_type, levels = c("unique", "shared")))
  dot_x <- bar_start + bar_span + 0.8
  label_x <- dot_x + 0.45
  x_limit_max <- label_x + 1.8
  y_limit_min <- min(tip_df$y, na.rm = TRUE) - 1
  y_limit_max <- max(tip_df$y, na.rm = TRUE) + 3
  sample_label_df <- tip_df %>%
    arrange(y) %>%
    transmute(x = label_x, y = y, seed_group = seed_group)
  
  p <- tree_plot +
    geom_tree(
      data = tree_plot$data,
      aes(color = I(branch_color)),
      linewidth = 0.65,
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    geom_label(
      data = node_text_df,
      aes(x = display_x, y = y, label = shared_mutation_n),
      inherit.aes = FALSE,
      size = 2.8,
      fontface = "bold",
      linewidth = 0.22,
      label.padding = grid::unit(0.1, "lines"),
      fill = "white",
      color = "#2f2f2f"
    ) +
    geom_rect(
      data = rect_df,
      aes(
        xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
        fill = seed_group, alpha = mutation_type
      ),
      inherit.aes = FALSE,
      color = NA
    ) +
    geom_point(
      data = tip_df,
      aes(x = dot_x, y = y, fill = seed_group),
      inherit.aes = FALSE,
      shape = 21,
      size = 2.1,
      stroke = 0.2,
      color = "#333333"
    ) +
    geom_text(
      data = sample_label_df,
      aes(x = x, y = y, label = seed_group, color = seed_group),
      inherit.aes = FALSE,
      hjust = 0,
      size = 2.3,
      fontface = "bold"
    ) +
    geom_text(
      data = tibble(
        x = c(0.05, bar_start, label_x),
        y = max(tip_df$y, na.rm = TRUE) + 2,
        label = c("Node shared mutations", "Unique + shared mutations", "Seed group")
      ),
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      size = 3.3,
      fontface = "bold",
      color = "#303030"
      ) +
      scale_fill_manual(values = seed_palette) +
      scale_alpha_manual(
        values = c(unique = 1, shared = 0.3),
        breaks = c("unique", "shared"),
        labels = c("Unique", "Shared"),
        name = NULL
      ) +
      scale_color_manual(values = seed_palette) +
      coord_cartesian(
        xlim = c(0, x_limit_max),
        ylim = c(y_limit_min, y_limit_max),
      clip = "off"
    ) +
    labs(title = "Rice lineage topology from shared mutations") +
    theme_tree2() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "none",
      plot.margin = margin(12, 28, 12, 12),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    )

  tryCatch(
    ggsave(paste0(output_base, ".pdf"), p, width = 15, height = 11, units = "in", bg = "white"),
    error = function(e) message("Skipping PDF export for ", output_base, ": ", conditionMessage(e))
  )
  ggsave(paste0(output_base, ".png"), p, width = 15, height = 11, units = "in", dpi = 300, bg = "white")
}

plot_rice_tree <- function(tree, rice_meta, counts_df, node_table, output_base) {
  seed_levels <- unique(rice_meta$seed_group)
  seed_palette <- c(
    "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
    "#e6ab02", "#a6761d", "#1f78b4", "#b2df8a", "#fb9a99",
    "#cab2d6", "#ff7f00", "#6a3d9a"
  )
  names(seed_palette) <- seed_levels

  tree_plot <- ggtree(tree, branch.length = "none", linewidth = 0.45, color = "#666666")
  branch_color_df <- make_branch_color_df(tree, rice_meta, seed_palette)
  tree_plot$data <- tree_plot$data %>%
    left_join(branch_color_df, by = "node") %>%
    mutate(branch_color = replace_na(branch_color, "#8d8d8d"))
  tip_df <- tree_plot$data %>%
    filter(isTip) %>%
    select(label, x, y) %>%
    left_join(rice_meta, by = c("label" = "new_id")) %>%
    left_join(counts_df, by = c("label" = "new_id")) %>%
    mutate(
      shared = replace_na(shared, 0L),
      unique = replace_na(unique, 0L),
      total = replace_na(total, 0L)
    )

  node_fill_df <- tree_plot$data %>%
    filter(!isTip) %>%
    select(node, x, y) %>%
    left_join(
      node_table %>%
        transmute(node = phylo_id, shared_mutation_n = replace_na(as.numeric(support_n), 0)),
      by = "node"
    ) %>%
    mutate(shared_mutation_n = replace_na(shared_mutation_n, 0))

  max_tree_x <- max(tree_plot$data$x, na.rm = TRUE)
  max_total <- max(tip_df$total, na.rm = TRUE)
  bar_start <- max_tree_x + 0.7
  bar_span <- 4.0
  x_scale <- if (max_total > 0) bar_span / max_total else 0
  rect_df <- make_rect_df(tip_df, x_start = bar_start, x_scale = x_scale)
  dot_x <- bar_start + bar_span + 0.8
  label_x <- dot_x + 0.45
  sample_label_df <- tip_df %>%
    arrange(y) %>%
    transmute(x = label_x, y = y, seed_group = seed_group)
  
  p <- tree_plot +
    geom_tree(
      data = tree_plot$data,
      aes(color = I(branch_color)),
      linewidth = 0.65,
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    geom_point(
      data = node_fill_df,
      aes(x = x, y = y, fill = shared_mutation_n),
      inherit.aes = FALSE,
      shape = 21,
      size = 3.0,
      stroke = 0.28,
      color = "#4a4a4a"
    ) +
    scale_fill_gradientn(
      colours = c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b")
    ) +
    ggnewscale::new_scale_fill() +
    geom_rect(
      data = rect_df,
      aes(
        xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
        fill = seed_group, alpha = mutation_type
      ),
      inherit.aes = FALSE,
      color = NA
    ) +
    geom_point(
      data = tip_df,
      aes(x = dot_x, y = y, fill = seed_group),
      inherit.aes = FALSE,
      shape = 21,
      size = 2.1,
      stroke = 0.2,
      color = "#333333"
    ) +
    geom_text(
      data = sample_label_df,
      aes(x = x, y = y, label = seed_group, color = seed_group),
      inherit.aes = FALSE,
      hjust = 0,
      size = 2.3,
      fontface = "bold"
    ) +
    geom_text(
      data = tibble(
        x = c(0.05, bar_start, label_x),
        y = max(tip_df$y, na.rm = TRUE) + 2,
        label = c("Node shared mutations", "Shared + unique mutations", "Seed group")
      ),
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      size = 3.3,
      fontface = "bold",
      color = "#303030"
    ) +
    scale_fill_manual(values = seed_palette) +
      scale_alpha_manual(
        values = c(shared = 1, unique = 0.35),
        breaks = c("shared", "unique"),
        labels = c("Shared", "Unique"),
        name = NULL
      ) +
      scale_color_manual(values = seed_palette) +
      coord_cartesian(
        xlim = c(0, label_x + 1.8),
        ylim = c(min(tip_df$y, na.rm = TRUE) - 1, max(tip_df$y, na.rm = TRUE) + 3),
      clip = "off"
    ) +
    labs(title = "Rice lineage topology from shared mutations") +
    theme_tree2() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "none",
      plot.margin = margin(12, 28, 12, 12),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    )

  tryCatch(
    ggsave(paste0(output_base, ".pdf"), p, width = 15, height = 11, units = "in", bg = "white"),
    error = function(e) message("Skipping PDF export for ", output_base, ": ", conditionMessage(e))
  )
  ggsave(paste0(output_base, ".png"), p, width = 15, height = 11, units = "in", dpi = 300, bg = "white")
}

write_clean_exports <- function() {
  export_list <- list()

  for (species_name in names(species_cfg)) {
    cfg <- species_cfg[[species_name]]
    mapping_df <- read_mapping_sheet(cfg$sheet)
    mutation_df <- read_mutation_sheet(cfg$sheet, cfg$exclude_ids, cfg$exclude_pattern)
    filtered_df <- join_with_mapping(mutation_df, mapping_df)
    export_list[[species_name]] <- filtered_df
  }

  write.xlsx(
    export_list,
    file = file.path(output_dir, "cleaned_mutation_tables.xlsx"),
    overwrite = TRUE
  )
}

rice_mapping <- read_mapping_sheet("rice") %>%
  mutate(
    seed_group = str_replace(new_id, "-[^-]+$", ""),
    sample_rank = as.numeric(str_extract(new_id, "\\d+$"))
  ) %>%
  arrange(factor(seed_group, levels = unique(seed_group)), sample_rank, new_id)

rice_mutation <- read_mutation_sheet("rice")
rice_filtered <- join_with_mapping(rice_mutation, rice_mapping)

write_clean_exports()

tip_order <- rice_mapping$new_id
candidate_sets <- collapse_set_support(rice_filtered, tip_order)
seed_lookup <- setNames(rice_mapping$seed_group, rice_mapping$new_id)
candidate_sets <- candidate_sets %>%
  mutate(
    seed_groups = map(descendants, ~ unique(seed_lookup[.x])),
    is_within_seed = lengths(seed_groups) == 1
  ) %>%
  filter(is_within_seed) %>%
  select(-seed_groups, -is_within_seed)

seed_group_sets <- rice_mapping %>%
  summarise(
    descendants = list(new_id),
    size = dplyr::n(),
    support_n = 0L,
    mutation_keys = NA_character_,
    .by = seed_group
  ) %>%
  mutate(sample_signature = map_chr(descendants, ~ paste(.x, collapse = ";"))) %>%
  select(sample_signature, descendants, size, support_n, mutation_keys)

candidate_sets <- bind_rows(candidate_sets, seed_group_sets) %>%
  arrange(desc(support_n), desc(size), sample_signature) %>%
  distinct(sample_signature, .keep_all = TRUE)

laminar_sets <- select_laminar_family(candidate_sets, tip_order)
tree_parts <- build_tree_structure(laminar_sets$selected, tip_order)
tip_counts <- compute_tip_counts(rice_filtered)

plot_rice_tree(
  tree = rotateConstr(tree_parts$tree, tip_order),
  rice_meta = rice_mapping,
  counts_df = tip_counts,
  node_table = tree_parts$node_table,
  output_base = file.path(output_dir, "rice_lineage_topology")
)

plot_rice_tree_numeric_nodes(
  tree = rotateConstr(tree_parts$tree, tip_order),
  rice_meta = rice_mapping,
  counts_df = tip_counts,
  node_table = tree_parts$node_table,
  output_base = file.path(output_dir, "rice_lineage_topology_numeric_nodes")
)

selected_export <- laminar_sets$selected %>%
  mutate(descendant_samples = map_chr(descendants, ~ paste(.x, collapse = ";"))) %>%
  select(-descendants)

rejected_export <- laminar_sets$rejected %>%
  mutate(descendant_samples = map_chr(descendants, ~ paste(.x, collapse = ";"))) %>%
  select(-descendants)

write.xlsx(
  list(
    rice_mapping = rice_mapping,
    rice_mutations = rice_filtered,
    rice_tip_mutation_counts = tip_counts,
    rice_selected_clades = selected_export,
    rice_rejected_clades = rejected_export,
    rice_tree_nodes = tree_parts$node_table
  ),
  file = file.path(output_dir, "rice_lineage_support_tables.xlsx"),
  overwrite = TRUE
)

message("Rice lineage analysis finished.")
