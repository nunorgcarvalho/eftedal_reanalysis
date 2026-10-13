# Figures for the main analysis (scripts/01-05). Run from anywhere:
#   Rscript R/make_figures.R
# Writes PNGs to results/figures/main/. Each block is skipped if its inputs are missing.
args <- commandArgs(trailingOnly = FALSE)
this <- sub("^--file=", "", args[grep("^--file=", args)])
REPO <- normalizePath(file.path(dirname(this), ".."))
RESULTS <- file.path(REPO, "results")
FIGURES <- file.path(RESULTS, "figures", "main")
source(file.path(REPO, "R", "plot_functions.R"))

rd <- function(...) as_tibble(fread(file.path(RESULTS, ...)))
has <- function(...) file.exists(file.path(RESULTS, ...))

WEIGHT_LABELS <- c(hybrid = "hybrid (type-balanced)", invvar = "inverse variance (1/SE²)",
                   hybrid_block = "hybrid, balanced over blocks", invSqErr = "Eftedal (log, 1/SE²)")
COHORT_LABELS <- c(bio_only = "Bio only", ALL_but_adopted_UZ = "Bio + in-law")
MODEL_ORDER <- c("A", "A+mat", "A+H", "A+H+mat", "A+F", "A+F+mat", "A+H+F", "A+H+F+mat",
                 "A+GxG", "A+GxG+mat", "A+GxG+H", "A+GxG+H+mat", "A+GxG+F", "A+GxG+F+mat",
                 "A+GxG+H+F", "A+GxG+H+F+mat")
CHOSEN <- c("A", "A+GxG", "A+F", "A+GxG+F", "A+F+mat")
FOLD_ORDER <- c("twins", "sibs", "half sibs", "1st cousins", "half cousins", "2nd cousins",
                "step sibs", "cocousins", "adoptive sibs", "adoptive cousins",
                "all in-laws", "all adoptive")
obs_from <- function(p) p %>% distinct(row, .keep_all = TRUE) %>% select(row, rel, side, sex, corr, se)
fmt_par <- function(d, ps) {
  d %>% select(any_of(ps)) %>% pivot_longer(everything()) %>% filter(!is.na(value)) %>%
    mutate(s = sprintf("%s=%.2f", sub("V_", "", name), value)) %>% pull(s) %>% paste(collapse = ", ")
}

# ---- 1. Eftedal et al.'s fit vs ours, in their bio + in-law cohort ----------------------
if (has("robustness", "contrast_predictions.csv")) {
  cs <- rd("robustness", "contrast.csv")
  cp <- rd("robustness", "contrast_predictions.csv")
  obs <- obs_from(cp)
  pick <- function(setting, w) cs %>% filter(setting == !!setting, weighting == !!w,
                                             cohort_name == "ALL_but_adopted_UZ") %>% pull(fit_id)
  ids <- c(pick("Eftedal settings", "invSqErr"), pick("main", "invvar"), pick("main", "hybrid"))
  labs <- sprintf(c("Eftedal et al. (log, 1/SE²): h²=%.2f", "ours, inverse variance: h²=%.2f",
                    "ours, hybrid: h²=%.2f"), cs$V_A[match(ids, cs$fit_id)])
  pr <- bind_rows(lapply(seq_along(ids), function(i) cp %>% filter(fit_id == ids[i]) %>% mutate(model = labs[i])))
  save_fig(plot_fits(obs, pr, colors = setNames(c("#F8766D", "#7CAE00", "#00BFC4"), labs)),
           "01_eftedal_vs_ours", height = 1850)

  # ---- 2. headline: bio-only vs bio + in-law -----------------------------------------
  d <- cs %>% mutate(cohort = factor(COHORT_LABELS[cohort_name], levels = COHORT_LABELS),
                     setting = factor(paste(setting, "|", weighting),
                                      levels = rev(unique(paste(setting, "|", weighting)))))
  g <- ggplot(d, aes(x = V_A, y = setting)) +
    geom_line(aes(group = setting), color = "grey60", linewidth = 1) +
    geom_point(aes(color = cohort), size = 3.5) +
    scale_color_manual(values = c("Bio only" = "#00BFC4", "Bio + in-law" = "#F8766D"), name = NULL) +
    labs(x = "h² (additive model)", y = NULL,
         title = "Adding in-laws: Eftedal et al.'s fitting vs. ours") +
    theme_bw() + theme(legend.position = "bottom")
  save_fig(g, "02_bio_vs_inlaw_contrast", width = 2400, height = 1200)
}
if (has("robustness", "phenotypes.csv")) {
  ph <- rd("robustness", "phenotypes.csv") %>%
    select(phenotype, setting, weighting, cohort_name, V_A) %>%
    pivot_wider(names_from = cohort_name, values_from = V_A) %>%
    mutate(gap = bio_only - ALL_but_adopted_UZ, setting = paste(setting, "|", weighting))
  g <- ggplot(ph, aes(x = gap, y = reorder(phenotype, gap), color = setting)) +
    geom_vline(xintercept = 0) + geom_point(size = 2.5) +
    labs(x = "h²(bio only) − h²(bio + in-law)", y = NULL, color = NULL) +
    theme_bw() + theme(legend.position = "bottom", legend.direction = "vertical")
  save_fig(g, "02b_contrast_all_phenotypes", width = 2200, height = 2000)
}

# ---- 3. block cross-validation ------------------------------------------------------------
if (has("cv", "cv_paired.csv")) {
  pc <- rd("cv", "cv_paired.csv") %>%
    left_join(rd("cv", "cv_summary.csv") %>% select(model, weights, cv_rmse), by = c("model", "weights")) %>%
    mutate(model = factor(model, levels = rev(MODEL_ORDER)), weights = WEIGHT_LABELS[weights])
  g <- ggplot(pc, aes(x = mean_diff, y = model)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_errorbar(aes(xmin = mean_diff - 1.96 * se_diff, xmax = mean_diff + 1.96 * se_diff),
                  width = 0.3, orientation = "y") +
    geom_point(size = 2.5) + facet_wrap(~weights, scales = "free_x") +
    labs(x = "Held-out score minus best model's (paired over 10 blocks; ±1.96 SE)", y = NULL,
         title = "Block cross-validation: which components help predict unseen relationship types") +
    theme_bw()
  save_fig(g, "03_cv_model_comparison", width = 3000, height = 1700)

  fsc <- rd("cv", "cv_fold_scores.csv") %>%
    mutate(model = factor(model, levels = rev(MODEL_ORDER)), fold = factor(fold, levels = FOLD_ORDER),
           weights = WEIGHT_LABELS[weights])
  g <- ggplot(fsc, aes(x = fold, y = model, fill = rmse)) + geom_tile() +
    geom_text(aes(label = sprintf("%.3f", rmse)), size = 2.2) +
    facet_wrap(~weights) + scale_fill_viridis_c(option = "magma", direction = -1, name = "RMSE") +
    labs(x = "Held-out block", y = NULL) + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_fig(g, "03b_cv_by_block", width = 3600, height = 1800)
}

# ---- 4. fits of the chosen models ---------------------------------------------------------
if (has("main", "fits_predictions.csv")) {
  fs <- rd("main", "fits_summary.csv")
  fp <- rd("main", "fits_predictions.csv")
  obs <- obs_from(fp)
  for (wt in c("hybrid", "invvar")) {
    for (mod in c("A+F", "A+GxG+F", "A+GxG")) {
      ia <- fs %>% filter(model == "A", weights == !!wt) %>% pull(fit_id)
      ib <- fs %>% filter(model == !!mod, weights == !!wt) %>% pull(fit_id)
      labs <- c(paste0("A: ", fmt_par(fs %>% filter(fit_id == ia), c("V_A", "m", "r"))),
                paste0(mod, ": ", fmt_par(fs %>% filter(fit_id == ib),
                                        c("V_A", "V_GxG", "V_family", "r_f", "m", "r"))))
      pr <- bind_rows(fp %>% filter(fit_id == ia) %>% mutate(model = labs[1]),
                      fp %>% filter(fit_id == ib) %>% mutate(model = labs[2]))
      save_fig(plot_fits(obs, pr, colors = setNames(c("#00BFC4", "#C77CFF"), labs)), height = 1800,
               sprintf("04_fit_A_vs_%s_%s", gsub("\\+", "", mod), wt))
    }
  }
  # variance components of every core model
  comp <- fs %>% select(model, weights, V_A, V_GxG, V_house, V_family) %>%
    pivot_longer(starts_with("V_"), names_to = "component") %>% filter(!is.na(value)) %>%
    mutate(component = factor(sub("V_", "", component), c("A", "GxG", "house", "family")),
           model = factor(model, levels = MODEL_ORDER), weights = WEIGHT_LABELS[weights])
  g <- ggplot(comp, aes(x = model, y = value, fill = component)) + geom_col() +
    facet_wrap(~weights, ncol = 1) + scale_fill_brewer(palette = "Set2") +
    labs(x = NULL, y = "Variance component", fill = NULL) + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_fig(g, "04b_components_all_models", width = 3000, height = 2000)
}

# ---- 5. parameter intervals ---------------------------------------------------------------
if (has("uncertainty_main", "intervals.csv")) {
  iv <- rd("uncertainty_main", "intervals.csv") %>%
    filter(param %in% c("V_A", "V_GxG", "V_family", "r_f", "delta", "m", "r")) %>%
    mutate(model = factor(model, levels = rev(CHOSEN)), weights = WEIGHT_LABELS[weights],
           method = factor(method, c("param_boot x1", "param_boot x1.3", "sandwich", "block jackknife")),
           lb = ifelse(param == "delta", -0.5, 0), lo = pmax(lo, lb), hi = pmin(hi, 1))
  g <- ggplot(iv, aes(x = estimate, xmin = lo, xmax = hi, y = model, color = method)) +
    geom_pointrange(position = position_dodge(width = 0.6), size = 0.25) +
    facet_grid(weights ~ param, scales = "free_x") +
    labs(x = "Estimate and 95% interval", y = NULL, color = NULL) +
    theme_bw() + theme(legend.position = "bottom")
  save_fig(g, "05_parameter_intervals", width = 3800, height = 1900)

  jk <- rd("uncertainty_main", "jackknife.csv") %>% filter(fold %in% FOLD_ORDER[1:10]) %>%
    select(model, weights, fold, V_A) %>%
    mutate(model = factor(model, levels = CHOSEN), fold = factor(fold, levels = rev(FOLD_ORDER)),
           weights = WEIGHT_LABELS[weights])
  g <- ggplot(jk, aes(x = V_A, y = fold, color = model)) + geom_point(size = 2.3) +
    facet_wrap(~weights) + labs(x = "V_A with the block left out", y = "Block left out", color = NULL) +
    theme_bw()
  save_fig(g, "05b_jackknife_VA", width = 2800, height = 1500)
}

# ---- 6. step-sibling kernel ---------------------------------------------------------------
if (has("robustness", "additive_predictions.csv")) {
  ar <- rd("robustness", "additive_rows.csv")
  ap <- rd("robustness", "additive_predictions.csv")
  ids <- ar %>% filter(cohort_name == "ALL", variant %in% c("main", "Eftedal step-sib kernel")) %>%
    select(fit_id, variant, weighting, V_A)
  st <- ap %>% filter(rel == "step sibs") %>% group_by(fit_id) %>% summarize(pred = first(pred)) %>%
    inner_join(ids, by = "fit_id") %>%
    mutate(kernel = ifelse(variant == "main", "m((1+r)/2)² [corrected]", "m((1+m)/2)² [Eftedal]"),
           weighting = WEIGHT_LABELS[weighting])
  so <- ap %>% filter(rel == "step sibs") %>% distinct(row, .keep_all = TRUE) %>%
    summarize(corr = weighted.mean(corr, 1 / se^2), se = sqrt(1 / sum(1 / se^2)))
  g <- ggplot(st, aes(x = pred, y = weighting, color = kernel)) +
    annotate("rect", xmin = so$corr - z95 * so$se, xmax = so$corr + z95 * so$se, ymin = -Inf, ymax = Inf, alpha = 0.15) +
    geom_vline(xintercept = so$corr, linetype = "dashed") +
    geom_point(size = 3.5) +
    geom_text(aes(label = sprintf("h²=%.3f", V_A), vjust = ifelse(grepl("Eftedal", kernel), -1.1, 2.1)),
              size = 3, show.legend = FALSE) +
    labs(x = "Predicted step-sibling correlation (additive model; shaded: observed 95% CI)", y = NULL, color = NULL) +
    theme_bw() + theme(legend.position = "bottom")
  save_fig(g, "06_stepsib_kernel", width = 2400, height = 1100)

  # ---- 7. robustness of h2 ----------------------------------------------------------------
  d <- ar %>% filter(!is.na(V_A)) %>%
    mutate(label = paste0(variant, " [", weighting, "]"),
           cohort_name = ifelse(cohort_name == "ALL", "All 82 categories", "Bio + in-law (Eftedal cohort)"))
  g <- ggplot(d, aes(x = V_A, y = reorder(label, V_A))) + geom_point(size = 2.5) +
    facet_wrap(~cohort_name) + labs(x = "h² (additive model)", y = NULL) + theme_bw()
  save_fig(g, "07_robustness_h2", width = 3000, height = 1500)
}
cat("figures written to", FIGURES, "\n")
