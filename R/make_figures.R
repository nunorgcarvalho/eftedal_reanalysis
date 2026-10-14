# Figures for the analysis (scripts/01-06). Run from anywhere:
#   Rscript R/make_figures.R
# Writes PNG + PDF to results/figures/: main figures m*, supplementary figures s*. Each block is
# skipped if its inputs are missing. All intervals are block-jackknife 95% intervals unless
# labelled otherwise.
args <- commandArgs(trailingOnly = FALSE)
this <- sub("^--file=", "", args[grep("^--file=", args)])
REPO <- normalizePath(file.path(dirname(this), ".."))
RESULTS <- file.path(REPO, "results")
FIGURES <- file.path(RESULTS, "figures")
source(file.path(REPO, "R", "plot_functions.R"))

rd <- function(...) as_tibble(fread(file.path(RESULTS, ...)))
has <- function(...) file.exists(file.path(RESULTS, ...))

WEIGHT_LABELS <- c(invvar = "inverse variance (1/SE²)", hybrid_block = "hybrid (block-balanced)",
                   invSqErr = "Eftedal (log, 1/SE²)")
CORE_ORDER <- c("A", "A+H", "A+F", "A>GxG", "A>GxG+H", "A>GxG+F")
JOINT_ORDER <- c("A", "A+H", "A+F", "A+H+F", "A+GxG", "A+GxG+H", "A+GxG+F", "A+GxG+H+F")
show_model <- function(x) gsub(">", " → ", x)
COMP_ORDER <- c("A", "GxG", "house", "family")
COMP_FILL <- c(A = "#66C2A5", GxG = "#FC8D62", house = "#8DA0CB", family = "#E5C494")
COMP_LINE <- c(A = "#1B7A5A", GxG = "#B0441A", house = "#3F5694", family = "#9C7A2E")
FOLD_ORDER <- c("twins", "sibs", "half sibs", "1st cousins", "half cousins", "2nd cousins",
                "step sibs", "cocousins", "adoptive sibs", "adoptive cousins",
                "all in-laws", "all adoptive")
obs_from <- function(p) p %>% distinct(row, .keep_all = TRUE) %>% select(row, rel, side, sex, corr, se)
fmt_par <- function(d, ps) {
  d %>% select(any_of(ps)) %>% pivot_longer(everything()) %>% filter(!is.na(value)) %>%
    mutate(s = sprintf("%s=%.2f", sub("V_", "", name), value)) %>% pull(s) %>% paste(collapse = ", ")
}

# Stacked variance components with a jackknife interval on each component (drawn at the top of
# its segment) and on the sum (black, right of the bar).
stacked_components <- function(jk, order, title) {
  d <- jk %>% filter(param %in% paste0("V_", COMP_ORDER)) %>%
    mutate(component = factor(sub("V_", "", param), COMP_ORDER),
           model = factor(model, order), x = as.numeric(model)) %>%
    filter(!is.na(estimate)) %>%
    arrange(weights, model, component) %>%
    group_by(weights, model) %>%
    mutate(top = cumsum(estimate), k = as.numeric(component),
           xe = x - 0.28 + 0.12 * (k - 1),
           ymin = top - (estimate - pmax(lo, 0)), ymax = top + (pmin(hi, 1) - estimate)) %>%
    ungroup()
  s <- jk %>% filter(param == "V_sum") %>% mutate(model = factor(model, order), x = as.numeric(model))
  ggplot() +
    geom_col(data = d, aes(x = x, y = estimate, fill = component), width = 0.7,
             position = position_stack(reverse = TRUE)) +
    geom_errorbar(data = d, aes(x = xe, ymin = ymin, ymax = ymax, color = component),
                  width = 0.08, linewidth = 0.7, show.legend = FALSE) +
    geom_errorbar(data = s, aes(x = x + 0.42, ymin = pmax(lo, 0), ymax = pmin(hi, 1.1)),
                  width = 0.1, linewidth = 0.9) +
    geom_point(data = s, aes(x = x + 0.42, y = estimate), shape = 18, size = 3) +
    facet_wrap(~weights, ncol = 1, labeller = as_labeller(WEIGHT_LABELS)) +
    scale_x_continuous(breaks = seq_along(order), labels = show_model(order)) +
    scale_fill_manual(values = COMP_FILL, name = NULL) +
    scale_color_manual(values = COMP_LINE) +
    coord_cartesian(ylim = c(0, 1.1)) +
    labs(x = NULL, y = "Variance component (proportion of phenotypic variance)", title = title,
         subtitle = "Coloured bars: 95% block-jackknife interval of each component; black diamond and bar: sum of the components") +
    theme_bw() + theme(legend.position = "bottom", panel.grid.major.x = element_blank())
}

# ---- M1. Eftedal's estimator -> ours ----------------------------------------------------------
if (has("corrections", "steps.csv")) {
  st <- rd("corrections", "steps.csv") %>%
    mutate(label = factor(label, rev(unique(label[order(order)]))),
           kind = ifelse(aside, "aside (not cumulative)", "cumulative step"),
           param = factor(param, c("V_A", "m", "r"), c("h² (V_A)", "m (mate genetic corr.)", "r (mate phenotypic corr.)")))
  g <- ggplot(st, aes(x = estimate, xmin = pmax(lo, 0), xmax = pmin(hi, 1), y = label, color = kind)) +
    geom_pointrange(size = 0.45, linewidth = 0.8) +
    facet_wrap(~param, scales = "free_x") +
    scale_color_manual(values = c("cumulative step" = "black", "aside (not cumulative)" = "grey55"), name = NULL) +
    labs(x = "Estimate (95% block-jackknife interval)", y = NULL,
         title = "Additive model: from Eftedal et al.'s estimator to ours, one correction at a time") +
    theme_bw() + theme(legend.position = "bottom")
  save_fig(g, "m1_corrections", width = 3300, height = 1350)
}
if (has("robustness", "additive_predictions.csv") && has("models", "core_predictions.csv")) {
  ar <- rd("robustness", "additive_rows.csv")
  ap <- rd("robustness", "additive_predictions.csv")
  cp <- rd("models", "core_predictions.csv") %>% filter(model == "A")
  cf <- rd("models", "core_fits.csv") %>% filter(model == "A")
  ie <- ar %>% filter(variant == "Eftedal settings", cohort_name == "ALL_but_adopted_UZ") %>% pull(fit_id)
  labs <- c(sprintf("Eftedal et al.: h²=%.2f", ar$V_A[ar$fit_id == ie]),
            sprintf("ours, inverse variance: h²=%.2f", cf$V_A[cf$weights == "invvar"]),
            sprintf("ours, hybrid (block-balanced): h²=%.2f", cf$V_A[cf$weights == "hybrid_block"]))
  pr <- bind_rows(ap %>% filter(fit_id == ie) %>% mutate(model = labs[1]),
                  cp %>% filter(weights == "invvar") %>% mutate(model = labs[2], in_cohort = TRUE),
                  cp %>% filter(weights == "hybrid_block") %>% mutate(model = labs[3], in_cohort = TRUE))
  save_fig(plot_fits(obs_from(cp), pr, colors = setNames(c("#F8766D", "#7CAE00", "#00BFC4"), labs)),
           "m1b_fits_eftedal_vs_ours", height = 1850)
}

# ---- M2. in-laws ----------------------------------------------------------------------------
if (has("corrections", "inlaw.csv")) {
  il <- rd("corrections", "inlaw.csv") %>%
    mutate(label = factor(label, rev(unique(label))),
           panel = ifelse(quantity == "difference", "Difference (bio only − with in-laws)", "h² (additive model)"),
           quantity = factor(quantity, c("bio only", "bio + in-law", "difference")))
  g <- ggplot(il, aes(x = estimate, xmin = lo, xmax = hi, y = label, color = quantity)) +
    geom_vline(data = tibble(panel = "Difference (bio only − with in-laws)", v = 0), aes(xintercept = v),
               linetype = "dashed") +
    geom_pointrange(position = position_dodge(width = 0.45), size = 0.45, linewidth = 0.8) +
    facet_wrap(~panel, scales = "free_x") +
    scale_color_manual(values = c("bio only" = "#00BFC4", "bio + in-law" = "#F8766D", difference = "black"), name = NULL) +
    labs(x = "Estimate (95% block-jackknife interval; difference: paired jackknife)", y = NULL,
         title = "Does adding in-laws change h²?") +
    theme_bw() + theme(legend.position = "bottom")
  save_fig(g, "m2_inlaw", width = 3000, height = 1100)
}

# ---- M3. core models ----------------------------------------------------------------------
if (has("models", "core_jackknife.csv")) {
  save_fig(stacked_components(rd("models", "core_jackknife.csv"), CORE_ORDER,
                              "Core models: variance components"),
           "m3_components_core", width = 2700, height = 2300)
}
if (has("models", "core_predictions.csv")) {
  fs <- rd("models", "core_fits.csv")
  fp <- rd("models", "core_predictions.csv")
  obs <- obs_from(fp)
  for (wt in c("hybrid_block", "invvar")) {
    for (mod in c("A+F", "A>GxG+F", "A+H")) {
      fa <- fs %>% filter(model == "A", weights == !!wt)
      fb <- fs %>% filter(model == !!mod, weights == !!wt)
      labs <- c(paste0("A: ", fmt_par(fa, c("V_A", "m", "r"))),
                paste0(show_model(mod), ": ", fmt_par(fb, c("V_A", "V_GxG", "V_house", "pi", "V_family", "r_f", "m", "r"))))
      pr <- bind_rows(fp %>% filter(model == "A", weights == !!wt) %>% mutate(model = labs[1]),
                      fp %>% filter(model == !!mod, weights == !!wt) %>% mutate(model = labs[2]))
      save_fig(plot_fits(obs, pr, colors = setNames(c("#00BFC4", "#C77CFF"), labs)) +
                 labs(title = paste("Weighting:", WEIGHT_LABELS[[wt]])), height = 1850,
               sprintf("m3b_fit_A_vs_%s_%s", gsub("[+>]", "", mod), wt))
    }
  }
}

# ---- S1. cross-validation -------------------------------------------------------------------
for (set in c("core", "joint")) {
  if (!has("models", paste0(set, "_cv_paired.csv"))) next
  order <- if (set == "core") CORE_ORDER else JOINT_ORDER
  pc <- rd("models", paste0(set, "_cv_paired.csv")) %>%
    left_join(rd("models", paste0(set, "_cv_summary.csv")) %>% select(model, weights, cv_rmse), by = c("model", "weights")) %>%
    mutate(is_best = model == best, model = factor(show_model(model), rev(show_model(order))),
           weights = WEIGHT_LABELS[weights])
  g <- ggplot(pc, aes(x = mean_diff, y = model)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_errorbar(aes(xmin = mean_diff - 1.96 * se_diff, xmax = mean_diff + 1.96 * se_diff),
                  width = 0.3, orientation = "y") +
    geom_point(aes(color = is_best, size = is_best)) +
    geom_text(data = pc %>% filter(is_best), aes(label = sprintf("best: CV RMSE = %.4f", cv_rmse)),
              color = "#D55E00", hjust = -0.08, vjust = -0.9, size = 3.2) +
    scale_color_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "black"), guide = "none") +
    scale_size_manual(values = c(`TRUE` = 4, `FALSE` = 2.5), guide = "none") +
    facet_wrap(~weights, scales = "free_x") +
    labs(x = "CV score minus the best model's (squared-correlation units; paired over the blocks, ±1.96 SE)",
         y = NULL, title = sprintf("Block cross-validation (%s models)", set)) +
    theme_bw()
  save_fig(g, sprintf("s1_cv_%s", set), width = 3000, height = 300 + 170 * length(order))
  fsc <- rd("models", paste0(set, "_cv_scores.csv")) %>%
    mutate(model = factor(show_model(model), rev(show_model(order))), fold = factor(fold, FOLD_ORDER),
           weights = WEIGHT_LABELS[weights])
  g <- ggplot(fsc, aes(x = fold, y = model, fill = rmse)) + geom_tile() +
    geom_text(aes(label = sprintf("%.3f", rmse)), size = 2.2) +
    facet_wrap(~weights) + scale_fill_viridis_c(option = "magma", direction = -1, name = "RMSE") +
    labs(x = "Held-out block", y = NULL, title = sprintf("Held-out error by block (%s models)", set)) +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_fig(g, sprintf("s1b_cv_by_block_%s", set), width = 3600, height = 400 + 160 * length(order))
}

# ---- S2. joint model set; two-step vs joint ---------------------------------------------------
if (has("models", "joint_jackknife.csv")) {
  jj <- rd("models", "joint_jackknife.csv")
  save_fig(stacked_components(jj, JOINT_ORDER, "All-joint models: variance components"),
           "s2_components_joint", width = 3000, height = 2300)
  if (has("models", "core_jackknife.csv")) {
    pairs <- c("A>GxG" = "A+GxG", "A>GxG+H" = "A+GxG+H", "A>GxG+F" = "A+GxG+F")
    tj <- bind_rows(rd("models", "core_jackknife.csv") %>% filter(model %in% names(pairs)) %>%
                      mutate(procedure = "two-step (additive part frozen)", pair = pairs[model]),
                    jj %>% filter(model %in% pairs) %>% mutate(procedure = "joint", pair = model)) %>%
      filter(param %in% c("V_A", "V_GxG", "V_house", "V_family", "V_sum")) %>%
      mutate(pair = factor(sub("A\\+GxG", "A + GxG", pair), rev(sub("A\\+GxG", "A + GxG", pairs))),
             weights = WEIGHT_LABELS[weights])
    g <- ggplot(tj, aes(x = estimate, xmin = pmax(lo, 0), xmax = pmin(hi, 1.1), y = pair, color = procedure)) +
      geom_pointrange(position = position_dodge(width = 0.5), size = 0.3) +
      facet_grid(weights ~ param) +
      scale_color_manual(values = c("joint" = "#00BFC4", "two-step (additive part frozen)" = "#F8766D"), name = NULL) +
      labs(x = "Estimate (95% block-jackknife interval)", y = NULL, title = "Two-step vs joint estimation of GxG") +
      theme_bw() + theme(legend.position = "bottom")
    save_fig(g, "s2b_twostep_vs_joint", width = 3300, height = 1500)
  }
}

# ---- S3. bootstrap vs jackknife; which block moves the estimates -----------------------------
if (has("bootstrap", "intervals.csv") && has("models", "core_jackknife.csv")) {
  iv <- bind_rows(rd("bootstrap", "intervals.csv") %>% mutate(method = "parametric bootstrap"),
                  rd("models", "core_jackknife.csv") %>% mutate(method = "block jackknife")) %>%
    filter(param %in% c("V_A", "V_GxG", "V_house", "V_family", "V_sum", "pi", "r_f", "r")) %>%
    mutate(model = factor(show_model(model), rev(show_model(CORE_ORDER))), weights = WEIGHT_LABELS[weights],
           lo = pmax(lo, 0), hi = pmin(hi, 1.1))
  g <- ggplot(iv, aes(x = estimate, xmin = lo, xmax = hi, y = model, color = method)) +
    geom_pointrange(position = position_dodge(width = 0.6), size = 0.25) +
    facet_grid(weights ~ param, scales = "free_x") +
    labs(x = "Estimate and 95% interval", y = NULL, color = NULL) +
    theme_bw() + theme(legend.position = "bottom")
  save_fig(g, "s3_bootstrap_vs_jackknife", width = 3800, height = 1900)
}
if (has("models", "core_cv_params.csv")) {
  jk <- rd("models", "core_cv_params.csv") %>% filter(fold %in% FOLD_ORDER[1:10]) %>%
    select(model, weights, fold, V_A, V_sum) %>% pivot_longer(c(V_A, V_sum), names_to = "param") %>%
    mutate(model = factor(show_model(model), show_model(CORE_ORDER)), fold = factor(fold, rev(FOLD_ORDER)),
           weights = WEIGHT_LABELS[weights])
  g <- ggplot(jk, aes(x = value, y = fold, color = model)) + geom_point(size = 2.2) +
    facet_grid(param ~ weights) + labs(x = "Estimate with the block left out", y = "Block left out", color = NULL) +
    theme_bw()
  save_fig(g, "s3b_jackknife_by_block", width = 3000, height = 2000)
}

# ---- S4. profile over V_A ----------------------------------------------------------------------
if (has("profile", "profile_VA.csv")) {
  pf <- rd("profile", "profile_VA.csv") %>% mutate(weights = WEIGHT_LABELS[weights])
  grid <- pf %>% filter(point == "grid")
  marks <- pf %>% filter(point != "grid")
  g1 <- ggplot(grid, aes(x = V_A, y = qlr)) +
    geom_hline(yintercept = qchisq(0.95, 1), linetype = "dashed", color = "grey40") +
    geom_line(linewidth = 0.9) +
    geom_point(data = marks, aes(color = point), size = 3) +
    facet_grid(weights ~ model) + coord_cartesian(ylim = c(0, 60)) +
    scale_color_manual(values = c(joint = "#00BFC4", `two-step` = "#F8766D"), name = NULL) +
    labs(x = "V_A (held fixed; everything else refitted)",
         y = "Loss increase / residual variance\n(quasi-likelihood ratio; dashed: 3.84)",
         title = "How well is V_A pinned down? Profile over V_A") +
    theme_bw() + theme(legend.position = "bottom")
  save_fig(g1, "s4_profile_VA_loss", width = 3300, height = 1700)
  comp <- grid %>% select(model, weights, V_A, V_GxG, V_family) %>%
    pivot_longer(c(V_GxG, V_family), names_to = "component") %>% filter(!is.na(value)) %>%
    mutate(component = sub("V_", "", component))
  g2 <- ggplot(comp, aes(x = V_A, y = value, color = component)) + geom_line(linewidth = 0.9) +
    geom_abline(intercept = 1, slope = -1, linetype = "dotted", color = "grey50") +
    facet_grid(weights ~ model) +
    labs(x = "V_A (held fixed)", y = "Refitted variance component", color = NULL,
         title = "What the other components do as V_A is moved (dotted: sum of V = 1)") +
    theme_bw() + theme(legend.position = "bottom")
  save_fig(g2, "s4b_profile_VA_components", width = 3300, height = 1700)
}

# ---- S5. robustness: step-sibling kernel, one-at-a-time variants, phenotypes --------------------
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
  save_fig(g, "s5_stepsib_kernel", width = 2400, height = 1100)
  d <- ar %>% filter(!is.na(V_A)) %>%
    mutate(label = paste0(variant, " [", weighting, "]"),
           cohort_name = ifelse(cohort_name == "ALL", "All 82 categories", "Bio + in-law (Eftedal cohort)"))
  g <- ggplot(d, aes(x = V_A, y = reorder(label, V_A))) + geom_point(size = 2.5) +
    facet_wrap(~cohort_name) + labs(x = "h² (additive model)", y = NULL) + theme_bw()
  save_fig(g, "s5b_robustness_h2", width = 3000, height = 1500)
}
if (has("robustness", "phenotypes.csv")) {
  ph <- rd("robustness", "phenotypes.csv") %>%
    select(phenotype, setting, weighting, part, V_A) %>%
    pivot_wider(names_from = part, values_from = V_A) %>%
    mutate(gap = `bio only` - `bio + in-law`, setting = paste(setting, "|", weighting))
  g <- ggplot(ph, aes(x = gap, y = reorder(phenotype, gap), color = setting)) +
    geom_vline(xintercept = 0) + geom_point(size = 2.5) +
    labs(x = "h²(bio only) − h²(bio + in-law)", y = NULL, color = NULL,
         title = "In-law effect on h², every test") +
    theme_bw() + theme(legend.position = "bottom", legend.direction = "vertical")
  save_fig(g, "s5c_inlaw_all_phenotypes", width = 2200, height = 2000)
}
cat("figures written to", FIGURES, "\n")
