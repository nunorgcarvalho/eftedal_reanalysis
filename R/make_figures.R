# Figures for the analysis (scripts/01-06). Run from anywhere:
#   Rscript R/make_figures.R
# Writes PNG + PDF to results/figures/: main figures m*, supplementary figures s*. Each block is
# skipped if its inputs are missing. Terminology: "marginal type" = relationship types as the unit
# (block-balanced weights; block-jackknife CI); "marginal pair" = pairs as the unit (1/SE^2 weights;
# parametric-bootstrap CI). Intervals are marginal-type 95% CIs unless labelled otherwise.
args <- commandArgs(trailingOnly = FALSE)
this <- sub("^--file=", "", args[grep("^--file=", args)])
REPO <- normalizePath(file.path(dirname(this), ".."))
RESULTS <- file.path(REPO, "results")
FIGURES <- file.path(RESULTS, "figures")
source(file.path(REPO, "R", "plot_functions.R"))

rd <- function(...) as_tibble(fread(file.path(RESULTS, ...)))
has <- function(...) file.exists(file.path(RESULTS, ...))

WEIGHT_LABELS <- c(invvar = "Marginal Pair weights (1/SE²)", hybrid_block = "Marginal Type weights (block-balanced)",
                   invSqErr = "Eftedal (log, 1/SE²)")
CI_LABELS <- c("block jackknife" = "Marginal Type CI (block jackknife)",
               "parametric bootstrap" = "Marginal Pair CI (parametric bootstrap)")
wrap_y <- scale_y_discrete(labels = scales::label_wrap(30))   # long labels span several lines
CORE_ORDER <- c("A", "A>GxG", "A+GxG", "A+F", "A+GxG>F", "A+F>GxG", "A>GxG+F", "A+GxG+F")
CI_COLORS <- c("Marginal Type CI (block jackknife)" = "#D55E00", "Marginal Pair CI (parametric bootstrap)" = "#0072B2")
# both CIs drawn from the same point: the wide marginal-type CI first, the narrow marginal-pair CI on top
overlay_ci <- function(d, x = "estimate") {
  list(geom_errorbar(data = d %>% filter(method == "block jackknife"), orientation = "y", width = 0,
                     aes(xmin = lo_c, xmax = hi_c, color = CI_LABELS[method]), linewidth = 0.9),
       geom_errorbar(data = d %>% filter(method == "parametric bootstrap"), orientation = "y", width = 0,
                     aes(xmin = lo_c, xmax = hi_c, color = CI_LABELS[method]), linewidth = 2.2),
       scale_color_manual(values = CI_COLORS, breaks = names(CI_COLORS), name = "95% interval"))
}
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
         subtitle = "Coloured whiskers: 95% marginal-type CI of each component (drawn at the top of its segment)\nBlack diamond and whisker: sum of the components with its 95% marginal-type CI") +
    theme_bw() + theme(legend.position = "bottom", panel.grid.major.x = element_blank())
}

# ---- M1. Eftedal's estimator -> ours ----------------------------------------------------------
if (has("corrections", "steps.csv")) {
  st <- rd("corrections", "steps.csv")
  lev <- rev(unique(st$label[order(st$order)]))
  PL <- c(V_A = "h² (V_A)", m = "m (mate genetic corr.)", r = "r (mate phenotypic corr.)")
  d <- st %>% mutate(label = factor(label, lev), panel = factor(PL[param], PL),
                     lo_c = pmax(lo, 0), hi_c = pmin(hi, 1))
  pts <- d %>% distinct(label, panel, estimate)
  g <- ggplot(d, aes(y = label)) + overlay_ci(d) +
    geom_point(data = pts, aes(x = estimate), size = 2.6) +
    facet_wrap(~panel, scales = "free_x", nrow = 1) + wrap_y +
    labs(x = "Estimate", y = NULL, title = "Additive model: Eftedal et al.'s estimator vs. the revised estimator",
         caption = "Both CIs drawn from the same point; the marginal-pair CI is often narrower than the point itself") +
    theme_bw() + theme(legend.position = "bottom")
  save_fig(g, "m1_corrections", width = 3300, height = 1250)
  ch <- st %>% distinct(label, chi2_per_df) %>% mutate(label = factor(label, lev))
  g <- ggplot(ch, aes(x = chi2_per_df, y = label)) + geom_col(width = 0.6, fill = "grey45") +
    geom_text(aes(label = sprintf("%.0f", chi2_per_df)), hjust = -0.2, size = 3.4) +
    geom_vline(xintercept = 1, linetype = "dashed") + wrap_y +
    scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(x = "\u03c7²/df (1 = fits within sampling error)", y = NULL, title = "Misfit of the additive model") +
    theme_bw()
  save_fig(g, "m1c_chi2", width = 2000, height = 1100)
}
if (has("robustness", "additive_predictions.csv") && has("models", "core_predictions.csv")) {
  ar <- rd("robustness", "additive_rows.csv")
  ap <- rd("robustness", "additive_predictions.csv")
  cp <- rd("models", "core_predictions.csv") %>% filter(model == "A")
  cf <- rd("models", "core_fits.csv") %>% filter(model == "A")
  ie <- ar %>% filter(variant == "Eftedal settings", cohort_name == "ALL_but_adopted_UZ") %>% pull(fit_id)
  labs <- c(sprintf("Eftedal et al.: h²=%.2f", ar$V_A[which(ar$fit_id == ie)]),
            sprintf("Revised (Marginal Pair): h²=%.2f", cf$V_A[cf$weights == "invvar"]),
            sprintf("Revised (Marginal Type): h²=%.2f", cf$V_A[cf$weights == "hybrid_block"]))
  pr <- bind_rows(ap %>% filter(fit_id == ie) %>% mutate(model = labs[1]),
                  cp %>% filter(weights == "invvar") %>% mutate(model = labs[2], in_cohort = TRUE),
                  cp %>% filter(weights == "hybrid_block") %>% mutate(model = labs[3], in_cohort = TRUE))
  cols <- setNames(c("#F8766D", "#7CAE00", "#00BFC4"), labs)
  save_fig(plot_fits(obs_from(cp), pr, colors = cols), "m1b_eftedalfit_eftedal_vs_revised", height = 1850)
  save_fig(plot_fits_condensed(obs_from(cp), pr, colors = cols, title = "Additive model"),
           "m1b_condensedfit_eftedal_vs_revised", width = 2600, height = 1500)
}

# ---- M2. in-laws ----------------------------------------------------------------------------
if (has("corrections", "inlaw.csv")) {
  PN <- c("h² (additive model)", "Difference (bio only \u2212 with in-laws)")
  il <- rd("corrections", "inlaw.csv") %>%
    mutate(label = factor(label, rev(unique(label))),
           panel = factor(ifelse(quantity == "difference", PN[2], PN[1]), PN),
           quantity = factor(quantity, c("bio only", "bio + in-law", "difference")),
           lo_c = lo, hi_c = hi)
  dodge <- position_dodge(width = 0.5)
  pts <- il %>% distinct(label, panel, quantity, estimate)
  g <- ggplot(il, aes(y = label, group = quantity)) +
    geom_vline(data = tibble(panel = factor(PN[2], PN), v = 0), aes(xintercept = v), linetype = "dashed") +
    geom_errorbar(data = il %>% filter(method == "block jackknife"), orientation = "y", width = 0,
                  aes(xmin = lo_c, xmax = hi_c, color = CI_LABELS[method]), linewidth = 0.9, position = dodge) +
    geom_errorbar(data = il %>% filter(method == "parametric bootstrap"), orientation = "y", width = 0,
                  aes(xmin = lo_c, xmax = hi_c, color = CI_LABELS[method]), linewidth = 2.2, position = dodge) +
    geom_point(data = pts, aes(x = estimate, shape = quantity), size = 2.8, position = dodge) +
    scale_color_manual(values = CI_COLORS, breaks = names(CI_COLORS), name = "95% interval") +
    facet_wrap(~panel, scales = "free_x") + wrap_y +
    scale_shape_manual(values = c("bio only" = 16, "bio + in-law" = 17, difference = 15), name = NULL) +
    labs(x = "Estimate (difference: CIs paired across the two fits)", y = NULL,
         title = "Does adding in-laws change h²?",
         caption = "Both CIs drawn from the same point; the marginal-pair CI is often narrower than the point itself") +
    theme_bw() + theme(legend.position = "bottom", legend.box = "vertical")
  save_fig(g, "m2_inlaw", width = 3000, height = 1250)
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
  COMPARISONS <- list(twostep_GxG = c("A", "A>GxG", "A+GxG"),
                      joint_GxG_F = c("A", "A+GxG", "A+F", "A+GxG+F"))
  PAL <- c("#00BFC4", "#F8766D", "#7CAE00", "#C77CFF")
  for (wt in c("hybrid_block", "invvar")) {
    for (cmp in names(COMPARISONS)) {
      ms <- COMPARISONS[[cmp]]
      labs <- sapply(ms, function(mod) paste0(show_model(mod), ": ", fmt_par(fs %>% filter(model == !!mod, weights == !!wt),
                                                    c("V_A", "V_GxG", "V_family", "r_f", "m", "r"))))
      pr <- bind_rows(lapply(seq_along(ms), function(i)
        fp %>% filter(model == ms[i], weights == !!wt) %>% mutate(model = labs[i])))
      cols <- setNames(PAL[seq_along(ms)], labs)
      save_fig(plot_fits(obs, pr, colors = cols) +
                 labs(title = paste("Weighting:", WEIGHT_LABELS[[wt]])), height = 1950,
               sprintf("m3b_eftedalfit_%s_%s", cmp, wt))
      save_fig(plot_fits_condensed(obs, pr, colors = cols, title = paste("Weighting:", WEIGHT_LABELS[[wt]])),
               sprintf("m3b_condensedfit_%s_%s", cmp, wt), width = 2600, height = 1500)
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
}

# ---- S2c. V_A vs m across the core models --------------------------------------------------
if (has("models", "core_jackknife.csv")) {
  vm <- rd("models", "core_jackknife.csv") %>% filter(param %in% c("V_A", "m")) %>%
    select(weights, model, param, estimate, lo, hi) %>%
    pivot_wider(names_from = param, values_from = c(estimate, lo, hi)) %>%
    mutate(model = factor(show_model(model), show_model(CORE_ORDER)), weights = WEIGHT_LABELS[weights])
  g <- ggplot(vm, aes(x = estimate_V_A, y = estimate_m, color = model)) +
    geom_errorbar(aes(ymin = pmax(lo_m, 0), ymax = hi_m), width = 0, alpha = 0.5) +
    geom_errorbar(aes(xmin = pmax(lo_V_A, 0), xmax = pmin(hi_V_A, 1)), width = 0, alpha = 0.5, orientation = "y") +
    geom_point(size = 3.2) +
    ggrepel::geom_text_repel(aes(label = model), size = 3, show.legend = FALSE, max.overlaps = 30) +
    facet_wrap(~weights) +
    labs(x = "V_A (h²)", y = "m (mate genetic correlation)", color = NULL,
         title = "V_A vs m across the core models (95% marginal-type CIs)",
         subtitle = "Models whose step 1 is A share A's V_A and m exactly") +
    theme_bw() + theme(legend.position = "bottom")
  save_fig(g, "s2c_VA_vs_m_core", width = 3000, height = 1600)
}

# ---- S3. bootstrap vs jackknife; which block moves the estimates -----------------------------
if (has("bootstrap", "intervals.csv") && has("models", "core_jackknife.csv")) {
  iv <- bind_rows(rd("bootstrap", "intervals.csv") %>% mutate(method = CI_LABELS[["parametric bootstrap"]]),
                  rd("models", "core_jackknife.csv") %>% mutate(method = CI_LABELS[["block jackknife"]])) %>%
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
    geom_point(data = marks %>% mutate(point = show_model(point)), aes(color = point), size = 3) +
    facet_grid(weights ~ model) + coord_cartesian(ylim = c(0, 60)) +
    scale_color_brewer(palette = "Set1", name = NULL) +
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
    facet_wrap(~cohort_name) + wrap_y + labs(x = "h² (additive model)", y = NULL) + theme_bw()
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
