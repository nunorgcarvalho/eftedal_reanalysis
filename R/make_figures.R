# Draw every figure from the CSVs in results/. Run from anywhere:
#   Rscript R/make_figures.R
args <- commandArgs(trailingOnly = FALSE)
this <- sub("^--file=", "", args[grep("^--file=", args)])
REPO <- normalizePath(file.path(dirname(this), ".."))
RESULTS <- file.path(REPO, "results")
FIGURES <- file.path(RESULTS, "figures")
source(file.path(REPO, "R", "plot_functions.R"))

rd <- function(...) as_tibble(fread(file.path(RESULTS, ...)))
has <- function(...) file.exists(file.path(RESULTS, ...))

SETTING_LABELS <- c(
  "log-invSqErr" = "log scale, 1/SE² weights (Eftedal et al.)",
  "log-invvar"   = "log scale, delta-method weights",
  "lin-invvar"   = "linear scale, 1/SE² weights",
  "log-hybrid"   = "log scale, hybrid weights",
  "lin-hybrid"   = "linear scale, hybrid weights (NHGRI)",
  "lin-uniform"  = "linear scale, uniform weights")
COHORT_LABELS <- c(ALL_but_adopted_UZ = "Bio + in-law", ALL_but_adopted = "Bio + in-law + UZ",
                   bio_only = "Bio only", ALL = "All 82")

split_label <- function(df) {
  df %>% separate(spec_label, c("model", "cohort_l", "setting", "mr_l", "step_l"), sep = "\\|",
                  remove = FALSE)
}

# ---- 1. NHGRI poster figures, redrawn from the new code's refits of the 2025 specs ------
if (has("legacy", "nhgri_refit_predictions.csv")) {
  p <- rd("legacy", "nhgri_refit_predictions.csv")
  obs <- p %>% filter(fit_id == 0) %>% select(row, rel, side, sex, corr, se)
  lab <- function(id, txt) p %>% filter(fit_id == id) %>% mutate(model = txt)
  save_fig(plot_fits(obs, bind_rows(lab(3, "A (lin,hybrid)"), lab(0, "A (log,invSqErr)")),
                     colors = c("A (lin,hybrid)" = "#00BFC4", "A (log,invSqErr)" = "#F8766D")),
           "nhgri_fig1_reproduced")
  save_fig(plot_fits(obs, bind_rows(lab(3, "A (lin,hybrid)"), lab(7, "A, GxG (lin,hybrid)")),
                     colors = c("A (lin,hybrid)" = "#00BFC4", "A, GxG (lin,hybrid)" = "#7CAE00")),
           "nhgri_fig2_reproduced")
  save_fig(plot_fits(obs, bind_rows(lab(3, "A (lin,hybrid)"), lab(23, "A, family (lin,hybrid)")),
                     colors = c("A (lin,hybrid)" = "#00BFC4", "A, family (lin,hybrid)" = "#C77CFF")),
           "nhgri_fig3_reproduced")
}

# ---- 2. Fisherian fits: Eftedal's weighting vs corrected weighting ----------------------
if (has("sensitivity", "fisher_predictions.csv")) {
  fs <- rd("sensitivity", "fisher_summary.csv") %>% split_label()
  fp <- rd("sensitivity", "fisher_predictions.csv")
  obs <- fp %>% filter(fit_id == 0) %>% select(row, rel, side, sex, corr, se)
  pick <- function(cohort, setting, mr = "free", step = "eftedal") {
    fs %>% filter(cohort == !!cohort, setting == !!setting, mr == !!mr, stepsib == !!step) %>% pull(fit_id)
  }
  mk <- function(id, txt) fp %>% filter(fit_id == id) %>% mutate(model = txt)
  h <- function(id) sprintf("%.2f", fs$V_A[fs$fit_id == id])
  i1 <- pick("ALL_but_adopted_UZ", "log-invSqErr"); i2 <- pick("ALL_but_adopted_UZ", "log-invvar")
  i3 <- pick("ALL_but_adopted_UZ", "lin-invvar")
  cols <- setNames(c("#F8766D", "#00BFC4", "#7CAE00"),
                   c(paste0("log, 1/SE² weights [Eftedal]: h²=", h(i1)),
                     paste0("log, delta-method weights: h²=", h(i2)),
                     paste0("linear, 1/SE² weights: h²=", h(i3))))
  save_fig(plot_fits(obs, bind_rows(mk(i1, names(cols)[1]), mk(i2, names(cols)[2]))),
           "fisher_eftedal_vs_corrected_weights")
  save_fig(plot_fits(obs, bind_rows(mk(i1, names(cols)[1]), mk(i2, names(cols)[2]), mk(i3, names(cols)[3])),
                     colors = cols), "fisher_three_weightings")

  # h2 / m / r across all settings
  d <- fs %>% filter(stepsib == "eftedal" | cohort == "bio_only") %>%
    mutate(setting = factor(setting, levels = rev(names(SETTING_LABELS)), labels = rev(SETTING_LABELS)),
           cohort = factor(COHORT_LABELS[cohort], levels = COHORT_LABELS),
           r_shown = r) %>%
    select(cohort, setting, mr, V_A, m, r) %>%
    pivot_longer(c(V_A, m, r), names_to = "param") %>%
    mutate(param = factor(param, c("V_A", "m", "r"), c("h² (V_A)", "m (mate genetic corr.)", "r (mate phenotypic corr.)")))
  g <- ggplot(d, aes(x = value, y = setting, color = cohort, shape = mr)) +
    geom_point(size = 2.6, position = position_dodge(width = 0.5)) +
    facet_wrap(~param, scales = "free_x") +
    scale_shape_manual(values = c(free = 16, equilibrium = 1),
                       labels = c(free = "m, r separate", equilibrium = "m = h² r"), name = NULL) +
    labs(x = NULL, y = NULL, color = NULL) + theme_bw() + theme(legend.position = "bottom")
  save_fig(g, "fisher_params_by_setting", width = 3300, height = 1500)
}

# ---- 3. Step-sibling formula --------------------------------------------------------------
if (has("sensitivity", "models_summary.csv")) {
  ms <- rd("sensitivity", "models_summary.csv") %>% split_label()
  mp <- rd("sensitivity", "models_predictions.csv")
  step_obs <- mp %>% filter(fit_id == 0, rel == "step sibs") %>%
    summarize(corr = weighted.mean(corr, 1 / se^2), se = sqrt(1 / sum(1 / se^2)))
  sp <- mp %>% filter(rel == "step sibs") %>% group_by(fit_id) %>% summarize(pred = first(pred)) %>%
    inner_join(ms %>% select(fit_id, model, setting, mr, stepsib, V_A), by = "fit_id")
  g <- ggplot(sp, aes(x = pred, y = model, color = stepsib, shape = mr)) +
    annotate("rect", xmin = step_obs$corr - z95 * step_obs$se, xmax = step_obs$corr + z95 * step_obs$se,
             ymin = -Inf, ymax = Inf, alpha = 0.15) +
    geom_vline(xintercept = step_obs$corr, linetype = "dashed") +
    geom_point(size = 2.5, position = position_dodge(width = 0.6)) +
    facet_wrap(~setting, nrow = 1) +
    scale_color_manual(values = c(eftedal = "#F8766D", corrected = "#00BFC4"),
                       labels = c(eftedal = "m((1+m)/2)² [Eftedal]", corrected = "m((1+r)/2)² [corrected]"),
                       name = "Step-sibling kernel") +
    scale_shape_manual(values = c(free = 16, equilibrium = 1), name = "m, r") +
    labs(x = "Predicted step-sibling correlation (shaded: observed 95% CI)", y = NULL) +
    theme_bw() + theme(legend.position = "bottom")
  save_fig(g, "stepsib_predictions", width = 3600, height = 1500)

  dv <- ms %>% select(model, setting, mr, stepsib, V_A, V_GxG, m, r) %>%
    pivot_wider(names_from = stepsib, values_from = c(V_A, V_GxG, m, r)) %>%
    mutate(dV_A = V_A_corrected - V_A_eftedal)
  g <- ggplot(dv, aes(x = dV_A, y = model, color = setting, shape = mr)) +
    geom_vline(xintercept = 0) + geom_point(size = 2.5, position = position_dodge(width = 0.6)) +
    scale_shape_manual(values = c(free = 16, equilibrium = 1), name = "m, r") +
    labs(x = "Change in V_A from correcting the step-sibling kernel", y = NULL, color = NULL) +
    theme_bw()
  save_fig(g, "stepsib_delta_VA", width = 2600, height = 1400)

  # ---- 4. variance components by model ----------------------------------------------------
  comp <- ms %>% filter(stepsib == "corrected") %>%
    select(model, setting, mr, starts_with("V_")) %>%
    pivot_longer(starts_with("V_"), names_to = "component", values_to = "value") %>%
    filter(!is.na(value)) %>%
    mutate(component = sub("V_", "", component),
           component = factor(component, c("A", "GxG", "house1", "house2", "cousins1", "cousins2", "family")))
  g <- ggplot(comp, aes(x = model, y = value, fill = component)) +
    geom_col() + facet_grid(mr ~ setting) +
    scale_fill_brewer(palette = "Set2") +
    labs(y = "Variance component (correlation scale)", x = NULL, fill = NULL) +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_fig(g, "components_by_model", width = 3600, height = 1800)

  # corrected re-make of the three poster comparisons (lin, 1/SE^2, m = h2 r, corrected step sibs)
  pick_m <- function(model, setting, mr = "equilibrium", step = "corrected") {
    ms %>% filter(model == !!model, setting == !!setting, mr == !!mr, stepsib == !!step) %>% pull(fit_id)
  }
  obs <- mp %>% filter(fit_id == 0) %>% select(row, rel, side, sex, corr, se)
  mk <- function(id, txt) mp %>% filter(fit_id == id) %>% mutate(model = txt)
  for (mod in c("M02", "M04b", "M05b")) {
    a <- pick_m("M01", "lin-invvar"); b <- pick_m(mod, "lin-invvar")
    txt <- c(sprintf("A: V_A=%.2f", ms$V_A[ms$fit_id == a]),
             sprintf("%s: %s", mod, ms %>% filter(fit_id == b) %>%
                       select(V_A, V_GxG, V_family) %>% pivot_longer(everything()) %>% filter(!is.na(value)) %>%
                       mutate(s = sprintf("%s=%.2f", name, value)) %>% pull(s) %>% paste(collapse = ", ")))
    save_fig(plot_fits(obs, bind_rows(mk(a, txt[1]), mk(b, txt[2])),
                       colors = setNames(c("#00BFC4", "#C77CFF"), txt)),
             paste0("corrected_A_vs_", mod))
  }
}

# ---- 5. leave-one-type-out ------------------------------------------------------------------
if (has("sensitivity", "loo_summary.csv")) {
  lo <- rd("sensitivity", "loo_summary.csv")
  base <- lo %>% filter(left_out == "(none)") %>% select(key_id, V_A0 = V_A)
  d <- lo %>% filter(left_out != "(none)") %>% left_join(base, by = "key_id") %>%
    mutate(left_out = factor(left_out, levels = rev(unique(lo$left_out[lo$left_out != "(none)"]))))
  g <- ggplot(d, aes(x = V_A, y = left_out)) +
    geom_vline(aes(xintercept = V_A0), linetype = "dashed") +
    geom_point() + facet_wrap(~key_label, ncol = 5) +
    labs(x = "V_A with one relationship type left out (dashed: all types)", y = NULL) +
    theme_bw() + theme(strip.text = element_text(size = 6))
  save_fig(g, "leave_one_type_out", width = 4000, height = 2200)
}

# ---- 6. uncertainty ------------------------------------------------------------------------
if (has("uncertainty", "cis.csv")) {
  ci <- rd("uncertainty", "cis.csv") %>% filter(param %in% c("V_A", "V_GxG", "V_family", "m", "r"))
  g <- ggplot(ci, aes(x = estimate, xmin = lo, xmax = hi, y = method, color = method)) +
    geom_pointrange(size = 0.2) + facet_grid(key_label ~ param, scales = "free_x") +
    labs(x = "95% interval", y = NULL) + theme_bw() +
    theme(legend.position = "none", strip.text.y = element_text(angle = 0, size = 6))
  save_fig(g, "intervals_by_method", width = 4000, height = 3000)
}
if (has("uncertainty", "hessian_scaling.csv")) {
  hs <- rd("uncertainty", "hessian_scaling.csv") %>% filter(param == "V_A")
  g <- ggplot(hs, aes(x = factor, y = se, color = method)) +
    geom_point(size = 2.5) + geom_line() + scale_x_log10() + scale_y_log10() +
    facet_wrap(~spec) +
    labs(x = "Constant multiplying every weight (estimate unchanged)", y = "Reported SE of V_A",
         color = NULL) + theme_bw()
  save_fig(g, "hessian_se_vs_weight_scale", width = 2800, height = 1300)
}
if (has("uncertainty", "profile_r.csv")) {
  pr <- rd("uncertainty", "profile_r.csv")
  if (!"df_best" %in% names(pr)) pr$df_best <- 50   # 53 categories - 3 free parameters
  # quasi-likelihood ratio: loss increase in units of the residual variance (~ chi2_1)
  pr <- pr %>% mutate(qlr = (loss - loss_best) / (loss_best / df_best))
  g1 <- ggplot(pr, aes(x = r, y = qlr, color = spec)) + geom_line() + geom_point(size = 1) +
    geom_hline(yintercept = qchisq(0.95, 1), linetype = "dashed") +
    labs(x = "r (fixed)", y = "Quasi-likelihood ratio (dashed: 95% cut-off)") + theme_bw() + theme(legend.position = "bottom", legend.direction = "vertical")
  g2 <- ggplot(pr, aes(x = r, y = V_A, color = spec)) + geom_line() + geom_point(size = 1) +
    labs(x = "r (fixed)", y = "V_A") + theme_bw() + theme(legend.position = "none")
  save_fig(g1, "profile_r_loss", width = 2200, height = 1600)
  save_fig(g2, "profile_r_VA", width = 2200, height = 1300)
}
cat("figures written to", FIGURES, "\n")

# ---- 7. all phenotypes ---------------------------------------------------------------------
if (has("phenotypes", "summary.csv")) {
  ph <- rd("phenotypes", "summary.csv") %>%
    mutate(cohort = factor(COHORT_LABELS[cohort], levels = COHORT_LABELS),
           setting = factor(setting, levels = names(SETTING_LABELS), labels = SETTING_LABELS))
  g <- ggplot(ph %>% filter(mr == "free"), aes(x = V_A, y = reorder(phenotype, V_A), color = setting)) +
    geom_point(size = 2.3) + facet_wrap(~cohort) +
    labs(x = "h² (additive model, m and r separate)", y = NULL, color = NULL) +
    theme_bw() + theme(legend.position = "bottom", legend.direction = "vertical")
  save_fig(g, "phenotypes_h2", width = 2600, height = 2000)
  g <- ggplot(ph %>% filter(mr == "equilibrium", cohort == "Bio + in-law"),
              aes(x = r, y = reorder(phenotype, r), color = setting)) +
    geom_point(size = 2.3) +
    labs(x = "Mate phenotypic correlation r implied by m = h² r (bio + in-law cohort)", y = NULL, color = NULL) +
    theme_bw() + theme(legend.position = "bottom", legend.direction = "vertical")
  save_fig(g, "phenotypes_equilibrium_r", width = 2400, height = 1800)
}
cat("phenotype figures done\n")
