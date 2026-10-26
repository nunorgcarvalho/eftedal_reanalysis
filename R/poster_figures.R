# Poster figures (ASHG 2026). Run from anywhere:
#   Rscript R/poster_figures.R
# Writes PDF (+ PNG) to results/figures/poster/. Reuses the helper definitions at the top of
# R/make_figures.R (everything before its "M1" section) and R/plot_functions.R, with poster-sized
# text. Sizes: each figure is drawn at roughly half its printed size, so text in points here is
# about doubled on the poster.
#   P4  additive model: condensed fit (Eftedal et al. vs revised, Marginal Type) with the
#       V_A / m / r estimates by cohort as an inset; one shared legend at the bottom
#   P5  GxG and family environment: per-type component fit of A+GxG+F with its parameter
#       estimates as an inset; one shared legend at the bottom
#   P6  variance components of four models
#   supplement/s_components_<model>: the P5 plot for every core model (supplemental slides)
args <- commandArgs(trailingOnly = FALSE)
this <- sub("^--file=", "", args[grep("^--file=", args)])
REPO <- normalizePath(file.path(dirname(this), ".."))
mf <- readLines(file.path(REPO, "R", "make_figures.R"))
eval(parse(text = mf[seq_len(grep("^# ---- M1\\.", mf)[1] - 1)]))   # helpers only
FIGURES <- file.path(RESULTS, "figures", "poster")

BASE <- 15          # base font size (pt) at drawing size
save_poster <- function(gg, name, width, height) {   # width/height in inches
  dir.create(dirname(file.path(FIGURES, name)), showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(FIGURES, paste0(name, ".pdf")), gg, width = width, height = height, device = cairo_pdf)
  ggsave(file.path(FIGURES, paste0(name, ".png")), gg, width = width, height = height, dpi = 150)
}

# one colour per estimator (as in m1b); bio-only cohort drawn in a lighter shade of the same colour
EST_LEVELS <- c("Eftedal et al.", "Revised (Marginal Type)")
EST_SHOW <- c("Eftedal et al. (Marginal Pair)", "Revised (Marginal Type)")
# not blue/orange/green, which mark the variance components elsewhere on the poster
EST_COLS <- setNames(c("#C2185B", "#6A3D9A"), EST_LEVELS)
lighten <- function(col, f = 0.5) { m <- col2rgb(col); rgb(t(m + (255 - m) * f), maxColorValue = 255) }
OBS_LABEL <- "observed correlation (95% CI)"
WHISKER <- 0.08     # half-width of the observed-correlation CI caps (x units)

# [0, 1] y axes: major ticks/grid every 0.1, minor every 0.05
y01 <- function(...) scale_y_continuous(breaks = seq(0, 1, 0.1), minor_breaks = seq(0, 1, 0.05),
                                        expand = c(0, 0), guide = guide_axis(minor.ticks = TRUE), ...)
minor_y <- theme(panel.grid.minor.x = element_blank(),
                 panel.grid.minor.y = element_line(color = "grey93", linewidth = 0.3),
                 axis.minor.ticks.length.y = rel(0.5))

# relative-type labels drawn inside the plot, just above each observed point / bar
type_labels <- function(po, top) {
  po %>% left_join(top, by = "rel") %>%
    mutate(y = pmax(hi, top) + 0.012, label = str_wrap(CONDENSED_LABELS[rel], 9))
}
condensed_theme <- function() {
  theme_bw(base_size = BASE) + minor_y +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.position = "bottom", legend.box = "horizontal", legend.text = element_text(size = BASE * 0.85),
          legend.key.width = unit(1.6, "lines"), legend.margin = margin(0, 6, 0, 6),
          legend.spacing.x = unit(1.2, "lines"))
}

# ---- P4. additive model -------------------------------------------------------------------------
if (has("robustness", "additive_predictions.csv") && has("models", "core_predictions.csv") &&
    has("poster", "additive_by_cohort.csv")) {
  ## main: condensed fit
  ar <- rd("robustness", "additive_rows.csv")
  ap <- rd("robustness", "additive_predictions.csv")
  cp <- rd("models", "core_predictions.csv") %>% filter(model == "A")
  ie <- ar %>% filter(variant == "Eftedal settings", cohort_name == "ALL_but_adopted_UZ") %>% pull(fit_id)
  pr <- bind_rows(ap %>% filter(fit_id == ie) %>% mutate(model = EST_SHOW[1]),
                  cp %>% filter(weights == "hybrid_block") %>% mutate(model = EST_SHOW[2]))
  obs <- obs_from(cp)
  po <- pool_by_type(obs) %>% mutate(lo = corr - z95 * se, hi = corr + z95 * se)
  w <- obs %>% transmute(row, w = 1 / se^2)
  pp <- pr %>% filter(rel %in% CONDENSED_ORDER) %>% inner_join(w, by = "row") %>%
    group_by(model, rel) %>% summarize(pred = sum(w * pred) / sum(w), .groups = "drop") %>%
    mutate(x = match(rel, CONDENSED_ORDER), model = factor(model, EST_SHOW))
  path <- pp %>% reframe(x = c(x - 0.38, x + 0.38), pred = c(pred, pred), .by = c(model, rel)) %>% arrange(model, x)
  lab <- type_labels(po, pp %>% group_by(rel) %>% summarize(top = max(pred)))
  ## inset: V_A / m / r by estimator x cohort
  PARAM_LEVELS <- c(V_A = "V[A]", m = "m", r = "r")
  COHORTS <- c("bio + in-laws" = "all relative types", "bio only" = "biological relatives only")
  DX <- 0.16
  d <- rd("poster", "additive_by_cohort.csv") %>% filter(param %in% names(PARAM_LEVELS), label %in% EST_LEVELS) %>%
    mutate(est = factor(label, EST_LEVELS), param = factor(PARAM_LEVELS[param], PARAM_LEVELS),
           cohort = factor(COHORTS[cohort_name], COHORTS),
           x = 1.5 * as.numeric(est) + ifelse(cohort == COHORTS[1], -DX, DX),   # estimators 1.5 apart
           col = EST_COLS[as.character(est)],                            # cohort shown by marker shape
           box_fill = ifelse(cohort == COHORTS[2], "white", col))
  mt <- d %>% filter(method == "block jackknife") %>%
    mutate(clip_lo = lo < 0, clip_hi = hi > 1, lo = pmax(lo, 0), hi = pmin(hi, 1))
  mp <- d %>% filter(method == "parametric bootstrap") %>%
    mutate(tcol = col,
           tx = ifelse(cohort == COHORTS[1], x - 0.1, x + 0.1), hj = ifelse(cohort == COHORTS[1], 1, 0))
  arrow_end <- arrow(length = unit(0.08, "in"), type = "closed")
  inset <- ggplot() +
    geom_segment(data = mt, aes(x = x, xend = x, y = lo, yend = hi, color = col), linewidth = 1.3) +
    geom_segment(data = mt %>% filter(clip_lo), aes(x = x, xend = x, y = lo + 0.02, yend = lo, color = col),
                 linewidth = 1.3, arrow = arrow_end) +
    geom_segment(data = mt %>% filter(clip_hi), aes(x = x, xend = x, y = hi - 0.02, yend = hi, color = col),
                 linewidth = 1.3, arrow = arrow_end) +
    geom_rect(data = mp, aes(xmin = x - 0.04, xmax = x + 0.04, ymin = lo, ymax = hi, fill = box_fill, color = col),
              linewidth = 0.5) +
    # all relative types: filled circle; biological relatives only: X
    geom_point(data = mp %>% filter(cohort == COHORTS[1]), aes(x = x, y = estimate, fill = col),
               shape = 21, size = 3.6, stroke = 0.6, color = "grey15") +
    geom_point(data = mp %>% filter(cohort == COHORTS[2]), aes(x = x, y = estimate, color = col),
               shape = 4, size = 3.8, stroke = 1.8) +
    geom_text(data = mp, aes(x = tx, y = estimate, label = sprintf("%.2f", estimate), color = tcol, hjust = hj),
              size = BASE * 0.22, fontface = "bold") +
    scale_color_identity() + scale_fill_identity() +
    facet_wrap(~param, nrow = 1, labeller = label_parsed) +
    scale_x_continuous(breaks = NULL, expand = expansion(add = 0.5)) +
    y01() +
    coord_cartesian(ylim = c(0, 1), clip = "off") +
    labs(x = NULL, y = "Estimate") +
    theme_bw(base_size = BASE) +
    minor_y + theme(panel.grid.major.x = element_blank(),
          strip.text = element_text(size = BASE * 1.1), panel.spacing.x = unit(0.7, "lines"),
          plot.background = element_rect(fill = "white", color = "grey40", linewidth = 0.6),
          plot.margin = margin(6, 8, 8, 6))
  ## legend keys: estimator lines (colour), observed points (shape), cohort shades (size)
  ckey <- tibble(cohort = factor(COHORTS, COHORTS), col = c("grey35", "grey35"))
  g <- ggplot() +
    geom_hline(yintercept = 0) +
    bio_split(y = 0.30, size = BASE * 0.2) +
    geom_path(data = path, aes(x = x, y = pred, color = model, group = model), linewidth = 0.5, alpha = 0.8) +
    geom_segment(data = pp, aes(x = x - 0.38, xend = x + 0.38, y = pred, yend = pred, color = model),
                 linewidth = 2.2, alpha = 0.85) +
    geom_errorbar(data = po, aes(x = x, ymin = lo, ymax = hi), width = 2 * WHISKER, linewidth = 0.7) +
    geom_point(data = po, aes(x = x, y = corr, shape = OBS_LABEL), size = 2.6, key_glyph = "pointrange") +
    geom_text(data = lab, aes(x = x, y = y, label = label), size = BASE * 0.24, vjust = 0, lineheight = 0.85) +
    geom_point(data = ckey, aes(x = 1, y = 0.5, size = cohort, fill = col), shape = 21, alpha = 0) +
    annotation_custom(ggplotGrob(inset), xmin = match("half sibs", CONDENSED_ORDER) - 0.42,
                      xmax = length(CONDENSED_ORDER) + 0.55, ymin = 0.40, ymax = 0.995) +
    scale_color_manual(values = setNames(unname(EST_COLS), EST_SHOW), name = NULL,
                       guide = guide_legend(order = 1, ncol = 1, override.aes = list(linewidth = 2.2, alpha = 1))) +
    scale_shape_manual(values = c(16), name = NULL, guide = guide_legend(order = 2, override.aes = list(size = 0.55))) +
    scale_size_manual(values = c(4.5, 4.5), name = NULL,
                      guide = guide_legend(order = 3, ncol = 1, override.aes = list(shape = c(21, 4), fill = c("grey35", NA),
                                                                                    stroke = c(0.6, 1.8), alpha = 1))) +
    scale_fill_identity() +
    scale_x_continuous(breaks = NULL, expand = expansion(add = 0.6)) +
    y01() +
    coord_cartesian(ylim = c(-0.02, 1)) +
    labs(x = NULL, y = "Correlation") +
    condensed_theme()
  save_poster(g, "p4_fit_and_params", width = 11.6, height = 8.85)
}

# ---- P5. GxG and family environment --------------------------------------------------------------
#   P5  (section g) per-type component fit of A+GxG+F with its parameter estimates as an inset,
#       coloured by component; one shared legend at the bottom
#   P6  (section f) variance components of four models (A+GxG+F label in bold)
P5_MODELS <- c("A", "A>GxG", "A>F", "A+GxG+F")
# parameter colours: additive (blue), GxG (orange), family environment (green); secondary
# parameters of a component (m, r; r_f) in a lighter shade of its colour
PAR_LEVELS <- c(V_A = "V[A]", V_A0 = "V[A]^0", m = "m", r = "r", V_GxG = "V[GxG]", V_family = "V[C]", r_f = "c")
show_c <- function(x) gsub("\\bF\\b", "C", gsub(">", " \u2192 ", x), perl = TRUE)   # model labels: F -> C
PAR_COLS <- c(V_A = COMP_FILL[["A"]], V_A0 = COMP_FILL[["A"]], m = lighten(COMP_FILL[["A"]], 0.35), r = lighten(COMP_FILL[["A"]], 0.35),
              V_GxG = COMP_FILL[["GxG"]], V_family = COMP_FILL[["family"]], r_f = lighten(COMP_FILL[["family"]], 0.35))
PAR_TEXT <- c(V_A = COMP_LINE[["A"]], V_A0 = COMP_LINE[["A"]], m = COMP_LINE[["A"]], r = COMP_LINE[["A"]], V_GxG = COMP_LINE[["GxG"]],
              V_family = COMP_LINE[["family"]], r_f = COMP_LINE[["family"]])
SEG_SHOW <- c(A0 = "A (random mating)", A_AM = "A (excess under AM)", GxG = "GxG",
              fam_ext = "C (beyond-sibling)", fam_sib = "C (within-sibling)")
PAR_GROUP <- c(V_A = 1, V_A0 = 1, m = 1, r = 1, V_GxG = 2, V_family = 3, r_f = 3)
# per-type component fit of one model with its parameter estimates as an inset (circle; box =
# Marginal Pair CI; line = Marginal Type CI). Only the model's own parameters and components are
# shown. Used for P5 (A+GxG+F) and for the supplemental slides (every core model).
# preds / est / bt default to the core model's predictions, Marginal Type jackknife and Marginal
# Pair bootstrap intervals; pass them to draw another fit (e.g. Eftedal et al.'s).
components_and_params <- function(model_name, weights_name = "hybrid_block", preds = NULL, est = NULL, bt = NULL) {
  ## main: per-type predicted components
  if (is.null(preds))
    preds <- rd("models", "core_predictions.csv") %>% filter(model == model_name, weights == weights_name)
  obs <- obs_from(preds)
  po <- pool_by_type(obs) %>% mutate(lo = corr - z95 * se, hi = corr + z95 * se)
  parts <- c(A0 = "part_A0", A_AM = "part_A_AM", GxG = "comp_GxG",
             fam_ext = "part_fam_ext", fam_sib = "part_fam_sib")   # familial env.: extended + sib-specific
  parts <- parts[sapply(parts, function(p) p %in% names(preds) && any(!is.na(preds[[p]])))]   # the model's own
  w <- obs %>% transmute(row, w = 1 / se^2)
  pp <- preds %>% filter(rel %in% CONDENSED_ORDER) %>% inner_join(w, by = "row") %>%
    select(rel, w, all_of(unname(parts))) %>%
    pivot_longer(all_of(unname(parts)), names_to = "col", values_to = "value") %>%
    mutate(value = replace_na(value, 0), segment = factor(names(parts)[match(col, parts)], SEG_ORDER)) %>%
    group_by(rel, segment) %>% summarize(value = sum(w * value) / sum(w), .groups = "drop") %>%
    mutate(x = match(rel, CONDENSED_ORDER)) %>% arrange(x, segment) %>%
    group_by(rel) %>% mutate(ymax = cumsum(value), ymin = ymax - value) %>% ungroup() %>%
    mutate(xmin = x - 0.36, xmax = x + 0.36, signif = TRUE, weights = "all")
  stripes <- seg_stripes(pp, spacing = 0.02, slope = 0.25)   # A_AM rising, fam_sib falling
  lab <- type_labels(po, pp %>% group_by(rel) %>% summarize(top = max(ymax)))
  ## inset: the model's parameter estimates, grouped by component
  if (is.null(est))
    est <- rd("models", "core_jackknife.csv") %>% filter(model == model_name, weights == weights_name)
  est <- est %>% filter(param %in% names(PAR_LEVELS), !is.na(estimate))
  pars <- names(PAR_LEVELS)[names(PAR_LEVELS) %in% est$param]
  est <- est %>% mutate(clip_lo = lo < 0, clip_hi = hi > 1, lo = pmax(lo, 0), hi = pmin(hi, 1))
  if (is.null(bt))
    bt <- rd("bootstrap", "intervals.csv") %>% filter(model == model_name, weights == weights_name)
  bt <- bt %>% filter(param %in% pars)
  est <- est %>% mutate(x = match(param, pars), col = PAR_COLS[param], tcol = PAR_TEXT[param])
  bt <- bt %>% mutate(x = match(param, pars), col = PAR_COLS[param])
  seps <- which(diff(PAR_GROUP[pars]) != 0) + 0.5      # thin lines between component groups
  arrow_end <- arrow(length = unit(0.08, "in"), type = "closed")
  inset <- ggplot() +
    geom_vline(xintercept = seps, color = "grey80", linewidth = 0.4) +
    geom_segment(data = est, aes(x = x, xend = x, y = lo, yend = hi, color = col), linewidth = 1.5) +
    geom_segment(data = est %>% filter(clip_lo), aes(x = x, xend = x, y = lo + 0.02, yend = lo, color = col),
                 linewidth = 1.5, arrow = arrow_end) +
    geom_segment(data = est %>% filter(clip_hi), aes(x = x, xend = x, y = hi - 0.02, yend = hi, color = col),
                 linewidth = 1.5, arrow = arrow_end) +
    geom_rect(data = bt, aes(xmin = x - 0.06, xmax = x + 0.06, ymin = lo, ymax = hi, fill = col),
              color = "grey15", linewidth = 0.3) +
    geom_point(data = est, aes(x = x, y = estimate, fill = tcol), shape = 21, size = 3.6, stroke = 0.6, color = "grey15") +
    # labels of estimates at/near 0 lifted off the axis
    geom_text(data = est, aes(x = x + 0.12, y = pmax(estimate, 0.035), label = sprintf("%.2f", estimate), color = tcol),
              hjust = 0, size = BASE * 0.24, fontface = "bold") +
    scale_color_identity() + scale_fill_identity() +
    scale_x_continuous(breaks = seq_along(pars), labels = parse(text = PAR_LEVELS[pars]),
                       expand = expansion(add = 0.45)) +
    y01() +
    coord_cartesian(ylim = c(0, 1), clip = "off") +
    labs(x = NULL, y = "Estimate") +
    theme_bw(base_size = BASE) +
    minor_y + theme(panel.grid.major.x = element_blank(),
          axis.text.x = element_text(size = BASE * 1.15, color = "black"),
          plot.background = element_rect(fill = "white", color = "grey40", linewidth = 0.6),
          plot.margin = margin(8, 10, 6, 6))
  ggplot() +
    bio_split(y = 0.30, size = BASE * 0.2) +
    geom_rect(data = pp, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = segment),
              key_glyph = draw_key_striped) +
    {if (nrow(stripes)) geom_segment(data = stripes, aes(x = x, xend = xend, y = y, yend = yend),
                                     color = "white", linewidth = 0.35)} +
    geom_hline(yintercept = 0) +
    geom_errorbar(data = po, aes(x = x, ymin = lo, ymax = hi), width = 2 * WHISKER, linewidth = 0.7) +
    geom_point(data = po, aes(x = x, y = corr, shape = OBS_LABEL), size = 2.6, key_glyph = "pointrange") +
    geom_text(data = lab, aes(x = x, y = y, label = label), size = BASE * 0.205, vjust = 0, lineheight = 0.85) +
    annotation_custom(ggplotGrob(inset), xmin = match("half sibs", CONDENSED_ORDER) - 0.45,
                      xmax = length(CONDENSED_ORDER) + 0.55, ymin = 0.40, ymax = 0.995) +
    scale_fill_manual(values = SEG_FILL, labels = SEG_SHOW[names(parts)], breaks = names(parts), name = NULL,
                      guide = guide_legend(order = 1, nrow = 2)) +
    scale_shape_manual(values = c(16), name = NULL, guide = guide_legend(order = 2, override.aes = list(size = 0.55))) +
    scale_x_continuous(breaks = NULL, expand = expansion(add = 0.6)) +
    y01() +
    coord_cartesian(ylim = c(-0.02, 1)) +
    labs(x = NULL, y = "Correlation") +
    condensed_theme() + theme(legend.box = "vertical", legend.spacing.y = unit(0.2, "lines"))
}
if (has("models", "core_jackknife.csv") && has("models", "core_predictions.csv") &&
    has("bootstrap", "intervals.csv")) {
  g <- components_and_params("A+GxG+F")
  save_poster(g, "p5_components_and_params", width = 10.0, height = 9.1)

  ## P6: variance components of four models (Marginal Type), A+GxG+F in bold
  g6 <- stacked_components(rd("models", "core_jackknife.csv") %>%
                             filter(weights == "hybrid_block", model %in% P5_MODELS), P5_MODELS, title = NULL) +
    labs(subtitle = NULL, y = "Variance component") +
    scale_x_continuous(breaks = seq_along(P5_MODELS), labels = show_c(P5_MODELS)) +
    theme_bw(base_size = BASE) +
    y01(limits = c(0, 1)) + minor_y +
    theme(legend.position = "none", panel.grid.major.x = element_blank(),
          strip.background = element_blank(), strip.text = element_blank(),
          axis.text.x = element_text(size = BASE * 0.95, color = "black",
                                     face = ifelse(P5_MODELS == "A+GxG+F", "bold", "plain")),
          plot.margin = margin(8, 10, 6, 6))
  save_poster(g6, "p6_components_four_models", width = 5.4, height = 4.75)
}

# ---- supplemental slides: the P5 plot for every core model (Marginal Type) ------------------------
CORE_MODELS <- c("A", "A>GxG", "A+GxG", "A>F", "A+F", "A+GxG>F", "A+F>GxG", "A>GxG+F", "A+GxG+F")
model_tag <- function(m) gsub("\\+", "", gsub(">", "_to_", m))   # file tag: '>' -> '_to_', '+' dropped
if (has("models", "core_jackknife.csv") && has("models", "core_predictions.csv") &&
    has("bootstrap", "intervals.csv")) {
  for (m in CORE_MODELS)
    save_poster(components_and_params(m), file.path("supplement", paste0("s_components_", model_tag(m))),
                width = 8.8, height = 8.0)   # smaller than P5, so text is larger on a slide
  # fit statistics (Marginal Type) as LaTeX snippets, \input by ASHG/supplemental_slides/supplemental_slides.tex
  fs <- rd("models", "core_fits.csv") %>% filter(weights == "hybrid_block") %>%
    left_join(rd("models", "core_cv_summary.csv") %>% select(model, weights, cv_rmse), by = c("model", "weights")) %>%
    mutate(cv_rank = rank(cv_rmse))
  for (m in CORE_MODELS) {
    s <- fs %>% filter(model == .env$m)   # fits also have a column m
    writeLines(sprintf("%s & %.1f \\\\\n%s & %.3f \\\\\n%s & %d of %d \\\\",
                       "$\\chi^2/\\mathrm{df}$", s$chi2_per_df, "CV RMSE", s$cv_rmse,
                       "CV rank", as.integer(s$cv_rank), length(CORE_MODELS)),
               file.path(FIGURES, "supplement", paste0("stats_", model_tag(m), ".tex")))
  }
}
# ---- supplemental slides: Eftedal et al.'s estimator -> revised, step by step (additive model) ------
# rows: the estimators of scripts/02_corrections.py; panels: V_A, V_A0, m, r. Same marks as the
# insets: circle = estimate, box = Marginal Pair CI (parametric bootstrap), line = Marginal Type CI
# (leave-one-type-out jackknife; arrow = truncated at 0 or 1). chi2/df under each row label.
if (has("corrections", "steps.csv")) {
  st <- rd("corrections", "steps.csv") %>% filter(!aside, param %in% c("V_A", "V_A0", "m", "r"))
  steps <- st %>% distinct(label, order, chi2_per_df) %>% arrange(order) %>%
    mutate(show = c("Eftedal et al.", "Eftedal et al., corrected\ninverse-variance weights",
                    "Eftedal et al., raw scale", "Revised (Marginal Pair)", "Revised (Marginal Type)")[order + 1],
           show = sprintf("%s\nχ²/df = %.0f", show, chi2_per_df))
  d <- st %>% left_join(steps %>% select(label, show), by = "label") %>%
    mutate(y = factor(show, rev(steps$show)), panel = factor(PAR_LEVELS[param], PAR_LEVELS[c("V_A", "V_A0", "m", "r")]),
           col = PAR_COLS[param], tcol = PAR_TEXT[param])
  jk <- d %>% filter(method == "block jackknife") %>%
    mutate(clip_lo = lo < 0, clip_hi = hi > 1, lo = pmax(lo, 0), hi = pmin(hi, 1))
  bt <- d %>% filter(method == "parametric bootstrap")
  arrow_end <- arrow(length = unit(0.08, "in"), type = "closed")
  g <- ggplot(jk, aes(y = y)) +
    geom_segment(aes(x = lo, xend = hi, yend = y, color = col), linewidth = 1.5) +
    geom_segment(data = jk %>% filter(clip_lo), aes(x = lo + 0.03, xend = lo, yend = y, color = col),
                 linewidth = 1.5, arrow = arrow_end) +
    geom_segment(data = jk %>% filter(clip_hi), aes(x = hi - 0.03, xend = hi, yend = y, color = col),
                 linewidth = 1.5, arrow = arrow_end) +
    geom_rect(data = bt, aes(xmin = lo, xmax = hi, ymin = as.numeric(y) - 0.12, ymax = as.numeric(y) + 0.12, fill = col),
              color = "grey15", linewidth = 0.3) +
    geom_point(aes(x = estimate, fill = tcol), shape = 21, size = 3.6, stroke = 0.6, color = "grey15") +
    geom_text(aes(x = estimate, label = sprintf("%.2f", estimate), color = tcol), vjust = -1.05,
              size = BASE * 0.22, fontface = "bold") +
    scale_color_identity() + scale_fill_identity() +
    facet_wrap(~panel, nrow = 1, labeller = label_parsed) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"),
                       minor_breaks = c(0.25, 0.75), expand = expansion(add = 0.04)) +
    scale_y_discrete(expand = expansion(add = 0.6)) +
    coord_cartesian(clip = "off") +
    labs(x = "Estimate", y = NULL) +
    theme_bw(base_size = BASE) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
          strip.text = element_text(size = BASE * 1.15), strip.background = element_rect(fill = "grey95"),
          axis.text.y = element_text(color = "black", lineheight = 0.95), panel.spacing.x = unit(1, "lines"))
  save_poster(g, file.path("supplement", "s_corrections"), width = 9.6, height = 5.0)
}

# ---- supplemental slides: leave-one-type-out cross-validation of the core models (Marginal Type) ---
# left: each model's CV RMSE (square root of the CV score: the precision-weighted squared miss of
# the held-out correlations, averaged within and then over the 10 relative types); right: the RMSE
# for each held-out relative type. Same model order as the slides.
FOLD_SHOW <- c(twins = "twins", sibs = "full-sibs", `half sibs` = "half-sibs", `1st cousins` = "1st cousins",
               `half cousins` = "half cousins", `2nd cousins` = "2nd cousins", `step sibs` = "step-sibs",
               cocousins = "co-cousins", `adoptive sibs` = "adoptive sibs", `adoptive cousins` = "adoptive cousins")
if (has("models", "core_cv_summary.csv") && has("models", "core_cv_scores.csv")) {
  ylev <- rev(show_c(CORE_MODELS))
  cvs <- rd("models", "core_cv_summary.csv") %>% filter(weights == "hybrid_block") %>%
    mutate(y = factor(show_c(model), ylev), best = cv_rmse == min(cv_rmse))
  # right panel: mean held-out residual of each relative type, observed - predicted, weighted as the
  # CV score (precision-weighted within each relationship type, then averaged over the block's types).
  # Hatched: the model predicts exactly 0 for every category of the type, so the residual is just the
  # observed correlation (e.g. adoptees under a genetics-only model).
  sc <- rd("models", "core_cv_heldout.csv") %>% filter(weights == "hybrid_block", fold %in% names(FOLD_SHOW)) %>%
    mutate(w = 1 / se^2, e = corr - pred) %>%
    group_by(model, fold, rel) %>%
    summarize(res = sum(w * e) / sum(w), zero = all(abs(pred) < 1e-12), .groups = "drop") %>%
    group_by(model, fold) %>% summarize(res = mean(res), zero = all(zero), .groups = "drop") %>%
    mutate(label = ifelse(round(res, 3) == 0, "0.000", sub("^-", "\u2212", sprintf("%+.3f", res))),
           y = factor(show_c(model), ylev), x = factor(FOLD_SHOW[fold], FOLD_SHOW))
  LIM <- ceiling(max(abs(sc$res)) * 100) / 100
  RES_COLS <- c("#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B")
  tile_col <- scales::gradient_n_pal(RES_COLS)(scales::rescale(sc$res, from = c(-LIM, LIM)))
  sc <- sc %>% mutate(dark = abs(res) > 0.55 * LIM, tile_col = tile_col)   # white text on the darker tiles
  # diagonal hatch lines clipped to each zero-prediction tile (tile = unit square around (x, y))
  hatch <- sc %>% filter(zero) %>% mutate(xi = as.numeric(x), yi = as.numeric(y)) %>%
    reframe(t = seq(-0.9, 0.9, by = 0.15), .by = c(xi, yi, dark)) %>%
    mutate(x0 = pmax(-0.5, -0.5 - t), x1 = pmin(0.5, 0.5 - t)) %>% filter(x1 > x0) %>%
    transmute(dark, x = xi + x0, xend = xi + x1, y = yi + x0 + t, yend = yi + x1 + t)
  ytheme <- theme(axis.text.y = element_text(color = "black", size = BASE * 0.95,
                                             face = ifelse(levels(cvs$y) == show_c("A+GxG+F"), "bold", "plain")))
  g1 <- ggplot(cvs, aes(x = cv_rmse, y = y)) +
    geom_segment(aes(x = 0, xend = cv_rmse, yend = y), color = "grey75", linewidth = 0.8) +
    geom_point(aes(color = best), size = 3.6) +
    geom_text(aes(label = sprintf("%.3f", cv_rmse), color = best), hjust = -0.35, size = BASE * 0.24, fontface = "bold") +
    scale_color_manual(values = c(`TRUE` = "#C2185B", `FALSE` = "#1F3A5F"), guide = "none") +
    scale_x_continuous(limits = c(0, 0.115), breaks = c(0, 0.025, 0.05, 0.075, 0.1),
                       labels = c("0", "", "0.05", "", "0.10"), expand = expansion(add = c(0, 0.002))) +
    coord_cartesian(clip = "off") +
    labs(x = "CV RMSE (all types)", y = NULL, title = "Overall") +
    theme_bw(base_size = BASE) + ytheme +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
          plot.margin = margin(5.5, 22, 5.5, 5.5),
          plot.title = element_text(size = BASE, face = "bold", hjust = 0.5))
  g2 <- ggplot(sc, aes(x = x, y = y)) +
    geom_tile(aes(fill = res), color = "white", linewidth = 0.6) +
    geom_segment(data = hatch, aes(x = x, xend = xend, y = y, yend = yend,
                                   color = ifelse(dark, "white", "grey10")), inherit.aes = FALSE,
                 alpha = 0.55, linewidth = 0.45) +
    scale_color_identity() +
    geom_tile(fill = NA, color = "white", linewidth = 0.6) +   # borders over the hatching
    new_scale_color() +
    # diverging: blue = negative residual (over-predicts), red = positive (under-predicts)
    scale_fill_gradientn(colours = RES_COLS,
                         limits = c(-LIM, LIM), name = NULL,
                         breaks = c(-0.1, 0, 0.1),
                         labels = c("\u22120.10\nover-predicts", "0", "+0.10\nunder-predicts"),
                         guide = guide_colorbar(barwidth = unit(20, "lines"), barheight = unit(0.8, "lines"))) +
    # labels on a patch of the tile's own colour, so the hatching never crosses the text
    new_scale_fill() +
    geom_label(aes(label = label, color = dark, fill = tile_col), size = BASE * 0.25, label.size = 0,
               label.padding = unit(0.12, "lines"), label.r = unit(0, "lines")) +
    scale_fill_identity() +
    scale_color_manual(values = c(`TRUE` = "white", `FALSE` = "grey10"), guide = "none") +
    scale_x_discrete(expand = expansion(add = 0)) + scale_y_discrete(expand = expansion(add = 0)) +
    labs(x = "Held-out relative type", y = NULL, title = "Mean residual by held-out relative type (observed \u2212 predicted)") +
    theme_bw(base_size = BASE) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid = element_blank(),
          legend.position = "bottom", legend.text = element_text(size = BASE * 0.8, lineheight = 0.9),
          legend.margin = margin(0, 0, 0, 0),
          axis.text.x = element_text(angle = 40, hjust = 1, color = "black"),
          plot.title = element_text(size = BASE, face = "bold", hjust = 0.5))
  g <- cowplot::plot_grid(g1, g2, nrow = 1, rel_widths = c(0.4, 1), align = "h", axis = "tb")
  save_poster(g, file.path("supplement", "s_cv"), width = 11.0, height = 6.3)
}

# ---- supplemental slides: parameter sampling correlations (leave-one-type-out refits), per core model -
# the S7 pairs plot (R/make_figures.R's plot_param_pairs) without heading or footnote, key on the right
if (has("models", "core_param_corr.csv") && has("models", "core_cv_params.csv")) {
  cvp_all <- rd("models", "core_cv_params.csv") %>% filter(fold %in% names(FOLD_FILL)) %>%
    mutate(fold = factor(fold, names(FOLD_FILL)))
  fits_all <- rd("models", "core_fits.csv"); jk_all <- rd("models", "core_jackknife.csv")
  pc_all <- rd("models", "core_param_corr.csv")
  for (m in CORE_MODELS) {
    sel <- function(d) d %>% filter(model == .env$m, weights == "hybrid_block")
    est <- sel(fits_all)
    k <- sum(!is.na(unlist(est[1, intersect(PAIR_PARAMS, names(est))])))
    # text scaled to the panel size (6.5 in / k), as in the notebook version (base 9 at ~1.4 in panels)
    g <- plot_param_pairs(sel(cvp_all), est, sel(jk_all), sel(pc_all),
                          title = NULL, foot = FALSE, legend = "right", base = min(14, 9 * (6.5 / k) / 1.4))
    save_poster(g, file.path("supplement", paste0("s_paircorr_", model_tag(m))), width = 8.6, height = 6.5)
  }
}

# ---- supplemental slides: where each core model's CV error comes from (Marginal Type row of m3e) ----
# left: the CV score (mean over the 10 held-out relative types of the within-type weighted MSE) stacked
# by held-out type, the adoptive types at the bottom; label = CV RMSE.
# Right: CV RMSE over all 10 types vs without the two adoptive types.
if (has("models", "core_cv_scores.csv")) {
  ADOPT <- c("adoptive sibs", "adoptive cousins")
  STACK <- c(ADOPT, setdiff(names(FOLD_FILL), ADOPT))
  xlev <- show_c(CORE_MODELS)
  bold_x <- ifelse(xlev == show_c("A+GxG+F"), "bold", "plain")
  sc <- rd("models", "core_cv_scores.csv") %>% filter(weights == "hybrid_block", fold %in% names(FOLD_FILL)) %>%
    group_by(model) %>% mutate(contrib = score / n() * 1000) %>% ungroup() %>%
    mutate(model = factor(show_c(model), xlev), fold = factor(fold, rev(STACK)))
  tot <- sc %>% group_by(model) %>%
    summarize(total = sum(contrib), adopt = sum(contrib[fold %in% ADOPT]),
              rmse_all = sqrt(mean(score)), rmse_noadopt = sqrt(mean(score[!fold %in% ADOPT])), .groups = "drop") %>%
    mutate(label = sprintf("%.3f", rmse_all))   # CV RMSE
  xtheme <- theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
                  axis.text.x = element_text(angle = 35, hjust = 1, color = "black", face = bold_x),
                  plot.title = element_text(size = BASE, face = "bold"), legend.position = "bottom")
  g1 <- ggplot(sc, aes(x = model, y = contrib, fill = fold)) +
    geom_col(width = 0.72, color = "white", linewidth = 0.25) +
    geom_text(data = tot, aes(x = model, y = total, label = label), inherit.aes = FALSE,
              vjust = -0.2, size = BASE * 0.24, lineheight = 0.9) +
    scale_fill_manual(values = FOLD_FILL, labels = FOLD_LABELS, breaks = STACK, name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    labs(x = NULL, y = "CV score (×10⁻³)", title = "CV score, stacked by held-out type") +
    theme_bw(base_size = BASE) + xtheme +
    theme(legend.key.size = unit(0.9, "lines"), legend.text = element_text(size = BASE * 0.8),
          legend.key.spacing.x = unit(18, "pt")) +
    guides(fill = guide_legend(nrow = 3, byrow = TRUE))
  db <- tot %>% select(model, rmse_all, rmse_noadopt) %>%
    pivot_longer(c(rmse_all, rmse_noadopt), names_to = "set", values_to = "rmse") %>%
    mutate(set = factor(set, c("rmse_all", "rmse_noadopt"), c("all 10 types", "without the adoptive types")))
  g2 <- ggplot(db, aes(x = model, y = rmse)) +
    geom_line(aes(group = model), color = "grey60", linewidth = 0.7) +
    geom_point(aes(color = set, shape = set), size = 3.2) +
    scale_color_manual(values = c("all 10 types" = "#eb6834", "without the adoptive types" = "#184f95"), name = NULL) +
    scale_shape_manual(values = c("all 10 types" = 16, "without the adoptive types" = 17), name = NULL) +
    labs(x = NULL, y = "CV RMSE", title = "CV RMSE with and without\nthe adoptive types") +
    theme_bw(base_size = BASE) + xtheme +
    theme(legend.direction = "vertical", legend.text = element_text(size = BASE * 0.8))
  g <- cowplot::plot_grid(g1, g2, nrow = 1, rel_widths = c(1.35, 1), align = "h", axis = "tb") +
    theme(plot.background = element_rect(fill = "white", color = NA))
  save_poster(g, file.path("supplement", "s_cv_by_type"), width = 13.0, height = 6.2)
}

# Eftedal et al.'s additive fit (their estimator; biological relatives + in-laws), same plot
if (has("poster", "eftedal_predictions.csv") && has("poster", "additive_by_cohort.csv")) {
  ac <- rd("poster", "additive_by_cohort.csv") %>% filter(label == "Eftedal et al.", cohort_name == "bio + in-laws")
  save_poster(components_and_params("Eftedal et al.", preds = rd("poster", "eftedal_predictions.csv"),
                                    est = ac %>% filter(method == "block jackknife"),
                                    bt = ac %>% filter(method == "parametric bootstrap")),
              file.path("supplement", "s_components_eftedal"), width = 8.8, height = 8.0)
  writeLines(sprintf("%s & %.1f \\\\", "$\\chi^2/\\mathrm{df}$", ac$chi2_per_df[1]),
             file.path(FIGURES, "supplement", "stats_eftedal.tex"))
}
cat("poster figures written to", FIGURES, "\n")

# provenance stamp (see src/eftedal/provenance.py): which commit produced these figures
local({
  out <- file.path(RESULTS, "PROVENANCE"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
  git <- function(...) suppressWarnings(system2("git", c("-C", REPO, ...), stdout = TRUE, stderr = FALSE))
  dirty <- length(git("status", "--porcelain", "--", "src", "scripts", "R")) > 0
  writeLines(c(paste0("commit=", git("rev-parse", "HEAD")), paste0("dirty=", as.integer(dirty)),
               paste0("time=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S"))),
             file.path(out, sprintf("%s.txt", "poster_figures")))
})
