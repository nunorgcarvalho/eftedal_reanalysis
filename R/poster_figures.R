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
components_and_params <- function(model_name, weights_name = "hybrid_block") {
  ## main: per-type predicted components
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
  est <- rd("models", "core_jackknife.csv") %>%
    filter(model == model_name, weights == weights_name, param %in% names(PAR_LEVELS), !is.na(estimate))
  pars <- names(PAR_LEVELS)[names(PAR_LEVELS) %in% est$param]
  est <- est %>% mutate(clip_lo = lo < 0, clip_hi = hi > 1, lo = pmax(lo, 0), hi = pmin(hi, 1))
  bt <- rd("bootstrap", "intervals.csv") %>%
    filter(model == model_name, weights == weights_name, param %in% pars)
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
  # fit statistics (Marginal Type) as LaTeX snippets, \input by poster/supplemental_slides.tex
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
cat("poster figures written to", FIGURES, "\n")
