# Poster figures (ASHG 2026). Run from anywhere:
#   Rscript R/poster_figures.R
# Writes PDF (+ PNG) to results/figures/poster/. Reuses the helper definitions at the top of
# R/make_figures.R (everything before its "M1" section) and R/plot_functions.R, with poster-sized
# text. Sizes: each figure is drawn at roughly half its printed size, so text in points here is
# about doubled on the poster.
#   P4  additive model: condensed fit (Eftedal et al. vs revised, Marginal Type) with the
#       V_A / m / r estimates by cohort as an inset; one shared legend at the bottom
#   P5  GxG and family environment: per-type component fit of A+GxG+F with the variance
#       components of four models as an inset; one shared legend at the bottom
args <- commandArgs(trailingOnly = FALSE)
this <- sub("^--file=", "", args[grep("^--file=", args)])
REPO <- normalizePath(file.path(dirname(this), ".."))
mf <- readLines(file.path(REPO, "R", "make_figures.R"))
eval(parse(text = mf[seq_len(grep("^# ---- M1\\.", mf)[1] - 1)]))   # helpers only
FIGURES <- file.path(RESULTS, "figures", "poster")

BASE <- 15          # base font size (pt) at drawing size
save_poster <- function(gg, name, width, height) {   # width/height in inches
  dir.create(FIGURES, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(FIGURES, paste0(name, ".pdf")), gg, width = width, height = height, device = cairo_pdf)
  ggsave(file.path(FIGURES, paste0(name, ".png")), gg, width = width, height = height, dpi = 150)
}

# one colour per estimator (as in m1b); bio-only cohort drawn in a lighter shade of the same colour
EST_LEVELS <- c("Eftedal et al.", "Revised (Marginal Type)")
EST_SHOW <- c("Eftedal et al. (Marginal Pair)", "Revised (Marginal Type)")
EST_COLS <- setNames(c("#F8766D", "#00BFC4"), EST_LEVELS)
lighten <- function(col, f = 0.5) { m <- col2rgb(col); rgb(t(m + (255 - m) * f), maxColorValue = 255) }
OBS_LABEL <- "observed correlation (95% CI)"
WHISKER <- 0.08     # half-width of the observed-correlation CI caps (x units)

# relative-type labels drawn inside the plot, just above each observed point / bar
type_labels <- function(po, top) {
  po %>% left_join(top, by = "rel") %>%
    mutate(y = pmax(hi, top) + 0.012, label = str_wrap(CONDENSED_LABELS[rel], 9))
}
condensed_theme <- function() {
  theme_bw(base_size = BASE) +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
          panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
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
  DX <- 0.13
  d <- rd("poster", "additive_by_cohort.csv") %>% filter(param %in% names(PARAM_LEVELS), label %in% EST_LEVELS) %>%
    mutate(est = factor(label, EST_LEVELS), param = factor(PARAM_LEVELS[param], PARAM_LEVELS),
           cohort = factor(COHORTS[cohort_name], COHORTS),
           x = as.numeric(est) + ifelse(cohort == COHORTS[1], -DX, DX),
           col = ifelse(cohort == COHORTS[2], sapply(EST_COLS[as.character(est)], lighten),
                        EST_COLS[as.character(est)]))
  mt <- d %>% filter(method == "block jackknife") %>%
    mutate(clip_lo = lo < 0, clip_hi = hi > 1, lo = pmax(lo, 0), hi = pmin(hi, 1))
  mp <- d %>% filter(method == "parametric bootstrap") %>%
    mutate(tcol = ifelse(cohort == COHORTS[2], sapply(EST_COLS[as.character(est)], lighten, f = 0.25),
                         EST_COLS[as.character(est)]),
           tx = ifelse(cohort == COHORTS[1], x - 0.07, x + 0.07), hj = ifelse(cohort == COHORTS[1], 1, 0))
  arrow_end <- arrow(length = unit(0.08, "in"), type = "closed")
  inset <- ggplot() +
    geom_segment(data = mt, aes(x = x, xend = x, y = lo, yend = hi, color = col), linewidth = 1.3) +
    geom_segment(data = mt %>% filter(clip_lo), aes(x = x, xend = x, y = lo + 0.02, yend = lo, color = col),
                 linewidth = 1.3, arrow = arrow_end) +
    geom_segment(data = mt %>% filter(clip_hi), aes(x = x, xend = x, y = hi - 0.02, yend = hi, color = col),
                 linewidth = 1.3, arrow = arrow_end) +
    geom_rect(data = mp, aes(xmin = x - 0.04, xmax = x + 0.04, ymin = lo, ymax = hi, fill = col),
              color = "grey15", linewidth = 0.3) +
    geom_point(data = mp, aes(x = x, y = estimate, fill = col), shape = 21, size = 3.4, stroke = 0.6, color = "grey15") +
    geom_text(data = mp, aes(x = tx, y = estimate, label = sprintf("%.2f", estimate), color = tcol, hjust = hj),
              size = BASE * 0.22, fontface = "bold") +
    scale_color_identity() + scale_fill_identity() +
    facet_wrap(~param, nrow = 1, labeller = label_parsed) +
    scale_x_continuous(breaks = NULL, expand = expansion(add = 0.5)) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, 1), clip = "off") +
    labs(x = NULL, y = "Estimate") +
    theme_bw(base_size = BASE) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          strip.text = element_text(size = BASE * 1.1), panel.spacing.x = unit(0.7, "lines"),
          plot.background = element_rect(fill = "white", color = "grey40", linewidth = 0.6),
          plot.margin = margin(6, 8, 8, 6))
  ## legend keys: estimator lines (colour), observed points (shape), cohort shades (size)
  ckey <- tibble(cohort = factor(COHORTS, COHORTS), col = c("grey35", "grey75"))
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
    annotation_custom(ggplotGrob(inset), xmin = match("1st cousins", CONDENSED_ORDER) - 0.2,
                      xmax = length(CONDENSED_ORDER) + 0.55, ymin = 0.40, ymax = 0.995) +
    scale_color_manual(values = setNames(unname(EST_COLS), EST_SHOW), name = NULL,
                       guide = guide_legend(order = 1, ncol = 1, override.aes = list(linewidth = 2.2, alpha = 1))) +
    scale_shape_manual(values = c(16), name = NULL, guide = guide_legend(order = 2, override.aes = list(size = 0.55))) +
    scale_size_manual(values = c(4.5, 4.5), name = NULL,
                      guide = guide_legend(order = 3, ncol = 1, override.aes = list(fill = c("grey35", "grey75"), alpha = 1))) +
    scale_fill_identity() +
    scale_x_continuous(breaks = NULL, expand = expansion(add = 0.6)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1), expand = c(0, 0)) +
    coord_cartesian(ylim = c(-0.02, 1)) +
    labs(x = NULL, y = "Correlation") +
    condensed_theme()
  save_poster(g, "p4_fit_and_params", width = 13.5, height = 8.8)
}

# ---- P5. GxG and family environment --------------------------------------------------------------
P5_MODELS <- c("A", "A>GxG", "A>F", "A+GxG+F")
if (has("models", "core_jackknife.csv") && has("models", "core_predictions.csv")) {
  ## main: per-type predicted components of A+GxG+F (Marginal Type)
  preds <- rd("models", "core_predictions.csv") %>% filter(model == "A+GxG+F", weights == "hybrid_block")
  obs <- obs_from(preds)
  po <- pool_by_type(obs) %>% mutate(lo = corr - z95 * se, hi = corr + z95 * se)
  parts <- c(A0 = "part_A0", A_AM = "part_A_AM", GxG = "comp_GxG",
             fam_ext = "part_fam_ext", fam_sib = "part_fam_sib")   # familial env.: extended + sib-specific
  w <- obs %>% transmute(row, w = 1 / se^2)
  pp <- preds %>% filter(rel %in% CONDENSED_ORDER) %>% inner_join(w, by = "row") %>%
    select(rel, w, all_of(unname(parts))) %>%
    pivot_longer(all_of(unname(parts)), names_to = "col", values_to = "value") %>%
    mutate(segment = factor(names(parts)[match(col, parts)], SEG_ORDER)) %>%
    group_by(rel, segment) %>% summarize(value = sum(w * value) / sum(w), .groups = "drop") %>%
    mutate(x = match(rel, CONDENSED_ORDER)) %>% arrange(x, segment) %>%
    group_by(rel) %>% mutate(ymax = cumsum(value), ymin = ymax - value) %>% ungroup() %>%
    mutate(xmin = x - 0.36, xmax = x + 0.36, signif = TRUE, weights = "all")
  stripes <- seg_stripes(pp, spacing = 0.02, slope = 0.25)   # A_AM rising, fam_sib falling
  lab <- type_labels(po, pp %>% group_by(rel) %>% summarize(top = max(ymax)))
  ## inset: variance components of four models (Marginal Type)
  inset <- stacked_components(rd("models", "core_jackknife.csv") %>%
                                filter(weights == "hybrid_block", model %in% P5_MODELS), P5_MODELS, title = NULL) +
    labs(subtitle = NULL, y = "Variance component") +
    scale_x_continuous(breaks = seq_along(P5_MODELS), labels = gsub(">", " → ", P5_MODELS)) +
    theme_bw(base_size = BASE) +
    theme(legend.position = "none", panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
          strip.background = element_blank(), strip.text = element_blank(),
          axis.text.x = element_text(size = BASE * 0.95, color = "black"),
          plot.background = element_rect(fill = "white", color = "grey40", linewidth = 0.6),
          plot.margin = margin(8, 10, 6, 6))
  g <- ggplot() +
    bio_split(y = 0.30, size = BASE * 0.2) +
    geom_rect(data = pp, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = segment),
              key_glyph = draw_key_striped) +
    {if (nrow(stripes)) geom_segment(data = stripes, aes(x = x, xend = xend, y = y, yend = yend),
                                     color = "white", linewidth = 0.35)} +
    geom_hline(yintercept = 0) +
    geom_errorbar(data = po, aes(x = x, ymin = lo, ymax = hi), width = 2 * WHISKER, linewidth = 0.7) +
    geom_point(data = po, aes(x = x, y = corr, shape = OBS_LABEL), size = 2.6, key_glyph = "pointrange") +
    geom_text(data = lab, aes(x = x, y = y, label = label), size = BASE * 0.24, vjust = 0, lineheight = 0.85) +
    # '(Sum)' key, for the black outlines in the inset
    geom_rect(data = tibble(a = 1), aes(xmin = 1, xmax = 1, ymin = 2, ymax = 2, linetype = "(Sum)"),
              fill = NA, color = "black") +
    annotation_custom(ggplotGrob(inset), xmin = match("half sibs", CONDENSED_ORDER) - 0.45,
                      xmax = length(CONDENSED_ORDER) + 0.55, ymin = 0.40, ymax = 0.995) +
    scale_fill_manual(values = SEG_FILL, labels = c(A0 = "A (random mating)", A_AM = "A (excess under AM)", GxG = "GxG",
                                                    fam_ext = "family env. (extended)", fam_sib = "family env. (sib-specific)"),
                      breaks = SEG_ORDER[SEG_ORDER != "house"], name = NULL,
                      guide = guide_legend(order = 1, nrow = 2)) +
    scale_shape_manual(values = c(16), name = NULL, guide = guide_legend(order = 2, override.aes = list(size = 0.55))) +
    scale_linetype_manual(values = c("(Sum)" = "solid"), name = NULL,
                          guide = guide_legend(order = 3, override.aes = list(fill = NA, color = "black"))) +
    scale_x_continuous(breaks = NULL, expand = expansion(add = 0.6)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1), expand = c(0, 0)) +
    coord_cartesian(ylim = c(-0.02, 1)) +
    labs(x = NULL, y = "Correlation") +
    condensed_theme()
  save_poster(g, "p5_components_and_models", width = 13.5, height = 9.05)
}
cat("poster figures written to", FIGURES, "\n")
