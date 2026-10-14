# Poster figures (ASHG 2026). Run from anywhere:
#   Rscript R/poster_figures.R
# Writes PDF (+ PNG) to results/figures/poster/. Reuses the helper definitions at the top of
# R/make_figures.R (everything before its "M1" section) and R/plot_functions.R, with poster-sized
# text. Sizes: each figure is drawn at roughly half its printed size, so text in points here is
# about doubled on the poster.
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
EST_LEVELS <- c("Eftedal et al.", "Revised (Marginal Pair)", "Revised (Marginal Type)")
EST_COLS <- setNames(c("#F8766D", "#7CAE00", "#00BFC4"), EST_LEVELS)
lighten <- function(col, f = 0.5) { m <- col2rgb(col); rgb(t(m + (255 - m) * f), maxColorValue = 255) }

# ---- P1. condensed fit: Eftedal et al. vs revised (additive model) ---------------------------
if (has("robustness", "additive_predictions.csv") && has("models", "core_predictions.csv")) {
  ar <- rd("robustness", "additive_rows.csv")
  ap <- rd("robustness", "additive_predictions.csv")
  cp <- rd("models", "core_predictions.csv") %>% filter(model == "A")
  cf <- rd("models", "core_fits.csv") %>% filter(model == "A")
  ie <- ar %>% filter(variant == "Eftedal settings", cohort_name == "ALL_but_adopted_UZ") %>% pull(fit_id)
  labs_m <- c(sprintf("'Eftedal et al. (Marginal Pair): '*V[A]=='%.2f'", ar$V_A[which(ar$fit_id == ie)]),
              sprintf("'Revised (Marginal Pair): '*V[A]=='%.2f'", cf$V_A[cf$weights == "invvar"]),
              sprintf("'Revised (Marginal Type): '*V[A]=='%.2f'", cf$V_A[cf$weights == "hybrid_block"]))
  pr <- bind_rows(ap %>% filter(fit_id == ie) %>% mutate(model = labs_m[1]),
                  cp %>% filter(weights == "invvar") %>% mutate(model = labs_m[2]),
                  cp %>% filter(weights == "hybrid_block") %>% mutate(model = labs_m[3]))
  obs <- obs_from(cp)
  po <- pool_by_type(obs) %>%
    mutate(lo = corr - z95 * se, hi = corr + z95 * se, label = str_wrap(CONDENSED_LABELS[rel], 9))
  w <- obs %>% transmute(row, w = 1 / se^2)
  pp <- pr %>% filter(rel %in% CONDENSED_ORDER) %>% inner_join(w, by = "row") %>%
    group_by(model, rel) %>% summarize(pred = sum(w * pred) / sum(w), .groups = "drop") %>%
    mutate(x = match(rel, CONDENSED_ORDER))
  path <- pp %>% reframe(x = c(x - 0.38, x + 0.38), pred = c(pred, pred), .by = c(model, rel)) %>% arrange(model, x)
  cols <- setNames(unname(EST_COLS), labs_m)
  g <- ggplot() +
    geom_hline(yintercept = 0) +
    bio_split(y = 0.55, size = BASE * 0.24) +
    geom_path(data = path, aes(x = x, y = pred, color = model, group = model), linewidth = 0.5, alpha = 0.8) +
    geom_segment(data = pp, aes(x = x - 0.38, xend = x + 0.38, y = pred, yend = pred, color = model),
                 linewidth = 2.2, alpha = 0.85) +
    geom_errorbar(data = po, aes(x = x, ymin = lo, ymax = hi), width = 0.18, linewidth = 0.7) +
    geom_point(data = po, aes(x = x, y = corr), size = 2.6) +
    scale_x_continuous(breaks = seq_along(CONDENSED_ORDER), labels = CONDENSED_LABELS[CONDENSED_ORDER],
                       expand = expansion(add = 0.6)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1), expand = c(0, 0)) +
    scale_color_manual(values = cols, breaks = labs_m, labels = parse(text = labs_m), name = NULL) +
    coord_cartesian(ylim = c(-0.02, 1)) +
    labs(x = NULL, y = "Correlation") +
    theme_bw(base_size = BASE) +
    theme(axis.ticks.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1, color = "black", size = BASE * 0.85),
          legend.position = "inside", legend.position.inside = c(0.99, 0.985),
          legend.justification = c("right", "top"), legend.key.width = unit(1.6, "lines"),
          legend.background = element_rect(fill = "white", color = "grey40"))
  save_poster(g, "p1_condensedfit_additive", width = 8.6, height = 9.6)
}

# ---- P2. additive-model parameters by estimator x cohort -------------------------------------
# Estimate vertical. Circle: point estimate. Narrow rectangle behind it: Marginal Pair CI (bootstrap;
# only visible where it reaches past the circle). Line: Marginal Type CI (jackknife), arrowheads where
# it leaves [0, 1]. The two cohorts of an estimator sit side by side, dark = all relative types,
# light = biological relatives only.
if (has("poster", "additive_by_cohort.csv")) {
  PARAM_LEVELS <- c(V_A = "V[A]", m = "m", r = "r")
  EST_LABELS <- c("Eftedal et al." = "Eftedal et al.\n(Marginal Pair)",
                  "Revised (Marginal Pair)" = "Revised\n(Marginal Pair)",
                  "Revised (Marginal Type)" = "Revised\n(Marginal Type)")
  COHORTS <- c("bio + in-laws" = "all relative types", "bio only" = "biological relatives only")
  DX <- 0.13   # half the gap between the two cohorts of an estimator
  d <- rd("poster", "additive_by_cohort.csv") %>% filter(param %in% names(PARAM_LEVELS)) %>%
    mutate(est = factor(label, EST_LEVELS), param = factor(PARAM_LEVELS[param], PARAM_LEVELS),
           cohort = factor(COHORTS[cohort_name], COHORTS),
           x = as.numeric(est) + ifelse(cohort == COHORTS[1], -DX, DX),
           col = ifelse(cohort == COHORTS[2], sapply(EST_COLS[as.character(est)], lighten),
                        EST_COLS[as.character(est)]))
  mt <- d %>% filter(method == "block jackknife") %>%
    mutate(clip_lo = lo < 0, clip_hi = hi > 1, lo = pmax(lo, 0), hi = pmin(hi, 1))
  mp <- d %>% filter(method == "parametric bootstrap")
  arrow_end <- arrow(length = unit(0.09, "in"), type = "closed")
  key <- tibble(cohort = factor(COHORTS, COHORTS), col = c("grey35", "grey75"))
  ekey <- tibble(est = factor(EST_LABELS, EST_LABELS))
  g <- ggplot() +
    geom_segment(data = mt, aes(x = x, xend = x, y = lo, yend = hi, color = col), linewidth = 1.3) +
    geom_segment(data = mt %>% filter(clip_lo), aes(x = x, xend = x, y = lo + 0.02, yend = lo, color = col),
                 linewidth = 1.3, arrow = arrow_end) +
    geom_segment(data = mt %>% filter(clip_hi), aes(x = x, xend = x, y = hi - 0.02, yend = hi, color = col),
                 linewidth = 1.3, arrow = arrow_end) +
    geom_rect(data = mp, aes(xmin = x - 0.035, xmax = x + 0.035, ymin = lo, ymax = hi, fill = col),
              color = "grey15", linewidth = 0.3) +
    geom_point(data = mp, aes(x = x, y = estimate, fill = col), shape = 21, size = 3.4, stroke = 0.6, color = "grey15") +
    # legend for the two shades (grey stand-ins)
    geom_point(data = key, aes(x = 1, y = -1, shape = cohort, fill = col), size = 4.5) +
    scale_shape_manual(values = c(21, 21), name = NULL,
                       guide = guide_legend(order = 2, override.aes = list(fill = c("grey35", "grey75"), size = 4.5))) +
    # legend for the estimator colours (as in the fit figure)
    geom_point(data = ekey, aes(x = 1, y = -1, size = est), shape = 22) +
    scale_size_manual(values = rep(4.5, 3), name = NULL,
                      guide = guide_legend(order = 1, override.aes = list(fill = unname(EST_COLS), color = "grey15"))) +
    scale_color_identity() + scale_fill_identity() +
    facet_wrap(~param, nrow = 1, labeller = label_parsed) +
    scale_x_continuous(breaks = NULL, expand = expansion(add = 0.45)) +
    scale_y_continuous(breaks = seq(0, 1, 0.1), expand = expansion(add = 0.02)) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = NULL, y = "Estimate") +
    theme_bw(base_size = BASE) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          strip.text = element_text(size = BASE * 1.1), legend.position = "bottom",
          legend.box = "vertical", legend.spacing.y = unit(0.1, "lines"), legend.margin = margin(0, 0, 0, 0),
          panel.spacing.x = unit(0.8, "lines"))
  save_poster(g, "p2_additive_params", width = 8.6, height = 9.6)
}

# ---- P3. core models, Marginal Type only ------------------------------------------------------
if (has("models", "core_jackknife.csv")) {
  g <- stacked_components(rd("models", "core_jackknife.csv") %>% filter(weights == "hybrid_block"),
                          CORE_ORDER, title = NULL) +
    labs(subtitle = NULL, y = "Variance component") +
    scale_x_continuous(breaks = seq_along(CORE_ORDER), labels = gsub(">", "\n\u2192 ", CORE_ORDER)) +
    theme_bw(base_size = BASE) +
    theme(legend.position = "bottom", panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
          strip.background = element_blank(), strip.text = element_blank(),
          axis.text.x = element_text(size = BASE * 0.9, color = "black"))
  save_poster(g, "p3_components_core_type", width = 8.6, height = 8.6)
}
cat("poster figures written to", FIGURES, "\n")
