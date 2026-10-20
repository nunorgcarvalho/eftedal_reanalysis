# Figures for the analysis (scripts/01-06). Run from anywhere:
#   Rscript R/make_figures.R
# Writes PNG + PDF to results/figures/: main figures m*, supplementary figures s*. Each block is
# skipped if its inputs are missing. Terminology (always in this order): "Marginal Pair" = pairs as
# the unit (1/SE^2 weights; parametric-bootstrap CI); "Marginal Type" = relationship types as the
# unit (type-balanced weights; leave-one-type-out jackknife CI). Intervals are Marginal Type CIs
# unless labelled otherwise.
args <- commandArgs(trailingOnly = FALSE)
this <- sub("^--file=", "", args[grep("^--file=", args)])
REPO <- normalizePath(file.path(dirname(this), ".."))
RESULTS <- file.path(REPO, "results")
FIGURES <- file.path(RESULTS, "figures")
source(file.path(REPO, "R", "plot_functions.R"))

rd <- function(...) as_tibble(fread(file.path(RESULTS, ...)))
has <- function(...) file.exists(file.path(RESULTS, ...))

WEIGHT_LABELS <- c(invvar = "Marginal Pair", hybrid_block = "Marginal Type", invSqErr = "Eftedal (log, 1/SE²)")
CI_LABELS <- c("parametric bootstrap" = "Marginal Pair CI", "block jackknife" = "Marginal Type CI")
wrap_y <- scale_y_discrete(labels = scales::label_wrap(30))   # long labels span several lines
CORE_ORDER <- c("A", "A>GxG", "A+GxG", "A>F", "A+F", "A+GxG>F", "A+F>GxG", "A>GxG+F", "A+GxG+F")
CI_COLORS <- c("Marginal Pair CI" = "#0072B2", "Marginal Type CI" = "#D55E00")
# both CIs drawn from the same point: the wide marginal-type CI first, the narrow marginal-pair CI on top
overlay_ci <- function(d, x = "estimate") {
  list(geom_errorbar(data = d %>% filter(method == "block jackknife"), orientation = "y", width = 0,
                     aes(xmin = lo_c, xmax = hi_c, color = CI_LABELS[method]), linewidth = 0.9),
       geom_errorbar(data = d %>% filter(method == "parametric bootstrap"), orientation = "y", width = 0,
                     aes(xmin = lo_c, xmax = hi_c, color = CI_LABELS[method]), linewidth = 2.2),
       scale_color_manual(values = CI_COLORS, breaks = names(CI_COLORS), name = "95% interval"))
}
JOINT_ORDER <- c("A", "A+H", "A+F", "A+H+F", "A+GxG", "A+GxG+H", "A+GxG+F", "A+GxG+H+F")
show_model <- function(x) gsub("\\bF\\b", "C", gsub(">", " \u2192 ", x), perl = TRUE)   # shared env. F shown as C
# parameter symbols as plotmath strings (for label_parsed / parse(text = ...))
PAR_MATH <- c(V_A = "V[A]", V_A0 = "V[A]^0", m = "m", r = "r", V_GxG = "V[GxG]", V_family = "V[C]",
              V_house = "V[H]", V_sum = "Sigma*V", r_f = "c", pi = "pi", V_fam_ext = "V[C]~c",
              V_fam_sib = "V[C]*(1-c)")
par_math <- function(p) unname(ifelse(p %in% names(PAR_MATH), PAR_MATH[p], p))
COMP_ORDER <- c("A", "GxG", "house", "family")
COMP_FILL <- c(A = "#4C8CB5", GxG = "#E3894B", house = "#9B7BB8", family = "#6FAF7A")
# bar segments: the additive component split into its random-mating part V_A0 = V_A(1 - m) and the
# excess due to assortative mating, V_A m; the shared environment C (V_family, decay r_f = c) split
# into its beyond-sibling part V_C c (shared by relatives of family degree >= 1) and the
# within-sibling part V_C (1 - c)
SEG_ORDER <- c("A0", "A_AM", "GxG", "house", "fam_ext", "fam_sib")
# A_AM / fam_sib: same colour as A0 / fam_ext, distinguished by white diagonal stripes (/ and \
# respectively); their hex differs by one unit so the legend key can recognise them
SEG_FILL <- c(A0 = "#4C8CB5", A_AM = "#4C8CB6", GxG = "#E3894B", house = "#9B7BB8",
              fam_ext = "#6FAF7A", fam_sib = "#6FAF7B")
STRIPE_DIR <- c("#4C8CB6" = 1, "#6FAF7B" = -1)   # fill -> stripe direction (slope sign)
SEG_LABELS <- c(A0 = "A (random mating)", A_AM = "A (excess under AM)", GxG = "GxG", house = "H (household)",
                fam_ext = "C (beyond-sibling)", fam_sib = "C (within-sibling)")
COMP_LINE <- c(A = "#1F4E6E", GxG = "#9A4A14", house = "#5B3F78", family = "#2F6B3A")
FOLD_ORDER <- c("twins", "sibs", "half sibs", "1st cousins", "half cousins", "2nd cousins",
                "step sibs", "cocousins", "adoptive sibs", "adoptive cousins",
                "all in-laws", "all adoptive")
# held-out blocks: blue ramp = biological, orange ramp = non-biological (darker = closer); shape
# separates the two groups as well
FOLD_FILL <- c(twins = "#0d366b", sibs = "#184f95", `half sibs` = "#256abf", `1st cousins` = "#3987e5",
               `half cousins` = "#6da7ec", `2nd cousins` = "#86b6ef",
               `step sibs` = "#9a3a12", cocousins = "#c94f1f", `adoptive sibs` = "#eb6834",
               `adoptive cousins` = "#f29466")
FOLD_SHAPE <- setNames(c(rep(21, 6), rep(24, 4)), names(FOLD_FILL))
FOLD_LABELS <- c(twins = "twins", sibs = "full-sibs", `half sibs` = "half-sibs", `1st cousins` = "1st cousins",
                 `half cousins` = "half cousins", `2nd cousins` = "2nd cousins", `step sibs` = "step-sibs",
                 cocousins = "co-cousins", `adoptive sibs` = "adoptive sibs", `adoptive cousins` = "adoptive cousins")
obs_from <- function(p) p %>% distinct(row, .keep_all = TRUE) %>% select(row, rel, side, sex, corr, se)
# parameters and their values as one plotmath string, e.g. list(V[A]=='0.42', m=='0.25')
fmt_par <- function(d, ps) {
  v <- d %>% select(any_of(ps)) %>% pivot_longer(everything()) %>% filter(!is.na(value))
  sprintf("list(%s)", paste(sprintf("%s=='%.2f'", par_math(v$name), v$value), collapse = ", "))
}
# "<model>:  <parameters>" as a plotmath string
model_par_label <- function(model, d, ps) sprintf("'%s:'~~%s", show_model(model), fmt_par(d, ps))

# Interval glyphs: vertical line, flat caps at ends inside the parameter range, arrowheads at
# ends that were clipped because the interval runs past the range (e.g. a negative variance).
ci_glyphs <- function(d, colour_aes = TRUE, color = "black", linewidth = 0.7, cap = 0.035) {
  # d: x, ymin, ymax, clip_lo, clip_hi (+ component if colour_aes)
  ae <- function(...) if (colour_aes) aes(..., color = component) else aes(...)
  fix <- if (colour_aes) list() else list(color = color)
  c(list(do.call(geom_segment, c(list(data = d, mapping = ae(x = x, xend = x, y = ymin, yend = ymax),
                                      linewidth = linewidth, show.legend = FALSE), fix)),
         do.call(geom_segment, c(list(data = d %>% filter(!clip_lo), mapping = ae(x = x - cap, xend = x + cap, y = ymin, yend = ymin),
                                      linewidth = linewidth, show.legend = FALSE), fix)),
         do.call(geom_segment, c(list(data = d %>% filter(!clip_hi), mapping = ae(x = x - cap, xend = x + cap, y = ymax, yend = ymax),
                                      linewidth = linewidth, show.legend = FALSE), fix)),
         do.call(geom_segment, c(list(data = d %>% filter(clip_lo), mapping = ae(x = x, xend = x, y = ymin + 0.03, yend = ymin),
                                      linewidth = linewidth, arrow = arrow(length = unit(0.07, "in"), type = "closed"),
                                      show.legend = FALSE), fix)),
         do.call(geom_segment, c(list(data = d %>% filter(clip_hi), mapping = ae(x = x, xend = x, y = ymax - 0.03, yend = ymax),
                                      linewidth = linewidth, arrow = arrow(length = unit(0.07, "in"), type = "closed"),
                                      show.legend = FALSE), fix))))
}

# Stacked variance components: black outline = the sum, with its marginal-type CI in the middle of
# the bar; each component's CI is drawn at the top of its segment, staggered away from the centre
# in the order the components are stacked (left, right, further left, further right).
# Diagonal stripes inside rectangles (data coordinates): lines y = slope * x + c, clipped to each
# rectangle; used for the AM-excess part of the additive component.
stripe_segments <- function(r, spacing = 0.035, slope = 0.6) {
  if (nrow(r) == 0) return(tibble())
  bind_rows(lapply(seq_len(nrow(r)), function(i) {
    xmin <- r$xmin[i]; xmax <- r$xmax[i]; ymin <- r$ymin[i]; ymax <- r$ymax[i]
    if (ymax - ymin <= 1e-6) return(NULL)
    cr <- range(c(ymin, ymax) - slope * rep(c(xmin, xmax), each = 2))   # intercepts touching the rectangle
    out <- lapply(seq(cr[1], cr[2], by = spacing), function(c0) {
      xy <- sort(c((ymin - c0) / slope, (ymax - c0) / slope))
      xa <- max(xmin, xy[1]); xb <- min(xmax, xy[2])
      if (xb <= xa) return(NULL)
      tibble(x = xa, xend = xb, y = slope * xa + c0, yend = slope * xb + c0)
    })
    bind_rows(out) %>% mutate(weights = r$weights[i], signif = r$signif[i])
  }))
}
# stripes for every striped segment: A_AM rising (/), fam_sib falling (\)
seg_stripes <- function(seg, spacing = 0.035, slope = 0.6) {
  bind_rows(stripe_segments(seg %>% filter(segment == "A_AM"), spacing, slope),
            stripe_segments(seg %>% filter(segment == "fam_sib"), spacing, -slope))
}
# legend key: plain filled square, with stripes for the AM-excess fill
draw_key_striped <- function(data, params, size) {
  rect <- grid::rectGrob(gp = grid::gpar(col = NA, fill = scales::alpha(data$fill %||% "grey35", data$alpha %||% 1)))
  dir <- STRIPE_DIR[toupper(data$fill %||% "")]
  if (is.na(dir)) return(rect)
  s <- seq(-1, 1, by = 0.28)
  y0 <- pmax(0, s); y1 <- pmin(1, 1 + s)
  if (dir < 0) { y0 <- 1 - y0; y1 <- 1 - y1 }
  grid::grobTree(rect, grid::segmentsGrob(x0 = pmax(0, -s), y0 = y0, x1 = pmin(1, 1 - s), y1 = y1,
                                          gp = grid::gpar(col = "white", lwd = 1.3)),
                 vp = grid::viewport(clip = "on"))
}

YTOP <- 1   # top of the component plots; CIs reaching past it are capped with an arrowhead
stacked_components <- function(jk, order, title) {
  d <- jk %>% filter(param %in% paste0("V_", COMP_ORDER)) %>%
    mutate(component = factor(sub("V_", "", param), COMP_ORDER),
           model = factor(model, order), xm = as.numeric(model)) %>%
    filter(!is.na(estimate)) %>%
    arrange(weights, model, component) %>%
    group_by(weights, model) %>%
    mutate(top = cumsum(estimate), k = row_number(),
           x = xm + ifelse(k %% 2 == 1, -1, 1) * 0.11 * ceiling(k / 2),
           ymin = top - (estimate - pmax(lo, 0)), ymax_full = top + (pmin(hi, 1) - estimate),
           clip_lo = lo < 0, clip_hi = hi > 1 | ymax_full > YTOP, ymax = pmin(ymax_full, YTOP)) %>%
    ungroup() %>% mutate(weights = factor(weights, c("invvar", "hybrid_block")))
  s <- jk %>% filter(param == "V_sum") %>%
    mutate(model = factor(model, order), x = as.numeric(model), clip_lo = lo < 0, clip_hi = hi > 1,
           ymin = pmax(lo, 0), ymax = pmin(hi, 1), weights = factor(weights, c("invvar", "hybrid_block")))
  d <- d %>% mutate(signif = lo > 0)          # CI excludes 0 -> full colour; otherwise faded
  # bar segments: A split into V_A0 and the AM excess when V_A0 is available (CI stays on total V_A)
  va0 <- jk %>% filter(param == "V_A0") %>% select(weights, model, V_A0 = estimate) %>%
    mutate(weights = factor(weights, c("invvar", "hybrid_block")), model = factor(model, order))
  # likewise the familial component into its extended and sib-specific parts (CI on total V_fam)
  fsib <- jk %>% filter(param == "V_fam_sib") %>% select(weights, model, V_fam_sib = estimate) %>%
    mutate(weights = factor(weights, c("invvar", "hybrid_block")), model = factor(model, order))
  seg <- d %>% select(weights, model, xm, component, estimate, signif) %>%
    left_join(va0, by = c("weights", "model")) %>% left_join(fsib, by = c("weights", "model"))
  seg <- bind_rows(
    seg %>% filter(component == "A") %>%
      reframe(segment = c("A0", "A_AM"),
              estimate = if (is.na(V_A0)) c(estimate, 0) else c(V_A0, estimate - V_A0),
              .by = c(weights, model, xm, signif)),
    seg %>% filter(component == "family") %>%
      reframe(segment = c("fam_ext", "fam_sib"),
              estimate = if (is.na(V_fam_sib)) c(estimate, 0) else c(estimate - V_fam_sib, V_fam_sib),
              .by = c(weights, model, xm, signif)),
    seg %>% filter(!component %in% c("A", "family")) %>%
      transmute(weights, model, xm, signif, segment = as.character(component), estimate)) %>%
    mutate(segment = factor(segment, SEG_ORDER))
  seg <- seg %>% arrange(weights, model, segment) %>% group_by(weights, model) %>%
    mutate(ymax = cumsum(estimate), ymin = ymax - estimate, xmin = xm - 0.35, xmax = xm + 0.35) %>% ungroup()
  stripes <- seg_stripes(seg)
  ggplot() +
    geom_rect(data = seg, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = segment, alpha = signif),
              key_glyph = draw_key_striped) +
    {if (nrow(stripes)) geom_segment(data = stripes, aes(x = x, xend = xend, y = y, yend = yend, alpha = signif),
                                     color = "white", linewidth = 0.45, show.legend = FALSE)} +
    geom_rect(data = s, aes(xmin = x - 0.35, xmax = x + 0.35, ymin = 0, ymax = estimate, linetype = "(Sum)"),
              fill = NA, color = "black", linewidth = 0.6) +
    ci_glyphs(d %>% select(weights, component, x, ymin, ymax, clip_lo, clip_hi)) +
    ci_glyphs(s %>% select(weights, x, ymin, ymax, clip_lo, clip_hi), colour_aes = FALSE, linewidth = 1.1, cap = 0.05) +
    geom_point(data = s, aes(x = x, y = estimate), shape = 18, size = 3.2) +
    facet_wrap(~weights, ncol = 1, labeller = as_labeller(WEIGHT_LABELS)) +
    scale_x_continuous(breaks = seq_along(order), labels = show_model(order)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1), expand = c(0, 0)) +
    scale_fill_manual(values = SEG_FILL, labels = SEG_LABELS, breaks = SEG_ORDER, name = NULL,
                      guide = guide_legend(order = 1, nrow = 2)) +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.35), guide = "none") +
    scale_linetype_manual(values = c("(Sum)" = "solid"), name = NULL,
                          guide = guide_legend(order = 2, override.aes = list(fill = NA, color = "black"))) +
    scale_color_manual(values = COMP_LINE) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = "Variance component (proportion of phenotypic variance)", title = title,
         subtitle = paste0("Black outline and diamond: sum of the components, with its 95% Marginal Type CI (black, centre)\n",
                           "Coloured: each component's 95% Marginal Type CI, at the top of its segment\n",
                           "A: random-mating part (solid) + excess under assortative mating (striped); CI for the total\n",
                           "C: beyond-sibling part (solid) + within-sibling part (striped); CI for the total\n",
                           "Faded: the component's CI includes 0. Arrowhead: the CI runs past [0, 1]")) +
    theme_bw() + theme(legend.position = "bottom", panel.grid.major.x = element_blank(),
                       panel.grid.minor.x = element_blank(), panel.spacing = unit(1.2, "lines"),
                       plot.margin = margin(10, 10, 5, 5))
}

# Condensed fit-component plot for one model: per relationship type, the inverse-variance pooled
# observed correlation (point, 95% CI) and the predicted correlation as a stacked bar of its parts
# (V_A0 * a_k under random mating, the AM excess of the additive part, GxG, house, familial env.),
# each pooled with the same weights as the observations.
plot_condensed_components <- function(preds, title = NULL) {
  obs <- obs_from(preds)
  po <- pool_by_type(obs) %>%
    mutate(lo = corr - z95 * se, hi = corr + z95 * se, label = str_wrap(CONDENSED_LABELS[rel], 10))
  parts <- c(A0 = "part_A0", A_AM = "part_A_AM", GxG = "comp_GxG", house = "comp_house", fam_ext = "part_fam_ext", fam_sib = "part_fam_sib")
  parts <- parts[parts %in% names(preds)]
  parts <- parts[sapply(parts, function(cl) any(!is.na(preds[[cl]])))]   # components in this model
  w <- obs %>% transmute(row, w = 1 / se^2)
  pp <- preds %>% filter(rel %in% CONDENSED_ORDER) %>% inner_join(w, by = "row") %>%
    select(rel, w, all_of(unname(parts))) %>%
    pivot_longer(all_of(unname(parts)), names_to = "col", values_to = "value") %>%
    mutate(segment = factor(names(parts)[match(col, parts)], SEG_ORDER)) %>%
    group_by(rel, segment) %>% summarize(value = sum(w * value) / sum(w), .groups = "drop") %>%
    mutate(x = match(rel, CONDENSED_ORDER)) %>% arrange(x, segment) %>%
    group_by(rel) %>% mutate(ymax = cumsum(value), ymin = ymax - value) %>% ungroup() %>%
    mutate(xmin = x - 0.35, xmax = x + 0.35, signif = TRUE, weights = "all")
  tot <- pp %>% group_by(rel, x) %>% summarize(pred = max(ymax), .groups = "drop")
  lab <- po %>% left_join(tot, by = c("rel", "x")) %>% mutate(y = pmax(hi, pred) + 0.02)
  stripes <- seg_stripes(pp, spacing = 0.02, slope = 0.25)
  ggplot() +
    bio_split() +
    geom_rect(data = pp, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = segment),
              key_glyph = draw_key_striped) +
    {if (nrow(stripes)) geom_segment(data = stripes, aes(x = x, xend = xend, y = y, yend = yend),
                                     color = "white", linewidth = 0.35)} +
    geom_hline(yintercept = 0) +   # drawn over the bars
    geom_errorbar(data = po, aes(x = x, ymin = lo, ymax = hi), width = 0.12) +
    geom_point(data = po, aes(x = x, y = corr), size = 1.9) +
    geom_text(data = lab, aes(x = x, y = y, label = label), size = 2.6, vjust = 0, lineheight = 0.85) +
    scale_fill_manual(values = SEG_FILL, labels = SEG_LABELS, breaks = SEG_ORDER, name = "Predicted from") +
    scale_x_continuous(breaks = NULL, expand = expansion(add = 0.6)) +
    scale_y_continuous(breaks = seq(-0.05, 1, by = 0.05), expand = c(0, 0)) +
    coord_cartesian(ylim = c(-0.05, 1)) +
    labs(x = "Relationship Type", y = "Correlation", title = title,
         caption = paste0("Points: inverse-variance pooled correlation of each relationship type (95% CI).\n",
                          "Bars: predicted correlation split by component; A = random-mating part + the excess under assortative mating;\n",
                          "C = beyond-sibling part + within-sibling part.")) +
    theme_bw() +
    theme(axis.ticks.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_line(linetype = "dashed"), panel.grid.minor.y = element_line(linetype = "dashed"),
          legend.position = "inside", legend.position.inside = c(0.985, 0.98), legend.justification = c("right", "top"),
          legend.box.background = element_rect(fill = "white", color = "black"), legend.box.margin = margin(1, 1, 1, 1))
}

# Component correlation decay (no observations): for given parameters, the share of each component
# that relatives of each type have in common, from results/decay/{name}.csv (scripts/07_decay.py),
# scaled so a pair whose kernels are all 1 (MZ twins) has correlation 1. Categories within a type are
# pooled with the observations' inverse-variance weights, as in the condensed fit plots.
plot_component_decay <- function(dec, title = NULL, subtitle = NULL) {
  parts <- c(A0 = "part_A0", A_AM = "part_A_AM", GxG = "comp_GxG", house = "comp_house", fam_ext = "part_fam_ext", fam_sib = "part_fam_sib")
  parts <- parts[parts %in% names(dec)]
  pp <- dec %>% filter(rel %in% CONDENSED_ORDER) %>% mutate(w = 1 / se^2) %>%
    select(rel, w, all_of(unname(parts))) %>%
    pivot_longer(all_of(unname(parts)), names_to = "col", values_to = "value") %>%
    mutate(segment = factor(names(parts)[match(col, parts)], SEG_ORDER)) %>%
    group_by(rel, segment) %>% summarize(value = sum(w * value) / sum(w), .groups = "drop") %>%
    mutate(x = match(rel, CONDENSED_ORDER)) %>% arrange(x, segment) %>%
    group_by(rel) %>% mutate(ymax = cumsum(value), ymin = ymax - value) %>% ungroup() %>%
    mutate(xmin = x - 0.35, xmax = x + 0.35, signif = TRUE, weights = "all")
  lab <- pp %>% group_by(rel, x) %>% summarize(y = max(ymax), .groups = "drop") %>%
    mutate(label = paste0(str_wrap(CONDENSED_LABELS[rel], 10), "\n", sprintf("%.2f", y)))
  stripes <- seg_stripes(pp, spacing = 0.025, slope = 0.25)
  ggplot() +
    bio_split(y = 0.8) +
    geom_rect(data = pp, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = segment),
              key_glyph = draw_key_striped) +
    {if (nrow(stripes)) geom_segment(data = stripes, aes(x = x, xend = xend, y = y, yend = yend),
                                     color = "white", linewidth = 0.35)} +
    geom_hline(yintercept = 0) +
    geom_text(data = lab, aes(x = x, y = y + 0.015, label = label), size = 2.6, vjust = 0, lineheight = 0.85) +
    scale_fill_manual(values = SEG_FILL, labels = SEG_LABELS, breaks = SEG_ORDER, name = "Shared via") +
    scale_x_continuous(breaks = NULL, expand = expansion(add = 0.6)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1), expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, 1.14)) +
    labs(x = "Relationship Type", y = "Correlation (scaled: all kernels 1 \u2192 1)", title = title, subtitle = subtitle,
         caption = paste0("Contributions divided by the summed variance components",
                          if ("part_A0" %in% names(dec)) "; A = random-mating part + the excess under assortative mating." else ".",
                          if ("part_fam_sib" %in% names(dec)) "\nC = beyond-sibling part + within-sibling part." else "")) +
    theme_bw() +
    theme(axis.ticks.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_line(linetype = "dashed"), panel.grid.minor.y = element_blank(),
          legend.position = "bottom")
}

# Parameter sampling correlations from the leave-one-type-out refits (scripts/08; figure S7).
# One pairs plot per core model x weighting. Lower triangle: the 10 single-block refit estimates
# (colour = block left out; blue ramp = biological, orange ramp = non-biological, darker = closer) and
# the full-data estimate (black cross); axes zoomed on the refits. Upper triangle: the jackknife
# correlation r. Diagonal: the estimate and its Marginal Type CI.
PAIR_PARAMS <- c("V_A", "m", "V_GxG", "V_family", "r_f")

pairs_grid_px <- function(k) max(420 * k + 300, 1900)   # side of the k x k grid, in pixels at 300 dpi

# title = NULL drops the heading, foot = FALSE the explanatory footnote; legend = "right" puts the
# key beside the grid (one column) instead of below it.
plot_param_pairs <- function(cvp, est, jk, pc, title = NULL, foot = TRUE, legend = c("bottom", "right"),
                             base = 9) {
  legend <- match.arg(legend)
  ps <- PAIR_PARAMS[PAIR_PARAMS %in% names(est)]
  ps <- ps[!is.na(unlist(est[1, ps]))]
  k <- length(ps)
  e <- setNames(as.numeric(est[1, ps]), ps)
  lo <- setNames(jk$lo, jk$param); hi <- setNames(jk$hi, jk$param)
  # axis range: the refits and the full-data estimate, padded
  lims <- lapply(ps, function(p) {
    r <- range(c(cvp[[p]], e[[p]]), na.rm = TRUE)
    if (diff(r) < 1e-6) r <- r + c(-0.02, 0.02)
    r + c(-1, 1) * 0.12 * diff(r)
  }) %>% setNames(ps)
  th <- theme_bw(base_size = base) +
    theme(panel.grid.minor = element_blank(), plot.margin = margin(4, 8, 2, 2),
          axis.title = element_blank(), legend.position = "none")
  panels <- list()
  for (i in seq_len(k)) for (j in seq_len(k)) {
    pi_ <- ps[i]; pj <- ps[j]
    if (i == j) {
      # symbol (plotmath) above the estimate and its Marginal Type CI
      g <- ggplot() + annotate("text", x = 0.5, y = 0.7, label = par_math(pi_), parse = TRUE, size = base * 0.5) +
        annotate("text", x = 0.5, y = 0.36, size = base * 0.38, lineheight = 1.1,
                 label = sprintf("%.3f\n[%.3f, %.3f]", e[[pi_]], lo[[pi_]], hi[[pi_]])) +
        xlim(0, 1) + ylim(0, 1) + theme_void() +
        theme(panel.background = element_rect(fill = "grey95", color = NA), plot.margin = margin(2, 2, 2, 2))
    } else if (i < j) {
      r <- pc %>% filter((p1 == pi_ & p2 == pj) | (p1 == pj & p2 == pi_))
      # r only: its interval (scripts/08) is too wide to be informative with 10 blocks
      lab <- if (nrow(r) == 0 || is.na(r$r)) "r undefined\n(a parameter does not move)" else sprintf("r = %.2f", r$r)
      rv <- if (nrow(r) == 0 || is.na(r$r)) 0 else r$r
      g <- ggplot() + annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,
                               fill = scales::gradient_n_pal(c("#2a78d6", "#f0efec", "#e34948"), c(-1, 0, 1))(rv), alpha = 0.55) +
        annotate("text", x = 0.5, y = 0.5, label = lab, size = base * 0.48, lineheight = 1.1) +
        scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
        theme_void() + theme(plot.margin = margin(2, 2, 2, 2))
    } else {
      g <- ggplot() +
        geom_point(data = cvp, aes(x = .data[[pj]], y = .data[[pi_]], fill = fold, shape = fold),
                   size = 2.4, stroke = 0.4, color = "white") +
        annotate("point", x = e[[pj]], y = e[[pi_]], shape = 3, size = 3.2, stroke = 1.1) +
        scale_fill_manual(values = FOLD_FILL) + scale_shape_manual(values = FOLD_SHAPE) +
        scale_x_continuous(breaks = scales::breaks_pretty(n = 3)) +
        scale_y_continuous(breaks = scales::breaks_pretty(n = 3)) +
        coord_cartesian(xlim = lims[[pj]], ylim = lims[[pi_]], expand = FALSE) + th +
        theme(axis.text.x = if (i == k) element_text(size = base * 0.78) else element_blank(),
              axis.text.y = if (j == 1) element_text(size = base * 0.78) else element_blank(),
              axis.ticks.x = if (i == k) element_line() else element_blank(),
              axis.ticks.y = if (j == 1) element_line() else element_blank())
    }
    panels[[length(panels) + 1]] <- g
  }
  grid <- cowplot::plot_grid(plotlist = panels, nrow = k, align = "hv", axis = "tblr")
  key <- ggplot(tibble(fold = factor(names(FOLD_FILL), names(FOLD_FILL)), x = seq_along(FOLD_FILL))) +
    geom_point(aes(x = x, y = x, fill = fold, shape = fold), size = 3, color = "white", stroke = 0.4) +
    scale_fill_manual(values = FOLD_FILL, labels = FOLD_LABELS, name = "Refit leaving out") +
    scale_shape_manual(values = FOLD_SHAPE, labels = FOLD_LABELS, name = "Refit leaving out") +
    theme_bw(base_size = if (legend == "right") max(base, 12) else base) + theme(legend.position = legend) +
    guides(fill = if (legend == "bottom") guide_legend(nrow = 3, title.position = "top") else guide_legend(ncol = 1))
  leg <- cowplot::get_plot_component(key, paste0("guide-box-", legend), return_all = TRUE)
  if (legend == "right") {
    return(cowplot::plot_grid(grid, leg, nrow = 1, rel_widths = c(1, 0.32)) +
             theme(plot.background = element_rect(fill = "white", color = NA), plot.margin = margin(4, 4, 4, 4)))
  }
  head <- cowplot::ggdraw() + cowplot::draw_label(title %||% "", x = 0.01, hjust = 0, size = 12)
  ft <- cowplot::ggdraw() + cowplot::draw_label(paste0(
    "Lower: leave-one-type-out refit estimates (axes zoomed on them); black cross = full-data estimate.\n",
    "Upper: jackknife correlation of the two estimates (correlation of the refits).\n",
    "Diagonal: estimate and 95% Marginal Type CI."),
    x = 0.01, hjust = 0, size = 7.5, lineheight = 1.1)
  gpx <- pairs_grid_px(k)
  cowplot::plot_grid(head, grid, leg, if (foot) ft else NULL, ncol = 1,
                     rel_heights = c(if (is.null(title)) 0 else 130, gpx, 330, if (foot) 280 else 0)) +
    theme(plot.background = element_rect(fill = "white", color = NA), plot.margin = margin(4, 10, 4, 14))
}

# ---- M1. Eftedal's estimator -> ours ----------------------------------------------------------
if (has("corrections", "steps.csv")) {
  st <- rd("corrections", "steps.csv")
  lev <- rev(unique(st$label[order(st$order)]))
  PL <- c(V_A = "V[A]", V_A0 = "V[A]^0==V[A]*(1-m)", m = "m~~'(mates\\' genetic correlation)'",
          r = "r~~'(mates\\' phenotypic correlation)'")
  d <- st %>% mutate(label = factor(label, lev), panel = factor(PL[param], PL),
                     lo_c = pmax(lo, 0), hi_c = pmin(hi, 1))
  pts <- d %>% distinct(label, panel, estimate)
  g <- ggplot(d, aes(y = label)) + overlay_ci(d) +
    geom_point(data = pts, aes(x = estimate), size = 2.6) +
    facet_wrap(~panel, scales = "free_x", nrow = 1, labeller = label_parsed) + wrap_y +
    labs(x = "Estimate", y = NULL, title = "Additive model: Eftedal et al.'s estimator vs. the revised estimator",
         caption = "Both CIs drawn from the same point; the Marginal Pair CI is often narrower than the point itself") +
    theme_bw() + theme(legend.position = "bottom")
  save_fig(g, "m1_corrections", width = 3900, height = 1300)
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
  labs <- c(sprintf("Eftedal et al. (Marginal Pair): h²=%.2f", ar$V_A[which(ar$fit_id == ie)]),
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
  il <- rd("corrections", "inlaw.csv") %>% filter(quantity != "difference") %>%
    mutate(label = factor(label, unique(label)),
           cohort = factor(ifelse(quantity == "bio only", "Bio only", "Bio + in-law"), c("Bio only", "Bio + in-law")),
           x = as.numeric(cohort), ci = factor(CI_LABELS[method], names(CI_COLORS)),
           clip_lo = lo < 0, clip_hi = hi > 1, ymin = pmax(lo, 0), ymax = pmin(hi, 1),
           lw = ifelse(method == "parametric bootstrap", 2.4, 0.9))
  arrow_end <- function(d, lo_end) {
    d <- d %>% filter(if (lo_end) clip_lo else clip_hi)
    geom_segment(data = d, aes(x = x, xend = x, y = if (lo_end) ymin + 0.04 else ymax - 0.04,
                               yend = if (lo_end) ymin else ymax, color = ci),
                 linewidth = 0.9, arrow = arrow(length = unit(0.09, "in"), type = "closed"), show.legend = FALSE)
  }
  jk <- il %>% filter(method == "block jackknife"); bt <- il %>% filter(method == "parametric bootstrap")
  g <- ggplot(il, aes(x = x)) +
    geom_line(data = il %>% distinct(label, x, estimate), aes(y = estimate, group = label), color = "grey60") +
    geom_segment(data = jk, aes(xend = x, y = ymin, yend = ymax, color = ci), linewidth = 0.9) +
    geom_segment(data = jk %>% filter(!clip_lo), aes(x = x - 0.04, xend = x + 0.04, y = ymin, yend = ymin, color = ci), linewidth = 0.9) +
    geom_segment(data = jk %>% filter(!clip_hi), aes(x = x - 0.04, xend = x + 0.04, y = ymax, yend = ymax, color = ci), linewidth = 0.9) +
    arrow_end(jk, TRUE) + arrow_end(jk, FALSE) +
    geom_segment(data = bt, aes(xend = x, y = ymin, yend = ymax, color = ci), linewidth = 2.4) +
    geom_point(data = il %>% distinct(label, x, estimate), aes(y = estimate), size = 2.6) +
    facet_grid(~label, labeller = label_wrap_gen(30)) +
    scale_x_continuous(breaks = 1:2, labels = c("Bio only", "Bio + in-law"), expand = expansion(add = 0.4)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1), expand = c(0, 0)) +
    scale_color_manual(values = CI_COLORS, breaks = names(CI_COLORS), name = "95% interval") +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = "h² (additive model)", title = "Does adding in-laws change h²?",
         caption = "Both CIs from the same point; arrowhead: CI runs past [0, 1]. The Marginal Pair CI is often narrower than the point.") +
    theme_bw() + theme(legend.position = "bottom", panel.grid.minor = element_blank())
  save_fig(g, "m2_inlaw", width = 2400, height = 1500)
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
      labs <- sapply(ms, function(mod) model_par_label(mod, fs %>% filter(model == !!mod, weights == !!wt),
                                                       c("V_A", "V_GxG", "V_family", "r_f", "m", "r")))
      pr <- bind_rows(lapply(seq_along(ms), function(i)
        fp %>% filter(model == ms[i], weights == !!wt) %>% mutate(model = labs[i])))
      cols <- setNames(PAL[seq_along(ms)], labs)
      save_fig(plot_fits(obs, pr, colors = cols, labels = parse(text = labs), legend = "bottom") +
                 labs(title = WEIGHT_LABELS[[wt]]), height = 1950,
               sprintf("m3b_eftedalfit_%s_%s", cmp, wt))
      save_fig(plot_fits_condensed(obs, pr, colors = cols, title = WEIGHT_LABELS[[wt]], labels = parse(text = labs)),
               sprintf("m3b_condensedfit_%s_%s", cmp, wt), width = 2600, height = 1500)
    }
  }
}

# ---- M3d. condensed fit-component plots ------------------------------------------------------
if (has("models", "core_predictions.csv")) {
  cpp <- rd("models", "core_predictions.csv")
  if ("part_A0" %in% names(cpp)) {
    cf <- rd("models", "core_fits.csv")
    for (wt in c("hybrid_block")) {   # Marginal Type only (user request)
      for (mod in CORE_ORDER) {
        ttl <- paste0(show_model(mod), " (", WEIGHT_LABELS[[wt]], ")")
        sub <- fmt_par(cf %>% filter(model == !!mod, weights == !!wt), c("V_A", "V_A0", "V_GxG", "V_family", "V_fam_ext", "V_fam_sib", "r_f", "m", "r"))
        save_fig(plot_condensed_components(cpp %>% filter(model == !!mod, weights == !!wt), title = ttl) +
                   labs(subtitle = parse(text = sub)),
                 sprintf("m3d_condensedcomp_%s_%s", gsub("\\+", "", gsub(">", "_to_", mod)), wt), width = 2600, height = 1500)
      }
    }
  }
}

# ---- M3e. where each core model's CV error comes from -------------------------------------------
# Left: the CV score (mean over the 10 held-out blocks of the within-block weighted MSE) stacked by
# block, the adoptive blocks at the bottom; label = CV RMSE and the adoptive blocks' share. Right: CV
# RMSE over all blocks vs without the two adoptive blocks.
if (has("models", "core_cv_scores.csv")) {
  ADOPT <- c("adoptive sibs", "adoptive cousins")
  STACK <- c(ADOPT, setdiff(names(FOLD_FILL), ADOPT))
  sc <- rd("models", "core_cv_scores.csv") %>% filter(fold %in% names(FOLD_FILL)) %>%
    group_by(weights, model) %>% mutate(contrib = score / n() * 1000) %>% ungroup() %>%
    mutate(model = factor(model, CORE_ORDER), weights = factor(weights, c("invvar", "hybrid_block")),
           fold = factor(fold, rev(STACK)))
  tot <- sc %>% group_by(weights, model) %>%
    summarize(total = sum(contrib), adopt = sum(contrib[fold %in% ADOPT]),
              rmse_all = sqrt(mean(score)), rmse_noadopt = sqrt(mean(score[!fold %in% ADOPT])), .groups = "drop") %>%
    mutate(label = sprintf("%.3f\n%.0f%%", rmse_all, 100 * adopt / total))
  g1 <- ggplot(sc, aes(x = model, y = contrib, fill = fold)) +
    geom_col(width = 0.72, color = "white", linewidth = 0.25) +
    geom_text(data = tot, aes(x = model, y = total, label = label), inherit.aes = FALSE,
              vjust = -0.25, size = 2.7, lineheight = 0.9) +
    facet_wrap(~weights, ncol = 1, labeller = as_labeller(WEIGHT_LABELS)) +
    scale_fill_manual(values = FOLD_FILL, labels = FOLD_LABELS, breaks = STACK, name = "Held-out block") +
    scale_x_discrete(labels = show_model) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
    labs(x = NULL, y = "CV score (×10⁻³), stacked by held-out block",
         subtitle = "Label: CV RMSE and the adoptive blocks' share of the score") +
    theme_bw(base_size = 10) +
    theme(panel.grid.major.x = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1),
          legend.position = "bottom", legend.key.spacing.x = unit(14, "pt")) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE))
  db <- tot %>% select(weights, model, rmse_all, rmse_noadopt) %>%
    pivot_longer(c(rmse_all, rmse_noadopt), names_to = "set", values_to = "rmse") %>%
    mutate(set = factor(set, c("rmse_all", "rmse_noadopt"),
                        c("all 10 blocks", "without the adoptive blocks")))
  g2 <- ggplot(db, aes(x = model, y = rmse)) +
    geom_line(aes(group = model), color = "grey60", linewidth = 0.6) +
    geom_point(aes(color = set, shape = set), size = 2.6) +
    facet_wrap(~weights, ncol = 1, labeller = as_labeller(WEIGHT_LABELS)) +
    scale_color_manual(values = c("all 10 blocks" = "#eb6834", "without the adoptive blocks" = "#184f95"), name = NULL) +
    scale_shape_manual(values = c("all 10 blocks" = 16, "without the adoptive blocks" = 17), name = NULL) +
    scale_x_discrete(labels = show_model) +
    labs(x = NULL, y = "CV RMSE (correlation units)", subtitle = "CV RMSE with and without the adoptive blocks") +
    theme_bw(base_size = 10) +
    theme(panel.grid.major.x = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1),
          legend.position = "bottom", legend.direction = "vertical")
  g <- cowplot::plot_grid(g1, g2, nrow = 1, rel_widths = c(1.35, 1), align = "h", axis = "tb")
  g <- cowplot::plot_grid(cowplot::ggdraw() + cowplot::draw_label(
    "Core models: where the cross-validation error comes from", x = 0.01, hjust = 0, size = 13),
    g, ncol = 1, rel_heights = c(0.05, 1)) + theme(plot.background = element_rect(fill = "white", color = NA))
  save_fig(g, "m3e_cv_by_block", width = 3600, height = 2400)
}

# ---- M4. component correlation decay across relationship types (scripts/07_decay.py) -----------
DECAY_TITLES <- c(genetics = "Shared genetic effects across relationship types",
                  environment = "Shared environment (C) across relationship types")
for (f in Sys.glob(file.path(RESULTS, "decay", "*.csv"))) {
  name <- tools::file_path_sans_ext(basename(f))
  dec <- as_tibble(fread(f))
  pars <- dec %>% slice(1) %>% select(any_of(c("V_A", "V_GxG", "V_house", "V_family", "pi", "r_f", "m", "r")))
  pars <- if ("V_A" %in% names(pars)) pars %>% mutate(V_A0 = V_A * (1 - m), .after = V_A) else pars %>% select(-m, -r)
  if (sum(startsWith(names(pars), "V_")) == 1) pars <- pars %>% select(-starts_with("V_"))   # scaled away
  ttl <- if (name %in% names(DECAY_TITLES)) DECAY_TITLES[[name]] else name
  save_fig(plot_component_decay(dec, title = ttl, subtitle = parse(text = model_par_label(dec$model[1], pars, names(pars)))),
           paste0("m4_decay_", name), width = 2600, height = 1500)
}

# ---- S1. cross-validation -------------------------------------------------------------------
for (set in c("core", "joint", "core_am2")) {
  if (!has("models", paste0(set, "_cv_paired.csv"))) next
  order <- if (set == "joint") JOINT_ORDER else CORE_ORDER
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
    labs(x = "CV score minus the best model's (squared-correlation units; paired over the held-out types, ±1.96 SE)",
         y = NULL, title = sprintf("Block cross-validation (%s models)", set)) +
    theme_bw()
  save_fig(g, sprintf("s1_cv_%s", set), width = 3000, height = 300 + 170 * length(order))
  fsc <- rd("models", paste0(set, "_cv_scores.csv")) %>%
    mutate(model = factor(show_model(model), rev(show_model(order))), fold = factor(fold, FOLD_ORDER),
           weights = WEIGHT_LABELS[weights])
  g <- ggplot(fsc, aes(x = fold, y = model, fill = rmse)) + geom_tile() +
    geom_text(aes(label = sprintf("%.3f", rmse)), size = 2.2) +
    facet_wrap(~weights) + scale_fill_viridis_c(option = "magma", direction = -1, name = "RMSE") +
    labs(x = "Held-out type", y = NULL, title = sprintf("Held-out error by type (%s models)", set)) +
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
    labs(x = expression(V[A]), y = "m (mates' genetic correlation)", color = NULL,
         title = expression(paste(V[A], " vs ", m, " across the core models (95% Marginal Type CIs)")),
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
           param = factor(par_math(param), par_math(c("V_A", "V_GxG", "V_house", "V_family", "V_sum", "pi", "r_f", "r"))),
           lo = pmax(lo, 0), hi = pmin(hi, 1.1))
  g <- ggplot(iv, aes(x = estimate, xmin = lo, xmax = hi, y = model, color = method)) +
    geom_pointrange(position = position_dodge(width = 0.6), size = 0.25) +
    facet_grid(weights ~ param, scales = "free_x", labeller = labeller(param = label_parsed)) +
    labs(x = "Estimate and 95% interval", y = NULL, color = NULL) +
    theme_bw() + theme(legend.position = "bottom")
  save_fig(g, "s3_bootstrap_vs_jackknife", width = 3800, height = 1900)
}
if (has("models", "core_cv_params.csv")) {
  jk <- rd("models", "core_cv_params.csv") %>% filter(fold %in% FOLD_ORDER[1:10]) %>%
    select(model, weights, fold, V_A, V_sum) %>% pivot_longer(c(V_A, V_sum), names_to = "param") %>%
    mutate(model = factor(show_model(model), show_model(CORE_ORDER)), fold = factor(fold, rev(FOLD_ORDER)),
           weights = WEIGHT_LABELS[weights], param = factor(par_math(param), par_math(c("V_A", "V_sum"))))
  g <- ggplot(jk, aes(x = value, y = fold, color = model)) + geom_point(size = 2.2) +
    scale_y_discrete(labels = FOLD_LABELS) +
    facet_grid(param ~ weights, labeller = labeller(param = label_parsed)) + labs(x = "Estimate with the type left out", y = "Type left out", color = NULL) +
    theme_bw()
  save_fig(g, "s3b_jackknife_by_block", width = 3000, height = 2000)
}

# ---- S4. profile over V_A ----------------------------------------------------------------------
if (has("profile", "profile_VA.csv")) {
  pf <- rd("profile", "profile_VA.csv") %>% mutate(weights = WEIGHT_LABELS[weights], model = show_model(model))
  grid <- pf %>% filter(point == "grid")
  marks <- pf %>% filter(point != "grid")
  g1 <- ggplot(grid, aes(x = V_A, y = qlr)) +
    geom_hline(yintercept = qchisq(0.95, 1), linetype = "dashed", color = "grey40") +
    geom_line(linewidth = 0.9) +
    geom_point(data = marks %>% mutate(point = show_model(point)), aes(color = point), size = 3) +
    facet_grid(weights ~ model) + coord_cartesian(ylim = c(0, 60)) +
    scale_color_brewer(palette = "Set1", name = NULL) +
    labs(x = expression(V[A]~"(held fixed; everything else refitted)"),
         y = "Loss increase / residual variance\n(quasi-likelihood ratio; dashed: 3.84)",
         title = expression(paste("How well is ", V[A], " pinned down? Profile over ", V[A]))) +
    theme_bw() + theme(legend.position = "bottom")
  save_fig(g1, "s4_profile_VA_loss", width = 3300, height = 1700)
  comp <- grid %>% select(model, weights, V_A, V_GxG, V_family) %>%
    pivot_longer(c(V_GxG, V_family), names_to = "component") %>% filter(!is.na(value)) %>%
    mutate(component = sub("V_", "", component))
  g2 <- ggplot(comp, aes(x = V_A, y = value, color = component)) + geom_line(linewidth = 0.9) +
    geom_abline(intercept = 1, slope = -1, linetype = "dotted", color = "grey50") +
    facet_grid(weights ~ model) +
    scale_color_manual(values = c(GxG = COMP_FILL[["GxG"]], family = COMP_FILL[["family"]]),
                       breaks = c("GxG", "family"), labels = expression(V[GxG], V[C])) +
    labs(x = expression(V[A]~"(held fixed)"), y = "Refitted variance component", color = NULL,
         title = expression(paste("What the other components do as ", V[A], " is moved (dotted: components sum to 1)"))) +
    theme_bw() + theme(legend.position = "bottom")
  save_fig(g2, "s4b_profile_VA_components", width = 3300, height = 1700)
}

# ---- S6. GxG kernel a^2 (square of the random-mating kernel) vs the AM kernel squared ----------
if (has("models", "core_am2_jackknife.csv") && has("models", "core_jackknife.csv")) {
  save_fig(stacked_components(rd("models", "core_am2_jackknife.csv"), CORE_ORDER,
                              "Core models with the earlier GxG kernel, (AM kernel)\u00b2: variance components"),
           "s6_components_core_am2", width = 2700, height = 2300)
  gx <- CORE_ORDER[grepl("GxG", CORE_ORDER)]
  kk <- bind_rows(rd("models", "core_jackknife.csv") %>% mutate(kernel = "a\u00b2 (random-mating kernel squared; main)"),
                  rd("models", "core_am2_jackknife.csv") %>% mutate(kernel = "(AM kernel)\u00b2 (earlier)")) %>%
    filter(model %in% gx, param %in% c("V_A", "V_GxG", "V_family", "V_sum")) %>%
    mutate(model = factor(show_model(model), rev(show_model(gx))),
           param = factor(par_math(param), par_math(c("V_A", "V_GxG", "V_family", "V_sum"))),
           weights = factor(WEIGHT_LABELS[weights], WEIGHT_LABELS[c("invvar", "hybrid_block")]))
  g <- ggplot(kk, aes(x = estimate, xmin = pmax(lo, 0), xmax = pmin(hi, 1), y = model, color = kernel)) +
    geom_pointrange(position = position_dodge(width = 0.55), size = 0.3) +
    facet_grid(weights ~ param, labeller = labeller(param = label_parsed)) +
    labs(x = "Estimate (95% Marginal Type CI)", y = NULL, color = "GxG kernel",
         title = "GxG kernel: a\u00b2 (main) vs the square of the assortative-mating kernel") +
    theme_bw() + theme(legend.position = "bottom")
  save_fig(g, "s6b_gxg_kernel_comparison", width = 3200, height = 1500)
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
# ---- S7. parameter sampling correlations from the leave-one-type-out refits (scripts/08) ------------
# (plot_param_pairs is defined with the helpers above)
if (has("models", "core_param_corr.csv")) {
  cvp_all <- rd("models", "core_cv_params.csv") %>% filter(fold %in% names(FOLD_FILL)) %>%
    mutate(fold = factor(fold, names(FOLD_FILL)))
  fits_all <- rd("models", "core_fits.csv"); jk_all <- rd("models", "core_jackknife.csv")
  pc_all <- rd("models", "core_param_corr.csv")
  for (wt in c("hybrid_block", "invvar")) for (mod in CORE_ORDER) {
    est <- fits_all %>% filter(model == !!mod, weights == !!wt)
    k <- sum(!is.na(unlist(est[1, intersect(PAIR_PARAMS, names(est))])))
    g <- plot_param_pairs(cvp_all %>% filter(model == !!mod, weights == !!wt), est,
                          jk_all %>% filter(model == !!mod, weights == !!wt),
                          pc_all %>% filter(model == !!mod, weights == !!wt),
                          paste0(show_model(mod), ": parameter sampling correlations (", WEIGHT_LABELS[[wt]], ")"))
    gpx <- pairs_grid_px(k)
    save_fig(g, sprintf("s7_paramcorr_%s_%s", gsub("\\+", "", gsub(">", "_to_", mod)), wt),
             width = gpx, height = gpx + 130 + 330 + 280)
  }
}

cat("figures written to", FIGURES, "\n")
