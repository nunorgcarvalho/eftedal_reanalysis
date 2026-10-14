# Plotting helpers. All inputs are the tidy CSVs written by the Python scripts.
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(ggnewscale)
  library(ggrepel)
  library(scales)
})

# REPO, RESULTS and FIGURES are set by the calling script (R/make_figures.R or
# R/exploratory_figures.R).

z95 <- qnorm(0.975)

side_colors <- c(maternal = "red", other = "#c7b267", paternal = "blue")
sex_shapes <- c(os = "triangle", ssf = "circle", ssm = "square")

# ---- observed correlations, laid out as in Eftedal et al. Fig. 2 ------------------------

# x position: categories in file order, with a one-unit gap between relationship types
layout_obs <- function(obs) {
  types <- unique(obs$rel)
  obs %>%
    arrange(row) %>%
    mutate(corr_low = corr - z95 * se, corr_upp = corr + z95 * se,
           type_index = match(rel, types),
           x = row + type_index)
}

gg_obs <- function(obs) {
  ggplot(obs, aes(x = x, y = corr)) +
    geom_hline(yintercept = 0) +
    geom_point(aes(color = side, shape = sex), size = 3) +
    geom_errorbar(aes(ymin = corr_low, ymax = corr_upp), width = 0.6) +
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = seq(-0.05, 1, by = 0.05), expand = c(0, 0)) +
    scale_color_manual(values = side_colors, labels = c("Maternal", "Other", "Paternal"),
                       name = "Sub-category of relative") +
    scale_shape_manual(values = sex_shapes,
                       labels = c("Opposite sex", "Same-sex female", "Same-sex male"),
                       name = "Gender composition") +
    coord_cartesian(ylim = c(-0.05, 1)) +
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_line(linetype = "dashed"),
          panel.grid.minor.y = element_line(linetype = "dashed"),
          legend.position = "inside", legend.position.inside = c(0.985, 0.98),
          legend.justification = c("right", "top"), legend.box = "horizontal",
          legend.box.background = element_rect(fill = "white", color = "black"),
          legend.box.margin = margin(1, 1, 1, 1)) +
    labs(y = "Correlation")
}

# One horizontal bar per relationship type (split by side when the prediction differs by
# side, e.g. half sibs under a household term). `preds`: rows of a *_predictions.csv for
# the fits to show, with a `model` column used for the legend.
pred_segments <- function(obs, preds) {
  preds %>%
    inner_join(obs %>% select(row, x), by = "row") %>%
    group_by(model, rel) %>%
    mutate(split = n_distinct(round(pred, 9)) > 1) %>%
    group_by(model, rel, grp = ifelse(split, side, "all")) %>%
    summarize(pred = first(pred), x_min = min(x), x_max = max(x),
              in_cohort = any(in_cohort), .groups = "drop")
}

# legend = "auto": the model legend sits inside the panel (top right, as in Eftedal Fig. 2)
# unless its labels are long enough to cover the MZ points; then it moves below the panel
# (the relative-category legends stay inside).
connect_path <- function(seg) {
  seg %>% arrange(model, x_min) %>%
    reframe(x = c(x_min - 0.5, x_max + 0.5), pred = c(pred, pred), .by = c(model, rel, grp, x_min)) %>%
    arrange(model, x_min, x)
}

plot_fits <- function(obs, preds, colors = NULL, R2 = NULL, legend = "auto") {
  obs <- layout_obs(obs)
  seg <- pred_segments(obs, preds)
  labs_df <- obs %>% group_by(rel) %>%
    summarize(x = (min(x) + max(x)) / 2, top = max(corr_upp), .groups = "drop") %>%
    left_join(seg %>% group_by(rel) %>% summarize(pmax = max(pred), .groups = "drop"), by = "rel") %>%
    mutate(y = pmax(top, pmax) + 0.025, label = str_wrap(rel, 15))
  if (is.null(colors)) colors <- setNames(hue_pal()(n_distinct(seg$model)), unique(seg$model))
  if (legend == "auto") legend <- if (max(nchar(unique(seg$model))) > 25) "bottom" else "inside"
  model_guide <- if (legend == "bottom") guide_legend(position = "bottom", ncol = 1) else "legend"
  gg <- gg_obs(obs) +
    new_scale_color() +
    # thin line joining each model's consecutive bars, so departures between models are visible
    geom_path(data = connect_path(seg), inherit.aes = FALSE, linewidth = 0.35, alpha = 0.8,
              aes(x = x, y = pred, color = model, group = model)) +
    geom_segment(data = seg, inherit.aes = FALSE, linewidth = 1.5, alpha = 0.6,
                 aes(x = x_min - 0.5, xend = x_max + 0.5, y = pred, yend = pred, color = model,
                     linetype = in_cohort)) +
    scale_color_manual(values = colors, breaks = names(colors), name = "Model", guide = model_guide) +
    scale_linetype_manual(values = c(`TRUE` = "solid", `FALSE` = "22"), guide = "none") +
    geom_text_repel(data = labs_df, inherit.aes = FALSE, size = 3, direction = "y", vjust = 1,
                    min.segment.length = 999, aes(x = x, y = y, label = label))
  if (!is.null(R2)) {
    gg <- gg + annotate("text", x = 0, y = 0.05, hjust = 0, vjust = 0, size = 5,
                        label = R2)
  }
  gg
}

save_fig <- function(gg, name, width = 3000, height = 1600) {
  # PNG for notebooks; vector PDF (cairo, for unicode glyphs) for the poster
  dir.create(FIGURES, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(FIGURES, paste0(name, ".png")), gg, width = width, height = height,
         units = "px", dpi = 300)
  ggsave(file.path(FIGURES, paste0(name, ".pdf")), gg, width = width / 300, height = height / 300,
         units = "in", device = cairo_pdf)
  invisible(gg)
}


# ---- condensed fit plot: one pooled point per relationship type --------------------------------
# Categories of a relationship type (maternal/paternal/other x gender composition) are pooled by a
# fixed-effect (inverse-variance) meta-analysis: rho = sum(w rho) / sum(w), SE = 1 / sqrt(sum(w)),
# w = 1/SE^2. Predictions are pooled with the same weights. UZ twins and relatives through UZ twins
# are omitted. A dashed line separates the biological relatives from the in-laws (step-sibs,
# co-cousins) and adoptees.
CONDENSED_ORDER <- c("mz", "dz", "sibs", "half sibs", "1st cousins", "half cousins", "2nd cousins",
                     "step sibs", "cocousins", "adoptive sibs", "adoptive 1st cousins",
                     "adoptive 2nd cousins", "adoptive cocousins")
CONDENSED_LABELS <- c(mz = "MZ twins", uz = "UZ twins", dz = "DZ twins", sibs = "full-sibs",
                      `half sibs` = "half-sibs", `1st cousins` = "1st cousins", `half cousins` = "half cousins",
                      `2nd cousins` = "2nd cousins", cocousins = "co-cousins", `step sibs` = "step-sibs",
                      `adoptive sibs` = "adoptive sibs", `adoptive 1st cousins` = "adoptive 1st cousins",
                      `adoptive 2nd cousins` = "adoptive 2nd cousins", `adoptive cocousins` = "adoptive co-cousins")

CONDENSED_SPLIT <- match("step sibs", CONDENSED_ORDER) - 0.5   # x of the bio | in-law dashed line

# the dashed line, with small faint vertical labels on either side, centred at height y
bio_split <- function(y = 0.4, size = 2.3) {
  list(geom_vline(xintercept = CONDENSED_SPLIT, linetype = "dashed", color = "grey55"),
       annotate("text", x = CONDENSED_SPLIT - 0.09, y = y, label = "biological relatives", angle = 90,
                vjust = 0, size = size, color = "grey55"),
       annotate("text", x = CONDENSED_SPLIT + 0.09, y = y, label = "non-biological relatives", angle = 90,
                vjust = 1, size = size, color = "grey55"))
}

pool_by_type <- function(obs) {
  obs %>% filter(rel %in% CONDENSED_ORDER) %>% mutate(w = 1 / se^2) %>%
    group_by(rel) %>% summarize(corr = sum(w * corr) / sum(w), se = sqrt(1 / sum(w)), .groups = "drop") %>%
    mutate(x = match(rel, CONDENSED_ORDER))
}

plot_fits_condensed <- function(obs, preds, colors, title = NULL) {
  po <- pool_by_type(obs) %>%
    mutate(lo = corr - z95 * se, hi = corr + z95 * se, label = str_wrap(CONDENSED_LABELS[rel], 10))
  w <- obs %>% transmute(row, w = 1 / se^2)
  pp <- preds %>% filter(rel %in% CONDENSED_ORDER) %>% inner_join(w, by = "row") %>%
    group_by(model, rel) %>% summarize(pred = sum(w * pred) / sum(w), .groups = "drop") %>%
    mutate(x = match(rel, CONDENSED_ORDER))
  path <- pp %>% reframe(x = c(x - 0.35, x + 0.35), pred = c(pred, pred), .by = c(model, rel)) %>%
    arrange(model, x)
  lab <- po %>% left_join(pp %>% group_by(rel) %>% summarize(pmax = max(pred)), by = "rel") %>%
    mutate(y = pmax(hi, pmax) + 0.02)
  ggplot() +
    geom_hline(yintercept = 0) +
    bio_split() +
    geom_path(data = path, aes(x = x, y = pred, color = model, group = model), linewidth = 0.3, alpha = 0.8) +
    geom_segment(data = pp, aes(x = x - 0.35, xend = x + 0.35, y = pred, yend = pred, color = model),
                 linewidth = 0.9, alpha = 0.8) +
    geom_errorbar(data = po, aes(x = x, ymin = lo, ymax = hi), width = 0.12) +
    geom_point(data = po, aes(x = x, y = corr), size = 1.9) +
    geom_text(data = lab, aes(x = x, y = y, label = label), size = 2.6, vjust = 0, lineheight = 0.85) +
    scale_x_continuous(breaks = NULL, expand = expansion(add = 0.6)) +
    scale_y_continuous(breaks = seq(-0.05, 1, by = 0.05), expand = c(0, 0)) +
    scale_color_manual(values = colors, breaks = names(colors), name = "Model") +
    coord_cartesian(ylim = c(-0.05, 1)) +
    labs(x = "Relationship Type", y = "Correlation", title = title,
         caption = "Points: inverse-variance pooled correlation of each relationship type (95% CI); bars: model predictions pooled with the same weights") +
    theme_bw() +
    theme(axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_line(linetype = "dashed"),
          panel.grid.minor.y = element_line(linetype = "dashed"),
          legend.position = "inside", legend.position.inside = c(0.985, 0.98),
          legend.justification = c("right", "top"),
          legend.box.background = element_rect(fill = "white", color = "black"),
          legend.box.margin = margin(1, 1, 1, 1))
}
