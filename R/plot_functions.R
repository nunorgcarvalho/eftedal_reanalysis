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

