# Figures for x_kernel.tex.
#   Rscript figure.R   (reads x_kernel_grid.csv and signatures.csv, writes fig/*.pdf/.png)
suppressPackageStartupMessages({ library(tidyverse); library(data.table) })
args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
dir.create(file.path(here, "fig"), showWarnings = FALSE)
save2 <- function(g, name, w, h) {
  ggsave(file.path(here, "fig", paste0(name, ".pdf")), g, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(here, "fig", paste0(name, ".png")), g, width = w, height = h, dpi = 200)
}

# ---- 1: the X correction D(m) = K^X(m) - K^auto(m) ----------------------------------------------
TYPES <- c(sibs = "full sibs", `half sibs` = "half sibs", `1st cousins` = "1st cousins",
           `2nd cousins` = "2nd cousins", `step sibs` = "step sibs", cocousins = "co-cousins")
grid <- as_tibble(fread(file.path(here, "x_kernel_grid.csv"))) %>%
  filter(rel %in% names(TYPES)) %>%
  mutate(type = factor(TYPES[rel], TYPES),
         side = factor(side, c("maternal", "paternal", "other")),
         sex = factor(sex, c("ssf", "ssm", "os")))
PAL <- c(maternal = "#C0392B", paternal = "#2471A3", other = "grey35")
g1 <- ggplot(grid, aes(m, D, color = side, linetype = sex, group = interaction(side, sex))) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.4) +
  geom_vline(xintercept = 0.43, color = "grey75", linewidth = 0.4, linetype = "dashed") +
  geom_line(linewidth = 0.8) +
  facet_wrap(~type, nrow = 2) +
  scale_color_manual(values = PAL, name = "side") +
  scale_linetype_manual(values = c(ssf = "solid", ssm = "22", os = "11"), name = "pair") +
  labs(x = "genotypic mate correlation m  (dashed line: 0.43, Eftedal et al. estimate)",
       y = expression("X correction  " * D == K^X - K^auto)) +
  theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), legend.position = "bottom")
save2(g1, "D_vs_m", 8, 5.2)

# ---- 2: predicted X-only signatures vs observed ------------------------------------------------
sg <- as_tibble(fread(file.path(here, "signatures.csv"))) %>%
  filter(!rel %in% c("dz", "step sibs")) %>%
  mutate(contrast = recode(contrast, `ss - os` = "same-sex minus opposite-sex",
                           `maternal - paternal` = "maternal minus paternal"),
         rel = recode(rel, cocousins = "co-cousins"),
         rel = factor(rel, rev(unique(rel))))
g2 <- ggplot(sg, aes(y = rel)) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.4) +
  geom_errorbar(aes(xmin = 0, xmax = `pred_VA0.7_s0.05`), width = 0, linewidth = 2.6,
                color = "#AED6F1", orientation = "y") +
  geom_errorbar(aes(xmin = `pred_VA0.4_s0.02`, xmax = `pred_VA0.7_s0.02`), width = 0, linewidth = 2.6,
                color = "#2471A3", orientation = "y") +
  geom_errorbar(aes(xmin = observed - 1.96 * obs_se, xmax = observed + 1.96 * obs_se), width = 0,
                linewidth = 0.6, color = "grey15", orientation = "y", position = position_nudge(y = 0.28)) +
  geom_point(aes(x = observed), shape = 21, fill = "white", color = "grey15", size = 2.4, stroke = 0.9,
             position = position_nudge(y = 0.28)) +
  facet_wrap(~contrast, ncol = 2, scales = "free_x") +
  labs(x = "difference in correlation", y = NULL,
       caption = paste("open circles: observed (95% CI);  dark bar: X only, s_X = 0.02, V_A 0.4-0.7;",
                       "light bar: up to s_X = 0.05, V_A = 0.7;  m = 0.43")) +
  theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), plot.caption = element_text(hjust = 0))
save2(g2, "signatures", 8, 3.6)
