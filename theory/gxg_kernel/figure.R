# Figure for gxg_kernel.tex: simulated GxG correlation between relatives vs the candidate kernels.
#   Rscript figure.R   (reads results/comparison.csv, writes results/kernel_check.pdf/.png)
suppressPackageStartupMessages({ library(tidyverse); library(data.table) })
args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
d <- as_tibble(fread(file.path(here, "results", "comparison.csv")))

ORDER <- c("parent-offspring", "full sibs", "half sibs", "grandparent", "avuncular",
           "first cousins", "step sibs", "co-cousins")
lab_regime <- c(strong = "strong AM (ρᵧ ≈ 0.86)", moderate = "moderate AM (ρᵧ ≈ 0.64)")
lab_arch <- c(random = "random pairwise epistasis", aligned = "scale epistasis (G ∝ A²)")

long <- d %>% filter(cls %in% ORDER, regime %in% c("strong", "moderate")) %>%
  transmute(regime, arch, cls, sim = K_G_sim, se = K_G_sim_se,
            `path model (this note)` = ifelse(arch == "random", K_G_path, NA),
            `K² (current kernel)` = K_A_sq, `a² (IBD only)` = a_sq) %>%
  pivot_longer(c(`path model (this note)`, `K² (current kernel)`, `a² (IBD only)`),
               names_to = "prediction", values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(cls = factor(cls, rev(ORDER)),
         regime = factor(lab_regime[regime], lab_regime), arch = factor(lab_arch[arch], lab_arch))
sims <- long %>% distinct(regime, arch, cls, sim, se)

PAL <- c(`path model (this note)` = "#0072B2", `K² (current kernel)` = "#D55E00", `a² (IBD only)` = "#009E73")
SHP <- c(`path model (this note)` = 17, `K² (current kernel)` = 15, `a² (IBD only)` = 4)
g <- ggplot(long, aes(y = cls)) +
  geom_errorbarh(data = sims, aes(xmin = sim - 1.96 * se, xmax = sim + 1.96 * se), height = 0, linewidth = 0.9,
                 color = "grey25") +
  geom_point(data = sims, aes(x = sim), shape = 21, fill = "white", color = "grey15", size = 3.2, stroke = 1) +
  geom_point(aes(x = value, color = prediction, shape = prediction), size = 2.3,
             position = position_nudge(y = 0.22)) +
  facet_grid(arch ~ regime) +
  scale_color_manual(values = PAL, name = NULL) + scale_shape_manual(values = SHP, name = NULL) +
  labs(x = "Correlation of GxG values between relatives  (open circles: simulation, 95% CI over 3 replicates)",
       y = NULL) +
  theme_bw(base_size = 10) + theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(file.path(here, "results", "kernel_check.pdf"), g, width = 9, height = 6.2, device = cairo_pdf)
ggsave(file.path(here, "results", "kernel_check.png"), g, width = 9, height = 6.2, dpi = 200)
