# Figure for the linkage section of gxg_kernel.tex: simulated GxG correlation between relatives with
# linked interacting loci vs predictions, by effect architecture (rows) and mating regime (columns).
#   Rscript figure_linked.R   (reads results/linked/comparison.csv, writes results/linked/kernel_linked.pdf/.png)
suppressPackageStartupMessages({ library(tidyverse); library(data.table) })
args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
d <- as_tibble(fread(file.path(here, "results", "linked", "comparison.csv")))

ORDER <- c("parent-offspring", "full sibs", "half sibs", "grandparent", "avuncular",
           "first cousins", "step sibs", "co-cousins")
lab_regime <- c(rm = "random mating", strong = "strong AM (ρᵧ ≈ 0.86)")
lab_arch <- c(uniform = "(1) uniform\nγ_kl",
              beta = "(2) β-scaled\nβ_k β_l δ_kl",
              beta_1mc = "(3a) × (1−c)\nβ_k β_l δ_kl (1−c)",
              `beta_prox_0.05` = "(3b) λ = 5 cM\nβ_k β_l δ_kl e^(−d/λ)",
              `beta_prox_0.001` = "(3c) λ = 0.1 cM\nβ_k β_l δ_kl e^(−d/λ)",
              `beta_prox_dir_0.001` = "(3d) λ = 0.1 cM,\ndirectional\nβ_k β_l e^(−d/λ)")

long <- d %>% filter(cls %in% ORDER, arch %in% names(lab_arch)) %>%
  transmute(regime, arch, cls, sim = K_G, se = K_G_se,
            `linkage path model (this note)` = linked_am_T, `K² (current kernel)` = K_A_sq,
            `a² (unlinked)` = a2) %>%
  pivot_longer(c(`linkage path model (this note)`, `K² (current kernel)`, `a² (unlinked)`),
               names_to = "prediction", values_to = "value") %>%
  mutate(cls = factor(cls, rev(ORDER)),
         regime = factor(lab_regime[regime], lab_regime), arch = factor(lab_arch[arch], lab_arch))
sims <- long %>% distinct(regime, arch, cls, sim, se)

PAL <- c(`linkage path model (this note)` = "#0072B2", `K² (current kernel)` = "#D55E00",
         `a² (unlinked)` = "#009E73")
SHP <- c(`linkage path model (this note)` = 17, `K² (current kernel)` = 15, `a² (unlinked)` = 4)
g <- ggplot(long, aes(y = cls)) +
  geom_errorbar(data = sims, aes(xmin = sim - 1.96 * se, xmax = sim + 1.96 * se), width = 0, orientation = "y",
                 linewidth = 0.9, color = "grey25") +
  geom_point(data = sims, aes(x = sim), shape = 21, fill = "white", color = "grey15", size = 2.8, stroke = 1) +
  geom_point(aes(x = value, color = prediction, shape = prediction), size = 2.1,
             position = position_nudge(y = 0.25)) +
  facet_grid(arch ~ regime) +
  scale_color_manual(values = PAL, name = NULL) + scale_shape_manual(values = SHP, name = NULL) +
  labs(x = "Correlation of GxG values between relatives  (open circles: simulation, 95% CI over 3 replicates)",
       y = NULL) +
  theme_bw(base_size = 9.5) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        strip.text.y = element_text(angle = 0, hjust = 0))
ggsave(file.path(here, "results", "linked", "kernel_linked.pdf"), g, width = 8.5, height = 10.5, device = cairo_pdf)
ggsave(file.path(here, "results", "linked", "kernel_linked.png"), g, width = 8.5, height = 10.5, dpi = 200)
