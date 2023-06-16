################################################################################
#                                  libraries                                   #
################################################################################

library(OptM)
library(tidyverse)
library(patchwork)

################################################################################
#                                  plot optm                                   #
################################################################################

# set working directory
setwd("~/documents/rcoimbra_phd/project_kenya/results/08_gene_flow/admixture_graphs/treemix_runs/")

# read treemix output files
optm <- optM(".", method = "Evanno", tsv = "optm.tsv")

# plot mean and standard deviation for composite likelihoods
p1 <- ggplot(optm, aes(x = m, group = m)) +
  geom_point(aes(y = `mean(Lm)`), alpha = 0.5, size = 2.5) +
  geom_errorbar(aes(ymin = `mean(Lm)` - `sd(Lm)`, ymax = `mean(Lm)` + `sd(Lm)`), alpha = 0.5, width = 0.5) +
  labs(
    x = expression(paste(italic("m"), "(migration edges)")),
    y = "Mean L(m) \u00B1 SD"
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank()
  )

# plot mean and standard deviation for explained variance
p2 <- ggplot(optm, aes(x = m, group = m)) +
  geom_point(aes(y = `mean(f)` * 100), color = "red", alpha = 0.5, size = 2.5) +
  geom_errorbar(aes(ymin = (`mean(f)` - `sd(f)`) * 100, ymax = (`mean(f)` + `sd(f)`) * 100),color = "red", alpha = 0.5, width = 0.5) +
  #geom_hline(aes(yintercept = 99.8), linewidth = 0.5,  linetype = 2) +
  labs(
    x = expression(paste(italic("m"), " (migration edges)")),
    y = "Variance explained (%) \u00B1 SD"
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank()
  )

# plot delta m
p3 <- ggplot(optm, aes(x = m, y = Deltam)) +
  geom_point(color = "blue", alpha = 0.5, size = 2.5) +
  geom_line(color = "blue", alpha = 0.5, linewidth = 0.5) +
  scale_x_continuous(breaks = seq(0, 5)) +
  labs(
    x = expression(paste(italic("m"), " (migration edges)")),
    y = expression(paste(Delta, italic("m")))
  ) +
  theme_minimal()

################################################################################
#                                    figure                                    #
################################################################################

# arrange figure layout
(p1 / p2 / p3) + plot_annotation(tag_levels = "a")

# save plot in '.png' format
ggsave(
  filename = "figureS5.png",
  path     = "~/documents/rcoimbra_phd/project_kenya/figures/revised/",
  width    = 170,
  height   = 170,
  units    = "mm",
  dpi      = 300
)