source("https://bitbucket.org/nygcresearch/treemix/raw/f38bfada3286027a09924d630efa3ad190bda380/src/plotting_funcs.R")

setwd("~/documents/rcoimbra_phd/project_kenya/results/08_gene_flow/admixture_graphs/orientagraph_runs/")

png(
  filename = "~/documents/rcoimbra_phd/project_kenya/figures/revised/figureS6.png",
  width    = 210,
  height   = 210,
  units    = "mm",
  res      = 300
)

par(mar = c(5,1,4,1) + 0.1)

layout(matrix(c(rep(c(1,2), 6), rep(c(3,4), 4)), nrow = 2, ncol = 10, byrow = F))

plot_tree("treemix.m1.r2", plus = 0.08, arrow = 0.075, cex = 0.8, flip = c(8, 14, 16, 22, 25, 28))
title(expression(paste(italic("m"), " = 1")))

plot_tree("treemix.m2.r9", plus = 0.08, arrow = 0.075, cex = 0.8, flip = c(7, 8, 10, 15, 16, 20, 27, 28))
title(expression(paste(italic("m"), " = 2")))

plot_resid(stem = "treemix.m1.r2", pop_order = "sample_sets.txt", cex = 0.45, min = -0.00144, max = 0.00144, usemax = F)
title(expression(paste("Residuals from ", italic("m"), " = 1")))

plot_resid(stem = "treemix.m2.r9", pop_order = "sample_sets.txt", cex = 0.45, min = -0.001425, max = 0.001425, usemax = F)
title(expression(paste("Residuals from ", italic("m"), " = 2")))

invisible(dev.off())