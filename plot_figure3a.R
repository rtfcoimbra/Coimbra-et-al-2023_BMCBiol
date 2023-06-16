source("https://bitbucket.org/nygcresearch/treemix/raw/f38bfada3286027a09924d630efa3ad190bda380/src/plotting_funcs.R")

setwd("~/documents/rcoimbra_phd/project_kenya/results/08_gene_flow/admixture_graphs/orientagraph_runs/")

png("figure3.png", width = 210, height = 297, units = "mm", res = 300)

par(mar = c(4, 3, 0, 0))

layout(matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2, byrow = T))

plot_tree("treemix.m2.r9", xmin = 0.3, arrow = 0.075, cex = 0.8, lwd = 1.5, flip = c(7, 8, 10, 15, 16, 20, 27, 28))

invisible(dev.off())