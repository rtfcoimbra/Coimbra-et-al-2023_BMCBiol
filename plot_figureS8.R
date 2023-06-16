# Script modified from: calculateDeviance.R; Patrick G. Meirmans (2013) Non-convergence in Bayesian estimation of migration rates. Mol Ecol Res; Supplementary Material

# This script will calculate the Bayesian Deviance from the output of BayesAss (Wilson & Rannala 2003)
# For more information on the use of the deviance, see Faubet et al. (2007)

# To use this script, run BayesAss version 3 with the -t flag to produce a trace-file
# Then set the working directory of R to the folder with the output from BayesAss

library(tidyverse)

# set working directory
setwd("~/documents/rcoimbra_phd/project_kenya/results/08_gene_flow/bayesass/")

# find input files
files <- list.files(pattern = "*.txt")

# burn-in used for the MCMC
burnin <- 2000000

# sampling interval used for the MCMC
sampling.interval <- 2000

# set factor levels
#dataset.lvls <- c("dataset_1", "dataset_2", "dataset_3")
run.lvls <- c("Run 1", "Run 2", "Run 3")

# create a tibble with the data from the trace files
trace <- files %>%
  set_names(str_remove(files, ".trace.txt")) %>%
  map_dfr(read_table, .id = "Run", col_names = TRUE) %>%
  select(Run, State, LogProb) %>%
  #separate(Dataset, c("Dataset", "Run"), "([.])") %>%
  mutate(
    #Dataset = fct_relevel(Dataset, dataset.lvls),
    Run = str_replace_all(Run, c("snps.sampled." = "", "run1" = "Run 1", "run2" = "Run 2", "run3" = "Run 3")),
    Run = fct_relevel(Run, run.lvls)
  )

# create basic plot
p <- ggplot(data = trace, mapping = aes(x = State, y = LogProb)) +
  # draw lineplot of log probabilities
  geom_line(linewidth = 0.25) +
  # draw a vertical line to indicate the end of the burnin
  geom_vline(aes(xintercept = burnin), color = "red", lty = 2) +
  # add matrix of panels defined by two column faceting variables
  facet_wrap(~Run) +
  #facet_grid(cols = vars(Run),
  #           rows = vars(Dataset)) +
  # rename axes' labels
  labs(
    x = "State",
    y = "Log probability"
  ) +
  # change title and axis font size and remove legend
  theme_bw() +
  theme(
    axis.title      = element_text(size = 10),
    axis.text.x     = element_text(size = 8, angle = 45, hjust = 1),
    axis.text.y     = element_text(size = 8),
    strip.text      = element_text(size = 10),
    legend.position = "none"
  )

# calculate the Bayesian deviance
bayes.dev <- trace %>%
  #group_by(Dataset, Run) %>%
  group_by(Run) %>%
  filter(State > burnin & State %% sampling.interval == 0) %>%
  summarise(D = -2 * mean(LogProb))

# add Bayesian deviance label
p + geom_text(
  data    = bayes.dev,
  mapping = aes(label = sprintf("D = %.1f", D)),
  x       = 22000000,
  y       = -200500,
  vjust   = "inward",
  hjust   = "inward",
  size    = 3
)

# save plot in '.png' format
ggsave(
  filename = "figureS8.png",
  path     = "~/documents/rcoimbra_phd/project_kenya/figures/revised/",
  width    = 210,
  height   = 74,
  units    = "mm",
  dpi      = 300
)
