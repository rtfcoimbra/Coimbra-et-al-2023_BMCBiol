library(scales)
library(tidyverse)

# set working directory
setwd("~/documents/rcoimbra_phd/project_kenya/results/09_demographic_inference/")

# find stairway plot files
files <- list.files(pattern = "*.final.summary")

# set factor levels
fct.lvls <- c("Nubian", "Reticulated", "Masai s. str.")

# set color palette
palette <- set_names(c("#d55e00", "#cc79a7", "#009e73"), fct.lvls)

# create a tibble of effective population sizes over time
tbl <- files %>%
  set_names(str_to_sentence(str_remove(files, ".final.summary"))) %>%
  map_dfr(read_table2, .id = "subspecies") %>%
  select(!(mutation_per_site:`theta_per_site_97.5%`)) %>%
  mutate(subspecies = replace(subspecies, subspecies == "Masai", "Masai s. str.")) %>%
  mutate(subspecies = fct_relevel(subspecies, fct.lvls))

# generate stairway plot with axis in log scale
p <- ggplot(data = tbl, mapping = aes(x = year, y = Ne_median)) +
  # added shaded vertical areas
  annotate("rect", xmin = 11700, xmax = Inf, ymin = 0, ymax = Inf, fill = "grey", alpha = 0.25) +
  # add median line
  geom_line(mapping = aes(color = subspecies)) +
  # add shaded area for the 95% confidence interval
  geom_ribbon(mapping = aes(ymin = `Ne_2.5%`, ymax = `Ne_97.5%`, fill = subspecies), alpha = 0.3) +
  # set facets for subspecies
  facet_wrap(~ subspecies, nrow = 3, strip.position = "right") +
  # set color legend
  scale_colour_manual(values = palette) +
  scale_fill_manual(values = palette) +
  # convert scales to log format
  scale_x_log10(
    breaks = breaks_log(n = 6, base = 10),
    labels = label_number()
  ) +
  scale_y_log10(
    breaks = breaks_log(n = 6, base = 10),
    labels = label_number(),
    limits = c(0.1, 500000)
  ) +
  # set axes labels
  labs(
    x = "Years ago",
    y = expression(italic("N"[e]))
  ) +
  # adjust appearance
  theme_minimal() +
  theme(
    axis.title       = element_text(size = 10),
    axis.text        = element_text(size = 8),
    strip.text       = element_text(size = 10),
    panel.grid.minor = element_blank(),
    legend.position  = "none"
  )

# add log scale ticks to plot
p <- p + annotation_logticks(size = 0.2)

# add epoch label
p + annotate(geom = "text", x = 9700, y = 500000, label = "Holocene", hjust = "inward", size = 2.5) +
  annotate(geom = "text", x = 13700, y = 500000, label = "Pleistocene", hjust = "outward", size = 2.5)

# save plot in '.tiff' format
ggsave(
  filename = "figure4.tiff",
  path     = "~/documents/rcoimbra_phd/project_kenya/figures/revised/",
  width    = 170,
  height   = 170,
  units    = "mm",
  #bg       = "white",
  dpi      = 600
)
