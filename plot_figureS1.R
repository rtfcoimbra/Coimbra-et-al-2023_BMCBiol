################################################################################
#                                  libraries                                   #
################################################################################

library(tidyverse)
library(ggasym)
library(ggcorrplot)
library(patchwork)

################################################################################
#                                  functions                                   #
################################################################################

# returns a tibble without observations that include the
# individual passed in the arguments
filter_tibble <- function(ind, dataset) {
  return(dataset %>% filter(ind1 != ind & ind2 != ind))
}

# counts the number of observations to be removed from tibble
# when filtering out the individual passed in the arguments
removal_counter <- function(ind, dataset) {
  clean_tbl <- dataset %>% filter(ind1 == ind | ind2 == ind)
  return(nrow(clean_tbl))
}

################################################################################
#                                 configurations                               #
################################################################################

# set location factor levels
lvls <- c(
  "Gambella NP",
  "Murchison Falls NP",
  "Ruma NP",
  "Kigio Wildlife Conservancy",
  "Lake Nakuru NP",
  "Soysambu Conservancy",
  "Mwea NR",
  "Aberdare Country Club",
  "El Karama Conservancy",
  "Loisaba Conservancy",
  "Mpala Conservancy",
  "Ol Pejeta Conservancy",
  "Solio GR",
  "Buffalo Springs NR",
  "Meru NP",
  "Samburu NR",
  "Garissa",
  "Ishaqbini Conservancy",
  "Selous GR",
  "Masai Mara",
  "Hell's Gate NP",
  "Naivasha (Private Ranches)",
  "Ngong",
  "Amboseli NP",
  "Mbirikani",
  "Nairobi NP",
  "Tsavo West NP",
  "Tsavo East NP"
)

################################################################################
#                                 preparation                                  #
################################################################################

# set working directory
setwd("~/documents/rcoimbra_phd/project_kenya/results/03_relatedness/")

# read bamlist used with ANGSD
bams <- read_table("bamlist", col_names = FALSE)

# import metadata
metadata <- read_csv("~/documents/rcoimbra_phd/project_kenya/metadata.csv")

# create tibble with samples information
sample.info <- pull(bams, 1) %>%
  str_replace_all(c("/.*/" = "", ".clean.bam" = "")) %>%
  as_tibble_col(column_name = "id") %>%
  left_join(metadata) %>%
  select(id, location) %>%
  mutate(
    index    = row_number(),
    location = fct_relevel(location, lvls)
  )

# read pairwise relatedness inferred by ngsremix
data <- read_tsv("ngsremix")

# replace individual number IDs with their respective sample names
data$ind1 <- sample.info$id[match(data$ind1, sample.info$index)]
data$ind2 <- sample.info$id[match(data$ind2, sample.info$index)]

# add missing comparisons and convert data to wide format (matrix)
tmp <- data %>%
  select(ind1, ind2, k0) %>%
  asymmetrise(ind1, ind2) %>%
  pivot_wider(
    names_from  = ind1,
    values_from = k0
  )

# arrange wide tibble by sampling location and convert to dataframe
df <- sample.info %>%
  rename(ind2 = id) %>%
  left_join(tmp) %>%
  arrange(location) %>%
  as.data.frame()

# reorder dataframe columns following the row order in column ind2
mat <- df[match(df$ind2, colnames(df))]

# set row names equal to column names
rownames(mat) <- colnames(mat)

# clean environment
rm(bams, metadata, sample.info, tmp, df)

################################################################################
#                                 relatedness                                  #
################################################################################

# plot pairwise k0 as a proxy for relatedness
p1 <- ggcorrplot(mat, type = "full", tl.cex = 3.5, tl.srt = 90) +
  geom_vline(xintercept = c(37.5, 76.5), linewidth = 0.3, color = "grey") +
  geom_hline(yintercept = c(37.5, 76.5), linewidth = 0.3, color = "grey") +
  scale_fill_viridis_c(option = "plasma", direction = -1, limit = c(0, 1)) +
  labs(fill = expression(italic("k"[0]))) +
  theme(
    legend.title = element_text(size = 8),
    legend.text  = element_text(size = 6),
    panel.grid   = element_blank()
  )

# create tibble to store the number and sample names
# of individuals removed per k0 threshold
removals.by.k0 <- tibble(
  k0_threshold   = double(),
  n_inds_removed = integer(),
  inds_removed   = list()
)

# iterate over k0 threshold values
for (i in seq(0, 1, 0.05)) {
  
  # get subset of data containing only observations for which k0 <= threshold
  data.subset <- data %>% filter(k0 <= i)
  
  # vector to store sample names of individuals removed for a given k0 threshold
  removed.inds <- c()
  
  # while there are observations in the tibble, perform the following actions
  while (nrow(data.subset) > 0) {
    
    # get individuals present in the filtered tibble
    candidates <- unique(append(data.subset$ind1, data.subset$ind2))
    
    # get the sample name of the individual that, when filtered out,
    # maximizes the number of observation removals from the tibble
    removed.ind <- candidates[
      which.max(map_int(candidates, removal_counter, dataset = data.subset))
    ]
    
    # append the selected individual to the vector of removed individuals
    removed.inds <- append(removed.inds, removed.ind)
    
    # remove observations containing the selected individual from the tibble
    data.subset <- filter_tibble(removed.ind, data.subset)
  }
  # append row to tibble of number and IDs of individuals removed per k0 threshold
  removals.by.k0 <- removals.by.k0 %>%
    add_row(
      k0_threshold   = i,
      n_inds_removed = length(removed.inds),
      inds_removed   = list(removed.inds)
    )
}

# plot number of individuals to remove by k0 threshold
p2 <- ggplot(removals.by.k0, aes(x = k0_threshold, y = n_inds_removed)) +
  geom_line(linewidth = 0.3) +
  geom_vline(xintercept = 0.45, color = "red", lty = 2, linewidth = 0.5) +
  labs(
    x = expression(italic("k"[0])),
    y = "No. of individuals to remove"
  ) +
  theme_minimal() +
  theme(
    axis.title      = element_text(size = 10),
    axis.text       = element_text(size = 8),
    legend.position = "none"
  )

# assign kinship status to original data tibble
data$status <- if_else(
  (data$ind1 %in% removals.by.k0$inds_removed[[10]]) | (data$ind2 %in% removals.by.k0$inds_removed[[10]]),
  "<= 0.45",
  "> 0.45"
)
# color palette
palette <- c(
  "<= 0.45" = "red",
  "> 0.45"  = "grey50"
)

# plot pairwise relatedness as k1 by k2
p3 <- ggplot(data = data, mapping = aes(x = k1, y = k2)) +
  geom_point(size = 1, aes(colour = status), alpha = 0.5) +
  scale_color_manual(name = expression(italic("k"[0])), values = palette) +
  labs(
    x = expression(italic("k"[1])),
    y = expression(italic("k"[2]))
  ) +
  lims(
    x = c(0, 1),
    y = c(0, 1)
  ) +
  theme_minimal() +
  theme(
    axis.title   = element_text(size = 10),
    axis.text    = element_text(size = 8),
    legend.title = element_text(size = 8, face = "bold"),
    legend.text  = element_text(size = 8)
  )

# save list of individuals to be removed from the data set
as_tibble(str_sort(removals.by.k0$inds_removed[[10]], numeric = TRUE)) %>%
  write_tsv("samples_to_remove.txt", col_names = FALSE)

################################################################################
#                                    figure                                    #
################################################################################

# arrange figure layout
p1 / (p2 | p3) + plot_layout(nrow = 2, heights = c(2/3, 1/3)) + plot_annotation(tag_levels = "a")

# save plot in '.png' format
ggsave(
  filename = "figureS1.png",
  path     = "~/documents/rcoimbra_phd/project_kenya/figures/revised/",
  width    = 210,
  height   = 270,
  units    = "mm",
  dpi      = 600
)
