################################################################################
#                                  libraries                                   #
################################################################################

library(tidyverse)
library(patchwork)

################################################################################
#                                  functions                                   #
################################################################################

# return tibble of ancestry proportions
make_ancestry_tibble <- function(file, tbl){
  return(
    read_table(file, col_names = FALSE) %>%
      select(!last_col()) %>%
      mutate(
        sample_id  = tbl$id,
        subspecies = tbl$taxonomy,
        location   = tbl$location
      ) %>%
      pivot_longer(
        cols      = starts_with("X"),
        names_to  = "ancestry",
        values_to = "probability"
      ) %>%
      arrange(subspecies, location)
  )
}

# return plot object of ancestry proportions
make_ancestry_plot <- function(tbl, palette, theme_opt){
  return(
    ggplot(tbl, aes(x = sample_id, y = probability, fill = ancestry)) +
      geom_col(color = "grey", linewidth = 0.025) +
      facet_grid(~fct_inorder(location), scales = 'free', space = 'free') +
      scale_x_discrete(expand = expansion(add = 1)) +
      scale_y_continuous(expand = c(0, 0), breaks = c(seq(0, 1, 0.2))) +
      scale_fill_manual(values = palette, guide = "none") +
      labs(y = paste("K =", length(unique(tbl$ancestry)))) +
      theme_minimal() +
      theme_opt
  )
}

################################################################################
#                                 configurations                               #
################################################################################

# set subspecies factor levels
lvls.1 <- c("Nubian", "Reticulated", "Masai s. str.")

# set location factor levels
lvls.2 <- c(
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

# set theme modifiers for middle plot
theme_middle <- theme(
  panel.spacing.x     = unit(0, "lines"),
  axis.title.x        = element_blank(),
  axis.title.y        = element_text(size = 7),
  axis.text.x         = element_blank(),
  axis.text.y         = element_text(size = 4),
  axis.ticks.length.y = unit(0, "lines"),
  strip.background    = element_blank(),
  strip.text          = element_blank(),
  panel.grid          = element_blank()
)

# set theme modifiers for top plot
theme_top <- theme_middle + theme(
  strip.text = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 0)
)

# set theme modifiers for bottom plot
theme_bottom <- theme_middle + theme(
  axis.text.x = element_text(size = 4, angle = 90, vjust = 0.5, hjust = 1)
)

################################################################################
#                                 preparation                                  #
################################################################################

# set working directory
setwd("~/documents/rcoimbra_phd/project_kenya/")

# read bamlist used with ANGSD
bams <- read_table("results/05_admixture/bamlist", col_names = FALSE)

# import metadata
metadata <- read_csv("metadata.csv")

# create tibble with samples information
sample.info <- pull(bams, 1) %>%
  str_replace_all(c("/.*/" = "", ".clean.bam" = "")) %>%
  as_tibble_col(column_name = "id") %>%
  left_join(metadata) %>%
  mutate(
    taxonomy = case_when(
      str_detect(taxonomy, regex("Giraffa tippelskirchi tippelskirchi"))   ~ "Masai s. str.",
      str_detect(taxonomy, regex("Giraffa camelopardalis camelopardalis")) ~ "Nubian",
      str_detect(taxonomy, regex("Giraffa reticulata"))                    ~ "Reticulated",
      TRUE                                                                 ~ NA_character_
    )
  ) %>%
  mutate(
    taxonomy = fct_relevel(taxonomy, lvls.1),
    location = fct_relevel(location, lvls.2)
  )

# clean environment
rm(bams, metadata)

################################################################################
#                             admixture for k=2-11                             #
################################################################################

# set working directory
setwd("~/documents/rcoimbra_phd/project_kenya/results/05_admixture/")

# get input file names
files <- str_sort(list.files(pattern = ".qopt"), numeric = TRUE)

# generate tibbles of ancestry proportions
qopts <- files %>% map(make_ancestry_tibble, sample.info)

# generate ancestry plots (color order determined by trial and error)
p.k2  <- make_ancestry_plot(qopts[[2]],  c("#d55e00", "#009f72"), theme_top)
p.k3  <- make_ancestry_plot(qopts[[3]],  c("#009f72", "#cc79a7", "#d55e00"), theme_middle)
p.k4  <- make_ancestry_plot(qopts[[4]],  c("#009f72", "#d55e00", "#cc79a7", "#66c5aa"), theme_middle)
p.k5  <- make_ancestry_plot(qopts[[5]],  c("#cc79a7", "#d55e00", "#009f72", "#954100", "#66c5aa"), theme_middle)
p.k6  <- make_ancestry_plot(qopts[[6]],  c("#b44582", "#009f72", "#66c5aa", "#954100", "#cc79a7", "#d55e00"), theme_middle)
p.k7  <- make_ancestry_plot(qopts[[7]],  c("#954100", "#009f72", "#66c5aa", "#d55e00", "#b44582", "#e18e4c", "#cc79a7"), theme_middle)
p.k8  <- make_ancestry_plot(qopts[[8]],  c("#cc79a7", "#66c5aa", "#009f72", "#e18e4c", "#b2e2d4", "#b44582", "#d55e00", "#954100"), theme_middle)
p.k9  <- make_ancestry_plot(qopts[[9]],  c("#cc79a7", "#b44582", "#e18e4c", "#006c4d", "#954100", "#d55e00", "#b2e2d4", "#66c5aa", "#009f72"), theme_middle)
p.k10 <- make_ancestry_plot(qopts[[10]], c("#66c5aa", "#e18e4c", "#893c00", "#b2e2d4", "#006c4d", "#009f72", "#cc79a7", "#b44582", "#e1b0cb", "#d55e00"), theme_middle)
p.k11 <- make_ancestry_plot(qopts[[11]], c("#893c00", "#b44582", "#006c4d", "#009f72", "#e1b0cb", "#b2e2d4", "#cc79a7", "#66c5aa", "red", "#e18e4c", "#d55e00"), theme_bottom)

# clean environment
rm(sample.info)

################################################################################
#                                  final plot                                  #
################################################################################

# generate composite plot
p.k2 / p.k3 / p.k4 / p.k5 / p.k6 / p.k7 / p.k8 / p.k9 / p.k10 / p.k11

# save plot in '.png' format
ggsave(
  filename = "figureS3.png",
  path     = "~/documents/rcoimbra_phd/project_kenya/figures/revised",
  width    = 210,
  height   = 230,
  units    = "mm",
  dpi      = 300
)
