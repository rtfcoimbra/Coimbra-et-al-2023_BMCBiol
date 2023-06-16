################################################################################
#                                  libraries                                   #
################################################################################

library(ggcorrplot)
library(patchwork)
library(tidyverse)

################################################################################
#                                  functions                                   #
################################################################################

# source modified r script from evaladmix
source("~/documents/rcoimbra_phd/project_kenya/results/05_admixture/visFuns_modified.R")

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

# return ancestry plot object
make_ancestry_plot <- function(tbl, palette){
  return(
    ggplot(tbl, aes(x = sample_id, y = probability, fill = ancestry)) +
      geom_col() +
      facet_grid(
        cols = vars(fct_inorder(location)),
        scales = 'free',
        space = 'free'
      ) +
      scale_y_continuous(expand = c(0, 0), breaks = c(seq(0, 1, 0.2))) +
      scale_fill_manual(values = palette, guide = "none") +
      labs(
        x = "Individuals",
        y = paste("K =", length(unique(tbl$ancestry)))
      ) +
      theme_minimal() +
      theme(
        panel.spacing.x     = unit(0, "lines"),
        axis.title          = element_text(size = 8),
        axis.text.x         = element_text(size = 3, angle = 90, vjust = 0.5, hjust = 1, margin = margin(t = -1)),
        axis.text.y         = element_text(size = 5),
        axis.ticks.length.y = unit(0, "lines"),
        strip.text          = element_blank(),
        strip.background    = element_blank(),
        panel.border        = element_rect(linewidth = 0.6, color = "grey80", fill = NA),
        panel.grid          = element_blank()
      )
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

# set color palette
palette <- set_names(c("#d55e00", "#cc79a7", "#009e73"), lvls.1)
                     
# set shapes
shapes <- set_names(c(18, 15, 17), lvls.1)

################################################################################
#                                 preparation                                  #
################################################################################

# set working directory
setwd("~/documents/rcoimbra_phd/project_kenya")

# read bamlist used with ANGSD
bams <- read_table("results/04_pca/bamlist", col_names = FALSE)

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
#                                     pca                                      #
################################################################################

# set working directory
setwd("~/documents/rcoimbra_phd/project_kenya/results/04_pca/")

# read covariance matrix generated with PCAngsd
cov.mat <- read_delim("snps.hwe_filter.cov", delim = " ", col_names = FALSE)

# perform PCA
pca <- prcomp(cov.mat, scale = TRUE)

# bind columns of PC coordinates to tibble with sample information
tbl <- sample.info %>% bind_cols(as_tibble(pca$x))

# caculate eigenvalues and percentage of variance explained
eigenvals <- pca$sdev ^ 2
eigen <- tibble(
  prin.comp   = c(seq(1, length(eigenvals))),
  var.percent = eigenvals / sum(eigenvals) * 100,
  cum.sum     = cumsum(var.percent)
)

# pca plot
panel.b <- ggplot(data = tbl, mapping = aes(x = PC1, y = PC2)) +
  geom_point(mapping = aes(color = taxonomy, shape = taxonomy), size = 1.5) +
  scale_shape_manual(values = shapes) +
  scale_color_manual(values = palette) +
  labs(
    x = paste0("PC1 (", round(eigen$var.percent[1], 2), "%)"),
    y = paste0("PC2 (", round(eigen$var.percent[2], 2), "%)")
  ) +
  theme_light() +
  theme(
    axis.title        = element_text(size = 8),
    axis.text         = element_text(size = 5),
    legend.title      = element_blank(),
    legend.text       = element_text(size = 7),
    legend.key.size   = unit(0.1, "lines"),
    legend.background = element_rect(fill = "transparent"),
    legend.position   = c(0.245, 0.135)
  ) +
  guides(shape = guide_legend(override.aes = list(size = 1.5)))

# clean environment
rm(sample.info, cov.mat)

################################################################################
#                                  admixture                                   #
################################################################################

# set working directory
setwd("~/documents/rcoimbra_phd/project_kenya/results/05_admixture/")

# get input file names
files <- str_sort(list.files(pattern = ".qopt"), numeric = TRUE)

# generate tibble of ancestry proportions
qopts <- files %>% map(make_ancestry_tibble, tbl = tbl)

# plot ancestry porportions
p.k3 <- make_ancestry_plot(qopts[[3]], c("#009f72", "#cc79a7", "#d55e00")) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank()
  )
p.k9  <- make_ancestry_plot(qopts[[9]],  c("#cc79a7", "#b44582", "#e18e4c", "#006c4d", "#954100", "#d55e00", "#b2e2d4", "#66c5aa", "#009f72"))

# arrange composite plot layout
panel.c <- p.k3 / (p.k9 + plot_layout(tag_level = 'new'))

################################################################################
#                               likelihoods per k                              #
################################################################################

# find input files
files <- str_sort(list.files(pattern = "likelihoods*"), numeric = TRUE)

# create a tibble of run likelihoods
tbl.likes <- files %>%
  map_dfr(read_table, col_names = FALSE, .id = "k") %>%
  rename(likelihoods = X1) %>%
  group_by(k) %>%
  summarize(
    n    = n(),
    mean = mean(likelihoods),
    sd   = sd(likelihoods)
  ) %>%
  mutate(
    k = fct_relevel(k, as.character(seq(1, length(k)))),
    se = sd / sqrt(n)
  )

# plot likelihoods across runs per k
panel.d <- ggplot(tbl.likes, aes(x = k, y = mean)) +
  geom_line(aes(group = 1), linewidth = 0.25, alpha = 0.5) +
  geom_point(size = 1.2, color = "red", alpha = 0.5) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.5, color = "red", alpha = 0.5) +
  labs(
    x = expression(italic("K")),
    y = expression(paste("Mean L(", italic("K"), ") \u00B1 SE"))
  ) +
  theme_light() +
  theme(
    axis.title      = element_text(size = 8),
    axis.text       = element_text(size = 5),
    axis.text.y     = element_text(angle = 90, hjust = 0.5),
    legend.position = "none"
  )

################################################################################
#                                  evaladmix                                   #
################################################################################

# find input files
files <- c(
  str_sort(list.files(pattern = ".qopt"), numeric = TRUE),
  str_sort(list.files(pattern = ".txt"), numeric = TRUE)
)

pop <- tbl %>% select(location, id) %>% as.data.frame()

for(i in c(3, 9)) {
  q <- read_table(files[i], col_names = FALSE) %>%
    select(!last_col()) %>%
    as.data.frame()
  ord <- orderInds(pop = as.vector(pop[, 1]), popord = lvls.2)
  r <- read.table(files[i + 14])
  assign(
    paste0("mat", i),
    plotCorRes(
      cor_mat        = r,
      pop            = as.vector(pop[, 1]),
      superpop       = as.vector(tbl$taxonomy),
      ord            = ord,
      title          = paste0("Correlation of residuals for K=", i),
      max_z          = 0.5,
      min_z          = -0.5,
      cex.main       = 1.2,
      cex.lab        = 0.4,
      cex.lab.2      = 1,
      cex.legend     = 1.2,
      rotatelabpop   = 90,
      adjlab         = 0.015,
      adjlabsuperpop = 0.165,
      pop_labels     = c(TRUE, TRUE),
      lineswidth     = 0.5,
      lineswidthsuperpop = 1.5
    )
  )
}

mat3 <- as.data.frame(mat3) %>%
  mutate_all(~ifelse(is.nan(.), NA, .))

mat9 <- as.data.frame(mat9) %>%
  mutate_all(~ifelse(is.nan(.), NA, .))

m3 <- ggcorrplot(mat3, type = "full", show.diag = TRUE, show.legend = FALSE) +
  scale_fill_gradientn(
    colors = c("#009ACD", "white", "#ED5504"),
    limit  = c(-0.5, 0.5)
  ) +
  labs(fill = "Population\nmean corr.") +
  geom_vline(xintercept = c(7.5, 18.5), linewidth = 0.5, color = "grey50") +
  geom_hline(yintercept = c(7.5, 18.5), linewidth = 0.5, color = "grey50") +
  theme(
    axis.text.x       = element_blank(),
    axis.text.y       = element_blank(),
    panel.background = element_rect(fill = "grey", colour = "transparent"),
    panel.grid       = element_blank()
  )

m9 <- ggcorrplot(mat9, type = "full", show.diag = TRUE) +
  scale_fill_gradientn(
    colors = c("#009ACD", "white", "#ED5504"),
    limit  = c(-0.5, 0.5)
  ) +
  labs(fill = "Population\nmean corr.") +
  geom_vline(xintercept = c(7.5, 18.5), linewidth = 0.5, color = "grey50") +
  geom_hline(yintercept = c(7.5, 18.5), linewidth = 0.5, color = "grey50") +
  theme(
    axis.text.x       = element_blank(),
    axis.text.y       = element_blank(),
    legend.title      = element_text(size = 8),
    legend.text       = element_text(size = 6),
    legend.key.width  = unit(1, "lines"),
    legend.key.height = unit(1, "lines"),
    legend.text.align = 1,
    panel.background  = element_rect(fill = "grey", colour = "transparent"),
    panel.grid        = element_blank()
  )

# arrange composite plot layout
panel.e <- (m3 | (m9 + plot_layout(tag_level = 'new'))) + plot_layout(guides = "collect")

################################################################################
#                                    figure                                    #
################################################################################

# set spacer plot that can be tagged
spacer <- ggplot() + theme_void()

# arrange figure layout
( spacer / (( panel.b | panel.c ) + plot_layout(ncol = 2, widths = c(1/4, 3/4))) / ((panel.d + panel.e + guide_area()) + plot_layout(ncol = 3, widths = c(1/3, 2/3, 0))) ) +
  plot_layout(nrow = 3, heights = c(2/3, 1/5, 1/5)) +
  plot_annotation(tag_levels = c("a"))

# save figure in '.pdf' format
ggsave(
  filename = "figure1.pdf",
  path     = "~/documents/rcoimbra_phd/project_kenya/figures/revised/",
  width    = 210,
  height   = 297,
  units    = "mm",
  dpi      = 300
)
