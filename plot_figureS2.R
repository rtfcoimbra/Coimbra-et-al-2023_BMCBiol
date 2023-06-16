################################################################################
#                                  libraries                                   #
################################################################################

library(tidyverse)
library(patchwork)

################################################################################
#                                  functions                                   #
################################################################################

# return PCA plot object
make_pca_plot <- function(tbl, x, y){
  return(
    ggplot(data = tbl, mapping = aes(x = {{x}}, y = {{y}})) +
      geom_point(mapping = aes(color = location, shape = taxonomy), size = 1.2) +
      scale_shape_manual(name = "Taxonomy", values = shapes) +
      scale_color_manual(name = "Location", values = palette) +
      theme_light() +
      theme(
        axis.title = element_text(size = 8),
        axis.text  = element_blank(),
        axis.ticks = element_blank()
      ) +
      guides(
        color = guide_legend(title.position = "top"),
        shape = guide_legend(title.position = "top", direction = "vertical")
      )
  )
}

################################################################################
#                                 configurations                               #
################################################################################

# set subspecies factor levels
lvls.1 <- c(
  "Nubian",
  "Reticulated",
  "Masai s. str."
)

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

# set shapes
shapes <- set_names(c(18, 15, 17), lvls.1)

# set color palette
nub.cols <- colorRampPalette(colors = c("#893c00", "#d55e00", "#e18e4c", "#ffaf70"))(7)
ret.cols <- colorRampPalette(colors = c("#b44582", "#cc79a7", "#e1b0cb", "#efd5e3"))(11)
mas.cols <- colorRampPalette(colors = c("#006c4d", "#009f72", "#66c5aa", "#b2e2d4"))(10)
palette <- set_names(c(nub.cols, ret.cols, mas.cols), lvls.2)

# set theme modifiers for pca plots
no.axis.title.x <- theme(axis.title.x = element_blank())
no.axis.title.y <- theme(axis.title.y = element_blank())
no.axis.title   <- theme(axis.title = element_blank())

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

# scree plot
p0 <- ggplot(eigen, aes(x = prin.comp, y = var.percent)) +
  geom_col(fill = heat.colors(length(eigenvals))) +
  labs(
    x = "Principal components",
    y = "Explained variance (%)"
  ) +
  theme_light() +
  theme(
    axis.title = element_text(size = 8),
    axis.text  = element_text(size = 6)
  )

# pca plots
p12 <- make_pca_plot(tbl, PC1, PC2) + no.axis.title.x
p13 <- make_pca_plot(tbl, PC1, PC3) + no.axis.title.x
p14 <- make_pca_plot(tbl, PC1, PC4)
p23 <- make_pca_plot(tbl, PC2, PC3) + no.axis.title
p24 <- make_pca_plot(tbl, PC2, PC4) + no.axis.title.y
p34 <- make_pca_plot(tbl, PC3, PC4) + no.axis.title.y

# clean environment
rm(sample.info, cov.mat)

################################################################################
#                                    figure                                    #
################################################################################

# arrange composite plot layout
wrap_plots(
  A = p12, B = p13, C = p14, D = p23, E = p24, F = p34, G = p0,
  guides = "collect",
  design = "A#G
            BD#
            CEF"
) & theme(
  legend.title    = element_text(size = 7, face = "bold"),
  legend.text     = element_text(size = 7),
  legend.key.size = unit(0.5, "lines"),
  legend.position = "bottom"
)

# save figure in '.png' format
ggsave(
  filename = "figureS2.png",
  path     = "~/documents/rcoimbra_phd/project_kenya/figures/revised",
  width    = 210,
  height   = 230,
  units    = "mm",
  dpi      = 300
)