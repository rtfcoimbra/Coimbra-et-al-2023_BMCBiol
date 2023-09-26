################################################################################
#                                  libraries                                   #
################################################################################

library(tidyverse)

################################################################################
#                                  functions                                   #
################################################################################

# source modified r script from evaladmix
source("~/documents/rcoimbra_phd/project_kenya/results/05_admixture/visFuns_modified.R")

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

################################################################################
#                                 preparation                                  #
################################################################################

# set working directory
setwd("~/documents/rcoimbra_phd/project_kenya")

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
#                                plot evaladmix                                #
################################################################################

# set working directory
setwd("~/documents/rcoimbra_phd/project_kenya/results/05_admixture/")

# get input file names
files <- c(
  str_sort(list.files(pattern = ".qopt"), numeric = TRUE),
  str_sort(list.files(pattern = ".txt"), numeric = TRUE)
)

pop <- sample.info %>% select(location, id) %>% as.data.frame()

png(height = 9, width = 10, units = "in", res = 300)

for(i in 1:11) {
  q <- read_table(files[i], col_names = FALSE) %>%
    select(!last_col()) %>%
    as.data.frame()
  ord <- orderInds(pop = as.vector(pop[, 1]), popord = lvls.2)
  r <- read.table(files[i + 11])
  assign(
    paste0("mat", i),
    plotCorRes(
      cor_mat        = r,
      pop            = as.vector(pop[, 1]),
      superpop       = as.vector(sample.info$taxonomy),
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

invisible(dev.off())
