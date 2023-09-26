################################################################################
#                                  libraries                                   #
################################################################################

library(tidyverse)
library(ape)
library(treeio)
library(ggtree)
library(patchwork)

################################################################################
#                                 configurations                               #
################################################################################

# set factor levels
lvls <- c(
  "West African",
  "Kordofan",
  "Nubian",
  "Reticulated",
  "Masai s. str.",
  "Luangwa",
  "South African",
  "Angolan"
)

# set color palette
palette <- c(
  "West African"  = "#F0E442",
  "Kordofan"      = "#E69F00",
  "Nubian"        = "#D55E00",
  "Reticulated"   = "#CC79A7",
  "Masai s. str." = "#009E73",
  "Luangwa"       = "#006046",
  "South African" = "#56B4E9",
  "Angolan"       = "#0072B2"
)

################################################################################
#                                 preparation                                  #
################################################################################

# set working directory
setwd("~/documents/rcoimbra_phd/project_kenya")

# import metadata
metadata <- read_csv("metadata.csv")

# create tibble with samples information
sample.info <- metadata %>%
  mutate(
    taxonomy = case_when(
      str_detect(taxonomy, regex("Giraffa camelopardalis peralta"))        ~ "West African",
      str_detect(taxonomy, regex("Giraffa camelopardalis antiquorum"))     ~ "Kordofan",
      str_detect(taxonomy, regex("Giraffa camelopardalis camelopardalis")) ~ "Nubian",
      str_detect(taxonomy, regex("Giraffa reticulata"))                    ~ "Reticulated",
      str_detect(taxonomy, regex("Giraffa tippelskirchi tippelskirchi"))   ~ "Masai s. str.",
      str_detect(taxonomy, regex("Giraffa tippelskirchi thornicrofti"))    ~ "Luangwa",
      str_detect(taxonomy, regex("Giraffa giraffa giraffa"))               ~ "South African",
      str_detect(taxonomy, regex("Giraffa giraffa angolensis"))            ~ "Angolan",
      TRUE                                                                 ~ NA_character_
    )
  ) %>%
  mutate(taxonomy = fct_relevel(taxonomy, lvls))

# clean environment
rm(metadata)

################################################################################
#                                mtdna phylogeny                               #
################################################################################

# set working directory
setwd("~/documents/rcoimbra_phd/project_kenya/results/07_phylogeny_mtdna/")

# import tree
t1 <- read.iqtree("partitions.nex.treefile")

# remove outgroup
t1.pruned <- drop.tip(t1, "JN632674")

# check the node labels
#plot.phylo(as.phylo(t1.pruned), no.margin = TRUE, cex = 0.75)
#nodelabels(cex = 0.75, frame = "circle")

# create the basic plot
t1.plot <- ggtree(tr = t1.pruned, layout = "roundrect", ladderize = TRUE)

# add taxonomy tibble to tree plot
t1.final <- t1.plot %<+% sample.info +
  # adjust x-axis plot limit
  xlim(0, 0.05) +
  # add scale bar
  geom_treescale(x = 0, y = 0, fontsize = 3) +
  # add tip labels
  geom_tiplab(mapping = aes(color = taxonomy), align = TRUE, size = 6 / .pt, linesize = 0.1, key_glyph = draw_key_label) +
  # add node points (support values)
  geom_point2(mapping = aes(subset = as.numeric(sub("/.*", "", label)) >= 80 & as.numeric(sub(".*/", "", label)) >= 95 & !isTip)) +
  # add branch labels (support values)
  #geom_text2(mapping = aes(x = branch,
  #                         label = label,
  #                         subset = as.numeric(sub("/.*", "", label)) >= 80 & as.numeric(sub(".*/", "", label)) >= 95 & !isTip),
  #           vjust = -0.5,
  #           size = 1.8) +
  # add legend
  scale_color_manual("Taxonomy", labels = lvls, values = palette) +
  theme(
    legend.title    = element_text(face = "bold", size = 10),
    legend.text     = element_text(size = 9), # face = "italic"
    legend.key.size = unit(0.8, "line"),
    legend.position = c(0.28, 0.92)
  ) +
  guides(color = guide_legend(override.aes = list(size = 4, linetype = 0)))

################################################################################
#                                 snp phylogeny                                #
################################################################################

# set working directory
setwd("~/documents/rcoimbra_phd/project_kenya/results/06_phylogeny_snp/")

# import tree
t2 <- read.iqtree("snps.sampled.min4.phy.varsites.phy.treefile")

# remove outgroup
t2.pruned <- drop.tip(t2, "WOAK")

# create the basic plot
t2.plot <- ggtree(tr = t2.pruned, layout = "roundrect", ladderize = TRUE)

# add taxonomy tibble to tree plot
t2.final <- t2.plot %<+% sample.info +
  # adjust x-axis plot limit
  xlim(0, 0.01) +
  # add scale bar
  geom_treescale(x = 0, y = 0, fontsize = 3) +
  # add tip labels
  geom_tiplab(mapping = aes(color = taxonomy), align = TRUE, size = 6 / .pt, linesize = 0.1) +
  # add node points (support values)
  geom_point2(mapping = aes(subset = as.numeric(sub("/.*", "", label)) >= 80 & as.numeric(sub(".*/", "", label)) >= 95 & !isTip)) +
  # add legend
  scale_color_manual("Taxonomy", labels = lvls, values = palette) +
  theme(legend.position = "none")

################################################################################
#                                  final plot                                  #
################################################################################

# create composite plot
t2.final + t1.final + plot_annotation(tag_levels = "a")

# save plot in '.tiff' format
ggsave(
  filename = "figure2.tiff",
  path     = "~/documents/rcoimbra_phd/project_kenya/figures/",
  width    = 8.5,
  height   = 11,
  units    = "in",
  dpi      = 600
)
