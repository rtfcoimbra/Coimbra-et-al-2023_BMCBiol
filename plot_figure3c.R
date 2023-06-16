library(tidyverse)
library(circlize)

# set working directory
setwd("~/documents/rcoimbra_phd/project_kenya/results/08_gene_flow/bayesass/")

# import migration matrix
data <- read.csv("run3.csv", row.names = 1, header = TRUE)

# transpose and convert migration matrix into long tibble
tbl <- t(data) %>%
  as_tibble(rownames = "origin") %>%
  pivot_longer(
    -origin,
    names_to  = "destination",
    values_to = "migration_rate"
  )

# set color palette
palette <- c(
  "Nubian"      = "#D55E00",
  "Reticulated" = "#CC79A7",
  "Masai"       = "#009E73"
)

# save plot in '.pdf' format
png(
  file   = "~/documents/rcoimbra_phd/project_kenya/figures/revised/figure3c.png",
  height = 210,
  width  = 210,
  units  = "mm",
  res    = 300
)

# set plot margins
par(mar = c(2, 0, 0, 0))

# set circos graphic parameters
circos.par(
  start.degree            = 90,
  gap.degree              = 4,
  track.margin            = c(-0.1, 0.1),
  cell.padding            = c(0, 0, 0, 0),
  points.overflow.warning = FALSE
)

# plot chord diagram
chordDiagram(
  x                     = tbl,
  order                 = c("Reticulated", "Masai", "Nubian"),
  grid.col              = palette,
  transparency          = 0.33,
  link.sort             = TRUE,
  link.decreasing       = TRUE,
  link.largest.ontop    = TRUE,
  link.arr.type         = "big.arrow",
  #link.visible          = tbl$migration_rate >= 0.01,
  directional           = 1,
  direction.type        = c("arrows", "diffHeight"),
  diffHeight            = -mm_h(1),
  annotationTrack       = "grid",
  annotationTrackHeight = c(0.05, 0.1)
)

# create plotting regions for cells in one single track
circos.track(
  track.index = 1,
  bg.border   = NA,
  panel.fun   = function(x, y) {
    xlim        <- get.cell.meta.data("xlim")
    sector.name <- get.cell.meta.data("sector.index")
    # add sector names
    circos.text(
      x          = mean(xlim),
      y          = 3.9,
      labels     = sector.name,
      facing     = "bending",
      niceFacing = TRUE,
      cex        = 1.8
    )
    # add axis
    circos.axis(
      h                 = "top",
      major.at          = seq(0, xlim[2] + 0.25, 0.25),
      minor.ticks       = 1,
      major.tick.length = mm_y(1.5),
      labels.niceFacing = FALSE,
      labels.cex        = 1.2,
      lwd               = 1.2
    )
  }
)

# reset graphic parameters and internal variables
circos.clear()

invisible(dev.off())