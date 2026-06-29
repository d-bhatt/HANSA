# HANSA — Homologous Antigen Similarity Analysis
# Visualization code to make Circos plot
# Author: Darshak K. Bhatt
# Affiliation: University Medical Center Groningen, University of Groningen
# License: MIT
# Version: 1.0.0
# Date: 2026-06-29

#==========================================================================

library(tidyverse)
library(circlize)
library(RColorBrewer)
library(viridisLite)
library(ComplexHeatmap)
library(grid)

# ----------------------------
# Read and summarize data if required
# ----------------------------
df <- read_tsv("crossreactivity_results.tsv")

links <- df %>%
  mutate(weight = 1)

#links <- df %>%                          # turn on if needing to summarize data
#  group_by(seqA_name, seqB_name) %>%
#  summarise(
#    n = n(),
#    mean_mismatch = mean(mismatch, na.rm = TRUE),
#    .groups = "drop"
#  )

# ----------------------------
# Nodes and sector colors
# ----------------------------
nodes <- sort(unique(c(links$seqA_name, links$seqB_name)))

grid.col <- setNames(
  colorRampPalette(brewer.pal(12, "Set2"))(length(nodes)),
  nodes
)

# ----------------------------
# Link colors = same as node colors
# ----------------------------
link_cols <- grid.col[links$seqA_name]


# ----------------------------
# Labels with total counts
# ----------------------------
totals <- bind_rows(
  links %>% transmute(node = seqA_name, weight),
  links %>% transmute(node = seqB_name, weight)
) %>%
  group_by(node) %>%
  summarise(total = sum(weight), .groups = "drop")

# For summary data
#totals <- bind_rows(
#  links %>% transmute(node = seqA_name, n),
#  links %>% transmute(node = seqB_name, n)
#) %>%
#  group_by(node) %>%
#  summarise(total = sum(n), .groups = "drop")

label_map <- setNames(
  paste0(totals$node, " (", totals$total, ")"),
  totals$node
)

# ----------------------------
# Viridis color mapping
# ----------------------------
#col_var <- links$mean_mismatch       # for summary data

#col_var <- links$pct_id

#gradient color
#link_cols_fun <- colorRamp2(
#  seq(min(col_var), max(col_var), length.out = 5),
#  viridisLite::viridis(5)
#)
#link_cols <- link_cols_fun(col_var)

#for categorical data
#col_var <- links$tier   # your categorical variable
#categories <- unique(col_var)
#cat_colors <- setNames(
#  viridisLite::viridis(length(categories), option = "D"),
#  categories
#)
#link_cols <- cat_colors[col_var]



# for summary data
#breaks <- seq(
#  min(col_var, na.rm = TRUE),
#  max(col_var, na.rm = TRUE),
#  length.out = 5
#)

#vir_cols <- viridis(
#  length(breaks),
#  option = "D"
#)

#link_cols_fun <- colorRamp2(
#  breaks,
#  vir_cols
#)

#link_cols <- link_cols_fun(col_var)


# ----------------------------
# Plot
# ----------------------------
circos.clear()

circos.par(
  gap.degree = 10,
  start.degree = 90
)

chordDiagram(
#  x = links[, c("seqA_name", "seqB_name", "n")],     #for summary data
   x = links[, c("seqA_name", "seqB_name", "weight")],
  grid.col = grid.col,
  col = link_cols,
  transparency = 0.5,
  annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.15),
  link.sort = TRUE,
  link.largest.ontop = TRUE
)

# ----------------------------
# Labels
# ----------------------------
spaced_label <- function(x, space = " ") {
  paste(strsplit(x, "")[[1]], collapse = space)
}

circos.trackPlotRegion(
  track.index = 1,
  bg.border = NA,
  panel.fun = function(x, y) {
    
    sector.name <- get.cell.meta.data("sector.index")
    xlim <- get.cell.meta.data("xlim")
    ylim <- get.cell.meta.data("ylim")
    
    circos.text(
      x = mean(xlim),
      y = ylim[1] + 0.1,
      #labels = label_map[sector.name],
      
      label <- spaced_label(label_map[sector.name], space = "\u00A0"),  # spaced labels
      
      facing = "bending.outside",
      niceFacing = TRUE,
      adj = c(0.5, 1),
      cex = 1
    )
  }
)

#title(
#  "Peptide Cross-Reactivity Chord Diagram\nChord color = Mean mismatch"
#)

# ----------------------------
# Legend
# ----------------------------
#lgd <- Legend(
#  title = "Variable",
#  col_fun = link_cols_fun,
#  at = round(pretty(range(col_var, na.rm = TRUE), n = 5), 2),
#  direction = "vertical"
#)

#draw(
#  lgd,
#  x = unit(0.9, "npc"),
#  y = unit(0.5, "npc"),
#  just = c("left", "center")
#)

# Optional cleanup after everything is drawn
 circos.clear()
