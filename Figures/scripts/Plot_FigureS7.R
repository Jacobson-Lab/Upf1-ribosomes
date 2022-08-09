# ---------------------------------------
# Figure S7 -- Metagene plot relative to start or stop codons
# ---------------------------------------

library(data.table)
library(dplyr)
library(ggplot2)
library(ggh4x)
library(ggpubr)

# Prepare Data
  # N-terminal data
dfn <- read.table("../data/Data_FigureS7_Nterm.txt", header = TRUE, sep = "\t")
dfn$strain <- factor(dfn$strain, levels = sort(unique(dfn$strain))[c(2, 1)])
dfn$size <- factor(dfn$size, levels = c("M", "S", "L")) # specify order of plotting
ylimit = 0.5

  # C-terminal data
dfc <- read.table("../data/Data_FigureS7_Cterm.txt", header = TRUE, sep = "\t")
dfc$strain <- factor(dfc$strain, levels = sort(unique(dfc$strain))[c(4, 2, 3, 1)])
dfc$size <- factor(dfc$size, levels = c("S", "M", "L")) # specify order of plotting
ylimit = 0.6

# Plot (the same code is used to plot N- and C-terminal data). Assign the desired data.frame to 'data' argument.
  # Adapted from riboWaltz's metaplots.R
cdsl <- 100

    # dashed line for in-frame nucleotide position
lines3nt <- data.table(reg = rep(c("Distance from start (nt)", "Distance from stop (nt)"), times = c(length(seq(3, cdsl, 3)), length(seq(-2, -cdsl, -3)))), 
                       line = c(seq(3, cdsl, 3), rev(seq(-2, -cdsl, -3))))
lines3nt2 <- bind_rows(list(IP = lines3nt, Total = lines3nt), .id = "ribo")

    # red line for start and stop codon
linered <- data.table(reg = rep(c("Distance from start (nt)", "Distance from stop (nt)"), times = 2), 
                      ribo = rep(c("IP", "Total"), each = 2),
                      line = rep(c(0, 1), times = 2))

A <- ggplot(data = dfn) +
  geom_vline(data = lines3nt2, aes(xintercept = line), linetype = "dotted", color = "grey50", size = 0.1) +
  geom_vline(data = linered, aes(xintercept = line), linetype = "solid", color = "grey50", size = 0.1) +
  geom_line(aes(x = distance, y = fraction_ave*100, color = size), size = 0.3) +
  facet_nested(strain+CHX~ribo+reg, switch = "x", scales = "free", nest_line = TRUE) +
  scale_color_manual(breaks = c("S", "M", "L"), values = c(S = "royalblue1", M = "tomato", L = "forestgreen"), name = "Footprint size") +
  ylab("% Footprint count") +
  coord_cartesian(ylim = c(0, ylimit)) + 
  theme_bw(base_size = 8) + 
  theme(panel.grid = element_blank(), axis.title.x = element_blank(), strip.text.y = element_text(face = "italic"),
        strip.background = element_blank(), strip.placement = "outside",
        legend.position = "top") +
  guides(color = guide_legend(order = 1, keywidth = unit(1.5, "cm"), label.position = "bottom",
                              override.aes = list(size = 1)))

# Combine panels
library(patchwork)
p <- (A / B) + 
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 9, face = "bold"), legend.position = "top", plot.tag.position = "topleft")

# Export plot
library(Cairo)
CairoFonts(
  regular = "Arial:style=Regular",
  bold = "Arial:style=Black",
  italic = "Arial:style=Italic",
  bolditalic = "Arial:style=Black Italic",
  symbol = "Symbol"
)
cairo_pdf(filename = "Metagene_7x10.pdf", family = "Arial", width = 7, height = 10) 
p
dev.off()
