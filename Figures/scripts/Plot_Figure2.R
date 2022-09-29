# ---------------------------------------
# Figure 2 -- footprint length distribution in ribosome profiling libraries
# ---------------------------------------

library(dplyr)
library(ggplot2)

# Prepare Data
  # N-terminal data
dfn <- read.table("../data/Data_Figure2_Nterm.txt", header = TRUE, sep = "\t")
dfn$strain <- recode(dfn$strain, `FLAG-UPF1 UPF2/UPF3 EE` = "FLAG-UPF1\nUPF2/UPF3 EE")
dfn <- dfn[, c("length", "fraction_average", "Ribosomes", "CHX", "strain")]
dfn$strain <- factor(dfn$strain, levels = c("WT + EV", "FLAG-UPF1", "FLAG-UPF1\nUPF2/UPF3 EE"))
dfn$CHX <- factor(dfn$CHX, levels = c("+", "-"))

  # C-terminal data
dfc <- read.table("../data/Data_Figure2_Cterm.txt", header = TRUE, sep = "\t")
dfc <- dfc[, c("length", "fraction_average", "Ribosomes", "CHX", "strain")]
dfc$strain <- factor(dfc$strain, levels = c("WT + EV", "UPF1-FLAG", "UPF1-FLAG/upf2\u0394", "DE572AA-FLAG"))


# Plot (the same code is used to plot N- and C-terminal data)
size_labels <- data.frame(strain = "WT + EV", lab = c("S", "M", "L"), xpos = c(21.5, 29.5, 40))
size_labels$strain <- factor(size_labels$strain, levels = c("WT + EV"))
df <- dfc # dfn or dfc
B <- ggplot(df) + # assign plot to A or B
  annotate("rect", ymin = -Inf, ymax = Inf, xmin = 20 , xmax = 23, fill = "grey60", alpha = 0.25) + # 20-23 nt (S)
  annotate("rect", ymin = -Inf, ymax = Inf, xmin = 27 , xmax = 32, fill = "grey60", alpha = 0.25) + # 27-32 nt (M)
  annotate("rect", ymin = -Inf, ymax = Inf, xmin = 37 , xmax = 43, fill = "grey60", alpha = 0.25) + # 37-43 nt (L)
  geom_text(data = size_labels, aes(x = xpos, y = 40, label = lab), size = 10/.pt, fontface = "bold", vjust = 0, nudge_y = 4) +
  geom_line(aes(x = as.numeric(as.character(length)), y = fraction_average*100, linetype = Ribosomes, color = CHX), alpha = 1, size = 0.5) +
  geom_label(data = data.frame(strain = unique(df[, c("strain")])),
             aes(label = strain, x = 45, y = 25), fontface = "italic", size = 10/.pt, hjust = 0, vjust = 0.5) +
  facet_grid(strain~.) +
  scale_linetype_manual(values = c(Total = "dashed", IP = "solid"), limits = c("Total", "IP")) +
  scale_color_manual(values = c(`+` = "#F8766D", `-` = "#00BFC4"), limits = c("+", "-")) +
  xlab("Footprint length (nt)") + ylab("% Footprint count") +
  coord_fixed(ylim = c(0, 40), ratio = 1/5, clip = "off") +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(), strip.background = element_rect(fill = "white"), strip.text = element_blank(),
        legend.position = "top", legend.box = "horizontal") +
  guides(color = guide_legend(title = "CHX", order = 1, keywidth = unit(1, "cm"), label.position = "bottom", override.aes = list(size = 1)),
         linetype = guide_legend(order = 2, keywidth = unit(1, "cm"), label.position = "bottom"), override.aes = list(size = 1))

# Combine panels
library(patchwork)
p <- (A / B) + 
  plot_layout(guides = "collect", heights = c(3, 4)) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 12, face = "bold"), legend.position = "top", plot.tag.position = "topleft")

# Export plot
library(Cairo)
CairoFonts(
  regular = "Arial:style=Regular",
  bold = "Arial:style=Black",
  italic = "Arial:style=Italic",
  bolditalic = "Arial:style=Black Italic",
  symbol = "Symbol"
)
cairo_pdf(filename = "Footprint_length_distribution_5x7.pdf", family = "Arial", width = 5, height = 7) 
p
dev.off()
