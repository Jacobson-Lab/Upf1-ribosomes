# ---------------------------------------
# Figure 4 -- Metagene bin across normalized CDS
# ---------------------------------------

library(dplyr)
library(ggplot2)

# Prepare Data for 4a-c
  # 4a: N-terminal (-CHX) data
dfa <- read.table("../data/Data_Figure4a_Nterm_CHXminus.txt", header = TRUE, sep = "\t")
dfa$size <- factor(dfa$size, levels = c("S", "M"))
dfa$strain2 <- paste0(dfa$strain, "\n(-CHX)")
dfa$strain2 <- factor(dfa$strain2, levels = c("WT + EV\n(-CHX)", "FLAG-UPF1\n(-CHX)"))

  # 4b: N-terminal (+CHX) data
dfb <- read.table("../data/Data_Figure4b_Nterm_CHXplus.txt", header = TRUE, sep = "\t")
dfb$size <- factor(dfb$size, levels = c("M", "L"))
dfb$strain2 <- paste0(dfb$strain, "\n(+CHX)")
dfb$strain2 <- factor(dfb$strain2, levels = c("WT + EV\n(+CHX)", "FLAG-UPF1\n(+CHX)"))

  # 4c: C-terminal (+CHX) data
dfc <- read.table("../data/Data_Figure4c_Cterm_CHXplus.txt", header = TRUE, sep = "\t")
dfc$size <- factor(dfc$size, levels = c("M", "L"))
dfc$strain2 <- paste0(dfc$strain, "\n(+CHX)")
dfc$strain2 <- factor(dfc$strain2, levels = c("WT + EV\n(+CHX)", "UPF1-FLAG\n(+CHX)", "UPF1-FLAG/upf2\u0394\n(+CHX)", "DE572AA-FLAG\n(+CHX)"))

# Plot for 4a-c (the same code is used to plot N- and C-terminal data). Assign the desired data.frame to 'data' argument. (Only turn on legend for c)
C <- ggplot(data = dfc, aes(x = bin, y = fraction_ave*100, color = Ribosomes, size = Ribosomes)) +
  geom_line() +
  geom_hline(yintercept = 1, color = "grey50", linetype = "longdash", size = 0.3) +
  facet_grid(size~strain2) +
  xlab("% CDS") + ylab("% Footprint count") +
  scale_color_manual(name = "Ribosomes", values = c(Total = "Purple", IP = "orange")) +
  scale_size_manual(values = c(IP = 0.3, Total = 0.5)) +
  coord_cartesian(ylim = c(0, 2.2)) +
  theme_bw(base_size = 10) + 
  theme(strip.background = element_blank(), strip.text.x = element_text(face = "italic"), strip.text.y = element_text(angle = 0),
        panel.grid = element_blank(), legend.position = "right", aspect.ratio = 1) +
  guides(color = guide_legend(order = 1, #keywidth = unit(1.5, "cm"), label.position = "bottom",
                              override.aes = list(size = 1)),
         size = "none")

# 4d: Main footprints in -CHX and +CHX
dfd <- bind_rows(list(`- (S)` = dfa[dfa$size == "S", ], `+ (M)` = dfb[dfb$size == "M", ]), .id = "CHX (size)")
dfd$Ribosomes <- factor(dfd$Ribosomes, levels = c("Total", "IP"))
dfd$`CHX (size)` <- factor(dfd$`CHX (size)`, levels = c("+ (M)", "- (S)"))
dfd$strain <- factor(dfd$strain, levels = c("WT + EV", "FLAG-UPF1"))

D <- ggplot(dfd, aes(x = bin, y = fraction_ave*100, color = `CHX (size)`)) +
  geom_line(size = 0.5) +
  geom_hline(yintercept = 1, color = "grey50", linetype = "longdash", size = 0.3) +
  facet_grid(Ribosomes~strain) +
  xlab("% CDS") + ylab("% Footprint count") +
  coord_cartesian(ylim = c(0, 2.2)) +
  theme_bw(base_size = 10) + 
  theme(strip.background = element_blank(), strip.text.x = element_text(face = "italic"), strip.text.y = element_text(angle = 0),
        panel.grid = element_blank(), legend.position = "right", aspect.ratio = 1) +
  guides(color = guide_legend(order = 1, #keywidth = unit(1.5, "cm"), label.position = "bottom",
                              override.aes = list(size = 1)))

# Combine panels
library(patchwork)
ABC <- ((A + B) / C) + plot_annotation(tag_levels = 'A')
p <- ABC / (D + plot_spacer() + plot_layout(widths = c(0.6, 0.4))) + 
  plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 12, face = "bold"))

# Export plot
library(Cairo)
CairoFonts(
  regular = "Arial:style=Regular",
  bold = "Arial:style=Black",
  italic = "Arial:style=Italic",
  bolditalic = "Arial:style=Black Italic",
  symbol = "Symbol"
)
cairo_pdf(filename = "Metagene_bin_7x10.pdf", family = "Arial", width = 7, height = 10) 
p
dev.off()

