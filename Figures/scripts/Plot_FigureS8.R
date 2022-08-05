# ---------------------------------------
# Figure S8
#   A-B: Fractions of P-site region
#   C-D: P-site reading frame fractions
#   Calculated with only relevant footprint lengths with corrected/verified P-site offsets (20-23nt, 27-32nt, 37-43nt)
# ---------------------------------------

library(dplyr)
library(ggplot2)
library(ggh4x)
library(ggpubr)

# A-B -----------------------------------
# Prepare Data
  # N-terminal data
dfa <- read.table("../data/Data_FigureS8a_Nterm.txt", header = TRUE, sep = "\t")
dfa$strain <- recode(dfa$strain, `FLAG-UPF1 UPF2/3 EE` = "FLAG-UPF1\nUPF2/3 EE")
dfa$strain <- factor(dfa$strain, levels = sort(unique(dfa$strain))[c(3, 1, 2)])
dfa$psite_region <- recode_factor(dfa$psite_region, `5utr` = "5'-UTR", cds = "CDS", `3utr` = "3'-UTR")
dfa$riborep <- sub("_.*_", " rep ", dfa$sample)

  # C-terminal data
dfb <- read.table("../data/Data_FigureS8b_Cterm.txt", header = TRUE, sep = "\t")
dfb$strain <- factor(dfb$strain, levels = sort(unique(dfb$strain))[c(4, 2, 3, 1)])
dfb$psite_region <- recode_factor(dfb$psite_region, `5utr` = "5'-UTR", cds = "CDS", `3utr` = "3'-UTR")
dfb$riborep <- sub("_.*_", " rep ", dfb$sample)
dfb$CHX <- "+CHX"

# Plot (the same code is used to plot N- and C-terminal data). Assign the desired data.frame to 'data' argument.
pb <- ggplot(data = dfb) +
  geom_tile(aes(x = psite_region, y = reorder(riborep, dplyr::desc(riborep)), fill = fraction*100), height = 0.9, width = 0.9, color = NA) +
  geom_text(aes(x = psite_region, y = reorder(riborep, dplyr::desc(riborep)), label = round(fraction*100, digits = 2), color = psite_region), 
            size = 1.5, show.legend = FALSE) +
  facet_nested(strain+CHX~., space = "free", scales = "free") +
  scale_fill_distiller(name = "% of footprints ", type = "seq", palette = "Blues", direction = 1, limit = c(0, 100)) +
  scale_color_manual(name = "", values = c("black", "white", "black")) +
  ylab("") + xlab("P-site region") +
  theme_bw(base_size = 8) + 
  theme(panel.grid = element_blank(), strip.text.y = element_text(face = "italic"),
        strip.background = element_rect(fill = "white"), strip.placement = "outside",
        legend.position = "top", legend.key.height = unit(0.25, "cm"))

# C-D -----------------------------------
# Prepare Data
  # N-terminal data
dfc <- read.table("../data/Data_FigureS8c_Nterm.txt", header = TRUE, sep = "\t")
dfc$strain <- recode(dfc$strain, `FLAG-UPF1 UPF2/3 EE` = "FLAG-UPF1\nUPF2/3 EE")
dfc$strain <- factor(dfc$strain, levels = sort(unique(dfc$strain))[c(3, 1, 2)])
dfc$psite_region <- recode_factor(dfc$psite_region, `5utr` = "5'-UTR", cds = "CDS", `3utr` = "3'-UTR")
dfc$riborep <- sub("_.*_", " rep ", dfc$sample)

  # C-terminal data
dfd <- read.table("../data/Data_FigureS8d_Cterm.txt", header = TRUE, sep = "\t")
dfd$strain <- factor(dfd$strain, levels = sort(unique(dfd$strain))[c(4, 2, 3, 1)])
dfd$psite_region <- recode_factor(dfd$psite_region, `5utr` = "5'-UTR", cds = "CDS", `3utr` = "3'-UTR")
dfd$riborep <- sub("_.*_", " rep ", dfd$sample)
dfd$CHX <- "+CHX"

# Plot (the same code is used to plot N- and C-terminal data). Assign the desired data.frame to 'data' argument.
pd <- ggplot(data = dfd) +
  geom_tile(aes(x = Frame, y = reorder(riborep, dplyr::desc(riborep)), fill = fraction_by_region*100), height = 0.9, width = 0.9, color = NA) +
  geom_text(aes(x = Frame, y = reorder(riborep, dplyr::desc(riborep)), label = round(fraction_by_region*100, digits = 2)),
            size = 1.5, show.legend = FALSE) +
  facet_nested(strain+CHX~psite_region, space = "free", scales = "free") +
  scale_fill_distiller(name = "% of footprints ", type = "seq", palette = "Blues", direction = 1) +
  ylab("") +
  theme_bw(base_size = 8) + 
  theme(panel.grid = element_blank(), strip.text.y = element_text(face = "italic"),
        strip.background = element_rect(fill = "white"), strip.placement = "outside",
        legend.position = "top", legend.key.height = unit(0.25, "cm"))

# ---------------------------------------
# Combine panels
library(patchwork)
s8 <- pa + pb + pc + pd +
  plot_layout(heights = c(16, 21), widths = c(3, 5), byrow = FALSE) +
  plot_annotation(tag_levels = 'A')  & theme(plot.tag = element_text(size = 9, face = "bold"))

# Export plot
library(Cairo)
CairoFonts(
  regular = "Arial:style=Regular",
  bold = "Arial:style=Black",
  italic = "Arial:style=Italic",
  bolditalic = "Arial:style=Black Italic",
  symbol = "Symbol"
)
cairo_pdf(filename = "psite_region_reading_frame_6x10.pdf", family = "Arial", width = 6, height = 10) 
s8
dev.off()
