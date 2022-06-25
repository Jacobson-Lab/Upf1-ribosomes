# ---------------------------------------
# Figure 7 -- Statistical comparison of 3'-UTR amounts in C-terminally FLAG-tagged strains
# ---------------------------------------

library(dplyr)
library(ggplot2)
library(ggh4x)
library(ggpubr)

# A: Compare percentage of 3'-UTR between IP and Total -----------------------------------
dfb <- read.table("../data/Data_FigureS8b_Cterm.txt", header = TRUE, sep = "\t")
dfb$strain <- factor(dfb$strain, levels = sort(unique(dfb$strain))[c(4, 2, 3, 1)])
dfb$psite_region <- recode_factor(dfb$psite_region, `5utr` = "5'-UTR", cds = "CDS", `3utr` = "3'-UTR")
dfb$riborep <- sub("_.*_", " rep ", dfb$sample)
dfb$CHX <- "+CHX"

utr3 <- dfb[which(dfb$psite_region == "3'-UTR"), c(1, 4, 5)]
utr3$ribo <- sub("_.*", "", utr3$sample)
utr3 <- utr3[which(utr3$strain != "WT + EV"), ] # WT doesn't have IP samples
utr3$ribo <- factor(utr3$ribo, levels = c("Total", "IP"))

utr3_p <- ggplot(utr3, aes(x = ribo, y = fraction*100, fill = ribo)) +
  stat_summary(fun = mean, geom = "bar", color = NA, position = "dodge") + # Add mean bar
  #stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 2)), size = 1.5, vjust = -0.3, color = "black") + # Add mean number above bar
  geom_point(color = "grey25", position = position_dodge2(width = 0.75), size = 0.75, alpha = 0.5, show.legend = FALSE) +
  stat_compare_means(method = "t.test", comparisons = list(c("IP", "Total")), size = 1.5, paired = TRUE) +
  facet_grid(.~strain) +
  scale_fill_manual(name = "", values = c(IP = "orange", Total = "purple"), drop = FALSE) +
  xlab("") + ylab("% Footprint in the 3'-UTR") +
  theme_bw(base_size = 8) + 
  theme(strip.background = element_rect(fill = "white"), panel.grid = element_blank(), strip.text.x = element_text(face = "italic"),
        legend.position = "none")

# B: Compare percentage of reading frame 0 in 3'-UTR between IP and Total -----------------------------------
dfd <- read.table("../data/Data_FigureS8d_Cterm.txt", header = TRUE, sep = "\t")
dfd$strain <- factor(dfd$strain, levels = sort(unique(dfd$strain))[c(4, 2, 3, 1)])
dfd$psite_region <- recode_factor(dfd$psite_region, `5utr` = "5'-UTR", cds = "CDS", `3utr` = "3'-UTR")
dfd$riborep <- sub("_.*_", " rep ", dfd$sample)
dfd$CHX <- "+CHX"

f0 <- dfd[which(dfd$psite_region == "3'-UTR" & dfd$Frame == 0), c(1, 7, 8)]
f0$ribo <- sub("_.*", "", f0$sample)
f0 <- f0[which(f0$strain != "WT + EV"), ] # WT doesn't have IP samples
f0$ribo <- factor(f0$ribo, levels = c("Total", "IP"))

f0_p <- ggplot(f0, aes(x = ribo, y = fraction_by_region*100, fill = ribo)) +
  stat_summary(fun = mean, geom = "bar", color = NA, position = "dodge") + # Add mean bar
  #stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 2)), size = 1.5, vjust = -0.3, color = "black") + # Add mean number above bar
  geom_point(color = "grey25", position = position_dodge2(width = 0.75), size = 0.75, alpha = 0.5, show.legend = FALSE) +
  geom_hline(yintercept = 100/3, linetype = "dashed", color = "grey50", size = 0.5) +
  stat_compare_means(method = "t.test", comparisons = list(c("IP", "Total")), size = 1.5, paired = TRUE) +
  facet_grid(.~strain) +
  scale_fill_manual(name = "", values = c(IP = "orange", Total = "purple"), drop = FALSE) +
  xlab("") + ylab("% Frame 0 footprint within the 3'-UTR") +
  theme_bw(base_size = 8) + 
  theme(strip.background = element_rect(fill = "white"), panel.grid = element_blank(), strip.text.x = element_text(face = "italic"),
        legend.position = "none")

# ---------------------------------------
# Combine panels
library(patchwork)
p <- utr3_p + f0_p + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 9, face = "bold"))

# Export plot
library(Cairo)
CairoFonts(
  regular = "Arial:style=Regular",
  bold = "Arial:style=Black",
  italic = "Arial:style=Italic",
  bolditalic = "Arial:style=Black Italic",
  symbol = "Symbol"
)
cairo_pdf(filename = "utr3_statistics_6x3.pdf", family = "Arial", width = 6, height = 3) 
p
dev.off()
