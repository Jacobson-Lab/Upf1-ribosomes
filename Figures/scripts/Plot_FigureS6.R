# ---------------------------------------
# Figure S6 -- Metagene bin (4 bins) across normalized CDS for individual transcript
# ---------------------------------------

library(dplyr)
library(ggplot2)
library(ggh4x)

# Prepare data
Nterm_bin4 <- read.table(file = "../data/Data_FigureS6_Nterm_bin4.txt", header = TRUE, sep = "\t")
Nterm_bin4$strain <- factor(Nterm_bin4$strain, levels = c("WT + EV", "FLAG-UPF1"))
Nterm_bin4$bin <- factor(Nterm_bin4$bin, levels = c("0-25", "25-50", "50-75", "75-100"))
Nterm_bin4$size <- factor(Nterm_bin4$size, levels = c("S", "M", "L"))
Nterm_bin4$Ribosomes <- factor(Nterm_bin4$Ribosomes, levels = c("Total", "IP"))
Nterm_bin4$CHX <- as.factor(Nterm_bin4$CHX)
Nterm_bin4$strain_CHX <- paste0(Nterm_bin4$strain, "\n(", Nterm_bin4$CHX, ")")

Nterm_bin4_diff <- reshape2::dcast(Nterm_bin4, CHX + strain + strain_CHX + size + transcript + bin ~ Ribosomes, value.var = "value")
Nterm_bin4_diff$diff <- Nterm_bin4_diff$IP - Nterm_bin4_diff$Total
Nterm_bin4_diff <- Nterm_bin4_diff[complete.cases(Nterm_bin4_diff), ]

Cterm_bin4 <- read.table(file = "../data/Data_FigureS6_Cterm_bin4.txt", header = TRUE, sep = "\t")
Cterm_bin4$strain <- factor(Cterm_bin4$strain, levels = c("WT + EV", "UPF1-FLAG", "UPF1-FLAG/upf2Δ", "DE572AA-FLAG"))
Cterm_bin4$bin <- factor(Cterm_bin4$bin, levels = c("0-25", "25-50", "50-75", "75-100"))
Cterm_bin4$size <- factor(Cterm_bin4$size, levels = c("S", "M", "L"))
Cterm_bin4$Ribosomes <- factor(Cterm_bin4$Ribosomes, levels = c("Total", "IP"))
Cterm_bin4$CHX <- as.factor(Cterm_bin4$CHX)
Cterm_bin4$strain_CHX <- paste0(Cterm_bin4$strain, "\n(", Cterm_bin4$CHX, ")")

Cterm_bin4_diff <- reshape2::dcast(Cterm_bin4, CHX + strain + strain_CHX + size + transcript + bin ~ Ribosomes, value.var = "value")
Cterm_bin4_diff$diff <- Cterm_bin4_diff$IP - Cterm_bin4_diff$Total
Cterm_bin4_diff <- Cterm_bin4_diff[complete.cases(Cterm_bin4_diff), ]
Cterm_bin4_diff$strain_CHX <- factor(Cterm_bin4_diff$strain_CHX, levels = sort(unique(Cterm_bin4_diff$strain_CHX))[c(2, 3, 1)])

# Plot
A <- ggplot(Nterm_bin4[which(Nterm_bin4$CHX == "-CHX"), ]) +
  geom_boxplot(aes(x = bin, y = value*100, color = Ribosomes), outlier.size = 0.5, outlier.alpha = 0.5, size = 0.5,
               position = position_dodge2(preserve = "single")) + # keep width & position of boxplot consistent
  geom_hline(yintercept = 25, linetype = "dashed", color = "grey50") +
  facet_nested(size~strain+CHX, space = "free") +
  scale_color_manual(name = "Ribosomes", values = c(Total = "purple", IP = "orange"), limits = c("Total", "IP")) +
  xlab("% CDS") + ylab("% Footprint Count") +
  theme_bw(base_size = 8) + 
  theme(legend.position = "none", panel.grid = element_blank(), strip.background = element_blank(), 
        strip.text.x = element_text(face = "italic"), strip.text.y = element_text(angle = 0)) 

B <- ggplot(Nterm_bin4[which(Nterm_bin4$CHX == "+CHX"), ]) +
  geom_boxplot(aes(x = bin, y = value*100, color = Ribosomes), outlier.size = 0.5, outlier.alpha = 0.5, size = 0.5,
               position = position_dodge2(preserve = "single")) + # keep width & position of boxplot consistent
  geom_hline(yintercept = 25, linetype = "dashed", color = "grey50") +
  facet_nested(size~strain+CHX, space = "free") +
  scale_color_manual(name = "Ribosomes", values = c(Total = "purple", IP = "orange"), limits = c("Total", "IP")) +
  xlab("% CDS") + ylab("% Footprint Count") +
  theme_bw(base_size = 8) + 
  theme(legend.position = "none", panel.grid = element_blank(), strip.background = element_blank(), 
        strip.text.x = element_text(face = "italic"), strip.text.y = element_text(angle = 0)) 

C <- ggplot(Cterm_bin4) +
  geom_boxplot(aes(x = bin, y = value*100, color = Ribosomes), outlier.size = 0.5, outlier.alpha = 0.5, size = 0.5,
               position = position_dodge2(preserve = "single")) + # keep width & position of boxplot consistent
  geom_hline(yintercept = 25, linetype = "dashed", color = "grey50") +
  facet_nested(size~strain+CHX, space = "free") +
  scale_color_manual(name = "Ribosomes", values = c(Total = "purple", IP = "orange"), limits = c("Total", "IP")) +
  xlab("% CDS") + ylab("% Footprint Count") +
  theme_bw(base_size = 8) + 
  theme(legend.position = "right", panel.grid = element_blank(), strip.background = element_blank(), 
        strip.text.x = element_text(face = "italic"), strip.text.y = element_text(angle = 0)) 

D <- ggplot(Nterm_bin4_diff[which(Nterm_bin4_diff$CHX == "-CHX"), ]) +
  geom_boxplot(aes(x = bin, y = diff, color = size), outlier.size = 0.5, outlier.alpha = 0.5, size = 0.5,
               position = position_dodge2(preserve = "single")) + # keep width & position of boxplot consistent
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  facet_nested(size~strain_CHX, space = "free") +
  scale_color_manual(name = "Footprint size", values = c(S = "royalblue1", M = "tomato"), limits = c("S", "M")) +
  xlab("% CDS") + ylab("% Footprint Count\nIP - Total") +
  theme_bw(base_size = 8) + 
  theme(legend.position = "none", panel.grid = element_blank(), strip.background = element_blank(), 
        strip.text.x = element_text(face = "italic"), strip.text.y = element_text(angle = 0)) 

E <- ggplot(Nterm_bin4_diff[which(Nterm_bin4_diff$CHX == "+CHX"), ]) +
  geom_boxplot(aes(x = bin, y = diff, color = size), outlier.size = 0.5, outlier.alpha = 0.5, size = 0.5,
               position = position_dodge2(preserve = "single")) + # keep width & position of boxplot consistent
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  facet_nested(size~strain_CHX, space = "free") +
  scale_color_manual(name = "Footprint size", values = c(S = "royalblue1", M = "tomato"), limits = c("S", "M")) +
  xlab("% CDS") + ylab("% Footprint Count\nIP - Total") +
  theme_bw(base_size = 8) + 
  theme(legend.position = "right", panel.grid = element_blank(), strip.background = element_blank(), 
        strip.text.x = element_text(face = "italic"), strip.text.y = element_text(angle = 0)) 

FF <- ggplot(Cterm_bin4_diff) +
  geom_boxplot(aes(x = bin, y = diff, color = size), outlier.size = 0.5, outlier.alpha = 0.5, size = 0.5,
               position = position_dodge2(preserve = "single")) + # keep width & position of boxplot consistent
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  facet_nested(size~strain_CHX, space = "free") +
  scale_color_manual(name = "Footprint size", values = c(S = "royalblue1", M = "tomato"), limits = c("S", "M")) +
  xlab("% CDS") + ylab("% Footprint Count\nIP - Total") +
  theme_bw(base_size = 8) + 
  theme(legend.position = "none", panel.grid = element_blank(), strip.background = element_blank(), 
        strip.text.x = element_text(face = "italic"), strip.text.y = element_text(angle = 0)) 

# Combine panels
library(patchwork)
ABC <- ((A + B) / C) + plot_annotation(tag_levels = 'A')  & theme(plot.tag = element_text(size = 9, face = "bold"))
EE <- (E + plot_spacer()) + plot_layout(widths = c(0.5, 0.5))
EF <- (EE / FF) + plot_annotation(tag_levels = list(c('E', 'F')))  & theme(plot.tag = element_text(size = 9, face = "bold"))
DEF <- (D + EF + plot_layout(widths = c(1, 3))) + 
  plot_annotation(tag_levels = list(c('D', 'E', 'F'))) & theme(plot.tag = element_text(size = 9, face = "bold"))

p <- ABC / DEF

# Export plot
library(Cairo)
CairoFonts(
  regular = "Arial:style=Regular",
  bold = "Arial:style=Black",
  italic = "Arial:style=Italic",
  bolditalic = "Arial:style=Black Italic",
  symbol = "Symbol"
)
cairo_pdf(filename = "Metagene_bin4_8x11.pdf", family = "Arial", width = 8, height = 11) 
p
dev.off()
