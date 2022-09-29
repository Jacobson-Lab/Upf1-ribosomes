# ---------------------------------------
# Figure 3 -- Footprint ends mapping relative to start or stop codons
# ---------------------------------------

library(ggplot2)
library(ggh4x)
library(scales)

df <- read.csv("../data/Data_Figure3.csv", header = TRUE)
df <- df[df$length %in% c(25:45), ]
df$strain2 <- paste0(df$strain, " (", df$CHX, ")")
df$strain2 <- factor(df$strain2, levels = sort(unique(df$strain2)))
df$end <- factor(df$end, levels = c("5' end", "3' end"))
size_labels <- data.frame(region = rep(c("Distance from start (nt)", "Distance from stop (nt)"), each = 2),
                          labs = c("M", "L", "M", "L"), ypos = c(29.5, 40, 29.5, 40), xpos = c(-32, -32, 32, 32))

mainp <- ggplot(df) +
  geom_tile(aes(x = dist, y = length, fill = fraction_average*100)) +
  geom_vline(xintercept = 0, linetype = 2, color = "red", size = 0.5) +
  geom_text(data = size_labels, aes(x = xpos, y = ypos, label = labs), size = 8/.pt, fontface = "bold") +
  facet_nested(strain2+end~region, switch = "x", scales = "free", nest_line = TRUE) +
  scale_fill_distiller(palette = "Spectral", type = "seq",
                       na.value = NA, trans = "log10",
                       breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  scale_x_continuous(expand = c(0.02, 0.02)) +
  xlab("") + ylab("Footprint length (nt)") +
  coord_cartesian(ylim = c(25, 45), clip = "off") +
  theme_bw(base_size = 10) + 
  theme(panel.grid = element_blank(), legend.position = "top", strip.text.y = element_text(face = "italic"), 
        strip.background.x = element_blank(), strip.background.y = element_blank(), strip.placement = "outside") +
  guides(fill = guide_colorbar(title = "% Footprint count    ", barheight = unit(0.3, "cm")))

size_fig <- data.frame(region = rep(c("Distance from start (nt)", "Distance from stop (nt)"), each = 3),
                       sizes = rep(c("", "L", "M"), times = 2),
                       x = c(-30, -26, -11, -50, -26, -11), xend = c(50, 14, 14, 30, 14, 14), 
                       y = rep(c(1.2, 1, 0), times = 2), yend = rep(c(1.2, 1, 0), times = 2),
                       left_labs = rep(c("", "5'", "5'"), times = 2), right_labs = rep(c("", "3'   L", "3'   M"), times = 2))
size_fig$left_pos <- size_fig$x - 2
size_fig$right_pos <- size_fig$xend + 2
                       
sp <- ggplot(size_fig) +
  geom_segment(aes(x = x, xend = xend, y = y, yend = yend, color = sizes)) +
  geom_text(aes(x = left_pos, y = y, label = left_labs), size = 10/.pt, fontface = "bold", hjust = 1) +
  geom_text(aes(x = right_pos, y = y, label = right_labs), size = 10/.pt, fontface = "bold", hjust = 0) +
  facet_grid(.~region, scales = "free") +
  scale_x_continuous(expand = c(0.02, 0.02)) +
  scale_color_manual(values = c("white", "black", "black")) +
  coord_cartesian(clip = "off") +
  theme_void(base_size = 10) +  
  theme(strip.placement = "outside", legend.position = "none", strip.text = element_blank())

# Combine panels
library(patchwork)
p <- (mainp / sp) + 
  plot_layout(heights = c(9.5, 0.5))

# Export plot
library(Cairo)
CairoFonts(
  regular = "Arial:style=Regular",
  bold = "Arial:style=Black",
  italic = "Arial:style=Italic",
  bolditalic = "Arial:style=Black Italic",
  symbol = "Symbol"
)
cairo_pdf(filename = "Footprint_ends_mapping_6x6.pdf", family = "Arial", width = 6, height = 6) 
p
dev.off()
