library(data.table)
library(ggplot2)

mytheme=theme(panel.grid.major = element_blank(),
			  panel.grid.minor = element_blank(),
			  axis.title=element_text(size=24),
			  axis.text=element_text(size=24),
			  strip.text.y=element_text(size=15),
			  plot.title=element_text(size=20),
			  #strip.text=element_text(color="white", face="bold", size=24),
			  #axis.ticks.x=element_blank(),
			  panel.spacing = unit(0,'lines'),
              strip.text=element_text(color="white", size=24),
              strip.background = element_rect(fill="#1B465A"))


lema.pcs <- fread("../data/derived_data/lema_pcs.csv")
lema.pcs <- lema.pcs[monitoring_score > 0.003]
clcgp.pcs <- fread("../data/derived_data/clcgp_pcs.csv")

mycolors <- c("#E69F00", "#56B4E9", "#A83253")
g <- ggplot(clcgp.pcs, aes(x = pc1, y = pc2, color = subtype, fill = subtype)) +
  geom_point(size = 3, alpha=0.8) +
  #geom_point(data=lema.pcs[subtype=="SCLC"], alpha = 0.6, size = 2) +
  stat_ellipse(level = 0.75, geom = "polygon", alpha = 0.15) +
  #stat_ellipse(level = 0.75, geom = "path", linewidth = 1) +
  scale_color_manual(values = c("LUAD" = "#E69F00",
                                 "LUSC" = "#56B4E9",
                                 "SCLC" = "#A83253")) +
  scale_fill_manual(values = c("LUAD" = "#E69F00",
                                "LUSC" = "#56B4E9",
                                "SCLC" = "#A83253")) +
  labs(title = "CLCGP",
       x = "Principal component 1",
       y = "Principal component 2",
       color = "Subtype",
       fill = "Subtype") +
  theme_classic() + 
  mytheme +
  theme(legend.text=element_text(size=18), 
        legend.title=element_text(size=24),
        legend.position="right",
        axis.title=element_text(size=34),
        axis.text=element_text(size=34),
        plot.title = element_text(size=36, hjust=0.5),
        plot.margin=margin(t=15, r=5, b=5, l=5, unit="pt"))
  ggsave("~/fig/clcgp-pca.pdf", g, width=12, height=12)

### LUAD projection
g <- ggplot(clcgp.pcs, aes(x = pc1, y = pc2, color = subtype, fill = subtype)) +
  geom_point(data=lema.pcs[subtype=="LUAD"], alpha=0.8, size = 2) +
  stat_ellipse(level = 0.75, geom = "polygon", alpha = 0.15) +
  #stat_ellipse(level = 0.75, geom = "path", linewidth = 1) +
  scale_color_manual(values = c("LUAD" = "#E69F00",
                                 "LUSC" = "#56B4E9",
                                 "SCLC" = "#A83253")) +
  scale_fill_manual(values = c("LUAD" = "#E69F00",
                                "LUSC" = "#56B4E9",
                                "SCLC" = "#A83253")) +
  labs(title = "LUAD",
       x = "Principal component 1",
       y = "Principal component 2",
       color = "Subtype",
       fill = "Subtype") +
  theme_classic() + 
  mytheme +
  theme(legend.text=element_text(size=18), 
        legend.title=element_text(size=24),
        legend.position="none",
        axis.title=element_text(size=30),
        axis.text=element_text(size=30),
        plot.title = element_text(size=30, hjust=0.5),
        plot.margin=margin(t=15, r=5, b=5, l=5, unit="pt"))
  ggsave("~/fig/lema-luad-pca.pdf", g, width=6, height=6)

### LUSC projection
g <- ggplot(clcgp.pcs, aes(x = pc1, y = pc2, color = subtype, fill = subtype)) +
  geom_point(data=lema.pcs[subtype=="LUSC"],  alpha=0.8, size = 2) +
  stat_ellipse(level = 0.75, geom = "polygon", alpha = 0.15) +
  #stat_ellipse(level = 0.75, geom = "path",  linewidth = 1) +
  scale_color_manual(values = c("LUAD" = "#E69F00",
                                 "LUSC" = "#56B4E9",
                                 "SCLC" = "#A83253")) +
  scale_fill_manual(values = c("LUAD" = "#E69F00",
                                "LUSC" = "#56B4E9",
                                "SCLC" = "#A83253")) +
  labs(title = "LUSC",
       x = "Principal component 1",
       y = "Principal component 2",
       color = "Subtype",
       fill = "Subtype") +
  theme_classic() + 
  mytheme +
  theme(legend.text=element_text(size=18), 
        legend.title=element_text(size=24),
        legend.position="none",
        axis.title=element_text(size=30),
        axis.text=element_text(size=30),
        plot.title = element_text(size=30, hjust=0.5),
        plot.margin=margin(t=15, r=5, b=5, l=5, unit="pt"))
  ggsave("~/fig/lema-lusc-pca.pdf", g, width=6, height=6)

### SCLC projection
g <- ggplot(clcgp.pcs, aes(x = pc1, y = pc2, color = subtype, fill = subtype)) +
  geom_point(data=lema.pcs[subtype=="SCLC"], alpha=0.8, size = 2) +
  stat_ellipse(level = 0.75, geom = "polygon", alpha = 0.15) +
  #stat_ellipse(level = 0.75, geom = "path", linewidth = 1) +
  scale_color_manual(values = c("LUAD" = "#E69F00",
                                 "LUSC" = "#56B4E9",
                                 "SCLC" = "#A83253")) +
  scale_fill_manual(values = c("LUAD" = "#E69F00",
                                "LUSC" = "#56B4E9",
                                "SCLC" = "#A83253")) +
  labs(title = "SCLC",
       x = "Principal component 1",
       y = "Principal component 2",
       color = "Subtype",
       fill = "Subtype") +
  theme_classic() + 
  mytheme +
  theme(legend.text=element_text(size=18), 
        legend.title=element_text(size=24),
        legend.position="none",
        axis.title=element_text(size=30),
        axis.text=element_text(size=30),
        plot.title = element_text(size=30, hjust=0.5),
        plot.margin=margin(t=15, r=5, b=5, l=5, unit="pt"))
  ggsave("~/fig/lema-sclc-pca.pdf", g, width=6, height=6)
