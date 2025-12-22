library(data.table)
library(ggplot2)
library(readxl)
library(EnvStats)
library(pROC)
library(MASS)

armlevels <- c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
               "7p","7q","8p","8q", "9p", "9q","10p","10q","11p","11q","12p",
               "12q","13q","14q","15q","16p","16q","17p","17q","18p","18q",
               "19p", "19q","20p","20q","21q","22q")

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


### LEMA metadata
master <- read_excel("../data/reference/lema_master_sheet.xlsx")
setDT(master)
ids <- master[,c("oc", "seqID4x", "CNCLBL", "LCNSTAGE", "LCNHIST", "LCNCAT")]

setnames(ids, c("patient", "id", "type", "stage", "subtype", "category"))
ids <- ids[type=="Lung Cancer" & (grepl("Adenocarcinoma|Squamous", subtype) |
                                  grepl("^SCLC$", category))]

ids[category=="SCLC", subtype:="SCLC"]

### Monitoring score
dtf <- fread("../data/reference/monitoring_score.4x.v3.1.tsv")
setnames(dtf, "entity_id", "id")

setkey(dtf, id)
setkey(ids, id)
ids <- ids[dtf, nomatch=NULL]
#ids <- ids[monitoring_score > 0.003]

### CLCGP data
clcgp.dt <- fread("../data/reference/clcgp_1mb_segments.csv")
binfilter <- fread("../data/reference/bin-cov-filter_1mb.csv")
setkey(binfilter, chr, start, end, bin)
setkey(clcgp.dt, chr, start, end, bin)

clcgp.dt <- clcgp.dt[binfilter, nomatch=NULL]
clcgp.dt[,signal:=scale(signal)]

### matrix for PCA
mat <- dcast(clcgp.dt, Sample + histology ~ bin, value.var="signal")
pca <- prcomp(mat[,-c(1,2)], center=TRUE)
pc.dt <- data.table(subtype=mat$histology, pc1=pca$x[,1], pc2=pca$x[,2],
                    pc3=pca$x[,3], pc4=pca$x[,4], pc5=pca$x[,5], pc6=pca$x[,6],
                    pc7=pca$x[,7], pc8=pca$x[,8])
pc.dt[,subtype:=factor(subtype, levels=c("AD", "SQ", "SCLC"), labels=c("LUAD", "LUSC", "SCLC"))]

#### PCA of CGLCP 
g <- ggplot(pc.dt, aes(x=pc1, y=pc2, color=subtype))  + geom_point(size=1.5) + theme_classic(base_size=15)
g <- g + ggtitle("PCA of CLCGP")
ggsave("~/fig/tmp.pdf", g, width=8)

### LEMA data
bins <- fread("../data/derived_data/lema-denoised-1mb.csv")
bins[,signal:=scale(svd.foreground), by=id]

lema.mat <- dcast(bins, id ~ bin, value.var="signal")
lema.proj <- predict(pca, newdata=lema.mat[,-1])
xx <- data.table(id=lema.mat$id, pc1=lema.proj[,1], pc2=lema.proj[,2],
                    pc3=lema.proj[,3], pc4=lema.proj[,4], pc5=lema.proj[,5],
                    pc6=lema.proj[,6], pc7=lema.proj[,7], pc8=lema.proj[,8])

setkey(ids, id)
setkey(xx, id)
xx <- xx[ids, nomatch=NULL]
xx[,subtype:=factor(subtype, levels=c("Adenocarcinoma", "Squamous cell carcinoma", "SCLC"), labels=c("LUAD", "LUSC", "SCLC"))]
fwrite(xx, "../data/derived_data/lema_pcs.csv")
fwrite(pc.dt, "../data/derived_data/clcgp_pcs.csv")

xx <- xx[monitoring_score > 0.003]

#### train model on TCGA
#model <- glm(subtype ~ pc1 + pc2 + pc3 + pc4 + pc5, data=pc.dt, family="binomial")
pc.dt[,NSCLC:=ifelse(subtype!="SCLC", "NSCLC", "SCLC")]
pc.dt[,NSCLC:=factor(NSCLC, levels=c("NSCLC", "SCLC"))]

model1 <- glm(NSCLC ~ pc1 + pc2 + pc3, data=pc.dt, family="binomial")
model2 <- glm(subtype ~ pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8, data=pc.dt[subtype %in% c("LUAD", "LUSC")], family="binomial")

### Predict on LEMA
xx[,NSCLC:=ifelse(subtype!="SCLC", "NSCLC", "SCLC")]
xx[,NSCLC:=factor(NSCLC, levels=c("NSCLC", "SCLC"))]

#xx2 <- xx[subtype %in% c("ADC", "SCC")]
xx[,predictions1:=predict(model1, newdata=xx)]
xx[,predictions2:=predict(model2, newdata=xx)]
fwrite(xx, "../data/derived_data/lema_preds.csv")


### Old PCA plot
mycolors <- c("#E69F00", "#56B4E9", "#A83253")
## PC1 vs PC2
g.pc1 <- ggplot(pc.dt, aes(x=pc1, y=pc2, color=subtype))  + 
  geom_point(alpha=0.5) + 
  theme_classic() + 
  mytheme +
  theme(legend.text=element_text(size=18), 
        legend.title=element_text(size=24),
        legend.position="none",
        plot.margin=margin(t=15, r=5, b=5, l=5, unit="pt")) +  # Increase left margin
  geom_point(data=xx, aes(x=pc1, y=pc2, color=subtype), size=2) + 
  scale_colour_manual(values=mycolors) +
  xlab("Principal Component 1\n(13.1% variance explained)") + 
  ylab("Principal Component 2\n(5.4% variance explained)")

ggsave("~/fig/tmp.pdf", g.pc1, width=9, height=8)

# Add a combined grouping variable
pc.dt$legend_group <- paste("CLCGP", pc.dt$subtype, sep=" - ")
xx$legend_group <- paste("LEMA", xx$subtype, sep=" - ")

# Combine datasets
combined <- rbind(pc.dt[,.(pc1, pc2, pc3, pc4, pc5, subtype, legend_group)], xx[,.(pc1, pc2, pc3, pc4, pc5, subtype, legend_group)])

# Create the plot
g.pc2 <- ggplot(combined, aes(x=pc2, y=pc3, color=legend_group)) + 
  geom_point(aes(alpha=legend_group), size=ifelse(combined$legend_group %in% paste("LEMA", c("LUAD", "LUSC", "SCLC"), sep=" - "), 2, 2)) +
  scale_colour_manual(values=c(
    "LEMA - LUAD" = mycolors[1],
    "LEMA - LUSC" = mycolors[2], 
    "LEMA - SCLC" = mycolors[3],
    "CLCGP - LUAD" = mycolors[1],
    "CLCGP - LUSC" = mycolors[2],
    "CLCGP - SCLC" = mycolors[3]
  )) +
  scale_alpha_manual(values=c(
    "LEMA - LUAD" = 1, "LEMA - LUSC" = 1, "LEMA - SCLC" = 1,
    "CLCGP - LUAD" = 0.2, "CLCGP - LUSC" = 0.2, "CLCGP - SCLC" = 0.2
  )) +
  guides(color=guide_legend(title="", override.aes=list(alpha=c(0.2, 0.2, 0.2, 1, 1, 1))), alpha="none") +
  theme_classic() + 
  mytheme +
  theme(legend.text=element_text(size=18),
        plot.margin=margin(t=15, r=5, b=5, l=5, unit="pt"),
        legend.title=element_text(size=20)) +
  xlab("Principal Component 2\n(5.4% variance explained)") + 
  ylab("Principal Component 3\n(4.8% variance explained)")

ggsave("~/fig/tmp.pdf", g.pc2, width=11, height=8)
