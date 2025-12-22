library(data.table)
library(ggplot2)
library(readxl)
library(scales)
library(cowplot)
library(grid)

slratios <- fread("../data/derived_data/lema_sl_bins.csv")
slratios <- slratios[,.(id=sample_id, bin, centered_sl_ratio)]
bincoords <- fread("../data/reference/hg19_5Mb_bin_coordinates.tsv")
setkey(slratios, bin)
setkey(bincoords, bin)
slratios <- slratios[bincoords]

### LEMA meta - classifications
master <- read_excel("~/projects/lema/data/reference/lema_master_sheet.xlsx")
setDT(master)
#ids <- master[,c("oc", "seqID4x", "CNCLBL", "LCNSTAGE", "NTTBAS", "LCNCAT")]
ids <- master[,c("oc", "seqID4x", "CNCLBL", "LCNSTAGE", "LCNHIST", "LCNCAT")]

setnames(ids, c("patient", "id", "type", "stage", "subtype", "category"))
#ids <- ids[type=="Lung Cancer" & grepl("^II$|^III$|^IV$", stage) &
#           grepl("Adenocarcinoma|Squamous", subtype)]
ids <- ids[type=="Lung Cancer" & (grepl("Adenocarcinoma|Squamous", subtype) |
                                  grepl("^SCLC$", category))]
ids[category=="SCLC", subtype:="SCLC"]

########
armlevels <- c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
               "7p","7q","8p","8q", "9p", "9q","10p","10q","11p","11q","12p",
               "12q","13q","14q","15q","16p","16q","17p","17q","18p","18q",
               "19p", "19q","20p","20q","21q","22q")

### Plot S/L ratios
setkey(slratios, id)
setkey(ids, id)
slratios <- ids[slratios, nomatch=NULL]
# adjust labels to include N
labelData <- unique(ids[,.(id, subtype)])[,.N,by=.(subtype)]
labelData[,Label := paste0(subtype, '\n(n=', N, ')')]
setkey(labelData, subtype)
setkey(slratios, subtype)
slratios <- slratios[labelData]

################################################################################
############# prepare data for plotting ########################################

# calculate correlation for each sample and factor so that high correlation samples are plotted last
#fragmentRatios[,correlation := cor(centered_sl_ratio, median_normal_sl_ratio), by=.(seqID)]
#sampleOrder <- unique(fragmentRatios[,.(seqID, correlation)])[order(correlation)]$seqID
#fragmentRatios[,seqID := factor(seqID, levels=sampleOrder)]

# factor arms to be in a sensible order
slratios[,ChromArm := factor(ChromArm, levels=armlevels)]

# manually add some padding to a couple arms for plotting

## 10p, 12p, 16p, 16q, 17p, 18p, 19p, 19q, 20p, 20q, 21q, 22q
## 10p
slratios[bin == 312, "Center"] <- slratios[bin == 312, "Center"] + 10000000
## 12p
slratios[bin == 360, "Center"] <- slratios[bin == 360, "Center"] + 25000000
## 16p
slratios[bin == 435, "Center"] <- slratios[bin == 435, "Center"] + 10000000
## 16q
slratios[bin == 443, "Center"] <- slratios[bin == 443, "Center"] + 55000000
## 17p
slratios[bin == 447, "Center"] <- slratios[bin == 447, "Center"] + 25000000
## 18p
slratios[bin == 460, "Center"] <- slratios[bin == 460, "Center"] + 50000000
## 19p
slratios[bin == 475, "Center"] <- slratios[bin == 475, "Center"] + 30000000
## 19q
slratios[bin == 481, "Center"] <- slratios[bin == 481, "Center"] + 50000000
## 20p
slratios[bin == 486, "Center"] <- slratios[bin == 486, "Center"] + 50000000
## 20q
slratios[bin == 492, "Center"] <- slratios[bin == 492, "Center"] + 50000000
## 21q
slratios[bin == 498, "Center"] <- slratios[bin == 498, "Center"] + 50000000
## 22q
slratios[bin == 504, "Center"] <- slratios[bin == 504, "Center"] + 50000000

slratios <- slratios[,ChromArm:=factor(ChromArm, levels=armlevels)]
slratios <- slratios[,Label:=factor(Label, levels=c("Adenocarcinoma\n(n=466)",
                                                    "Squamous cell carcinoma\n(n=155)",
                                                    "SCLC\n(n=42)"))]
p1 <- ggplot() +
  geom_line(data=slratios, mapping=aes(x=Center, y=centered_sl_ratio, group=id), alpha=.5, color='#14478E', linewidth=.8) +
  facet_grid(Label ~ ChromArm, scales="free_x", space="free_x", switch="x") +
  scale_y_continuous(limits=c(-.45, .45), oob=squish, breaks=c(-.4, -.3, -.15, 0, .15, .3, .4)) +
  ylab("cfDNA Fragmentation Ratio") +
  theme_bw() +
  theme(axis.text.x=element_blank(), axis.text.y=element_text(size=24)) +
  theme(axis.ticks.x=element_blank(), axis.title.x=element_blank()) +
  theme(axis.title.y=element_text(size=32)) +
  theme(panel.grid.minor = element_blank()) +
  theme(strip.text.y=element_blank()) +  # Hide original y strip labels
  theme(strip.text.x.bottom=element_text(color="black", face="bold", size=19)) +
  theme(strip.text.x = element_text(angle=45)) +
  theme(strip.background = element_rect(fill="white", color="white")) +
  theme(panel.border=element_blank()) +
  theme(panel.grid=element_blank())

# Add labels on top using cowplot with better positioning and larger font
p_final <- ggdraw(p1) +
  # Main cancer type labels
  draw_label("Adenocarcinoma", x = 0.05, y = 0.97, hjust = 0, vjust = 0.5, size = 28, fontface = "bold") +
  draw_label("Squamous cell carcinoma", x = 0.05, y = 0.65, hjust = 0, vjust = 0.5, size = 28, fontface = "bold") +
  draw_label("SCLC", x = 0.05, y = 0.33, hjust = 0, vjust = 0.5, size = 28, fontface = "bold") +
  # Sample size labels underneath
  draw_label("(n=457)", x = 0.05, y = 0.95, hjust = 0, vjust = 0.5, size = 22) +
  draw_label("(n=149)", x = 0.05, y = 0.63, hjust = 0, vjust = 0.5, size = 22) +
  draw_label("(n=42)", x = 0.05, y = 0.31, hjust = 0, vjust = 0.5, size = 22)

ggsave("~/fig/LEMA_fragmentation_profiles.pdf", p_final, height=18, width=32)


###### ZSCORES
zscores <- fread("../data/derived_data/lema_zscores.csv")
zscores <- zscores[,.(id=sample_id, arm, z)]
zscores[,arm:=factor(arm, levels=paste0("z", sprintf("%02d", 1:39)), labels=armlevels)]

setkey(zscores, id)
setkey(ids, id)
zscores <- ids[zscores, nomatch=NULL]
zscores[,label:=ifelse(z < 0, "del", "dup")]
signed_log <- function(x) ifelse(x == 0, 0, sign(x) * sqrt(abs(x)))
zscores[,logz:=signed_log(z)]
zscores[,subtype:=ifelse(subtype=="Squamous cell carcinoma", "Squamous cell\ncarcinoma", subtype)]
zscores[,subtype:=factor(subtype, levels=c("Adenocarcinoma", "Squamous cell\ncarcinoma", "SCLC"))]

mytheme <-  theme(
  panel.grid=element_blank(),
  axis.text=element_text(size=10),
  axis.title=element_text(size=20),
  axis.text.y=element_blank(),
  axis.ticks.y=element_blank(),
  plot.caption = element_text(size = 30, face = "bold", hjust = 0.5),
  legend.position='none',
  panel.spacing = unit(0, "cm"),
  panel.spacing.x = unit(1, "cm"),
  strip.placement = "outside",
  strip.text.y=element_text(color="black", face="bold", size=15),
  strip.text.y.left = element_text(angle = 0),
  strip.text.x=element_text(color="black", face="bold", size=15),
  strip.background = element_rect(fill="white", color="white"))

g.zscores <- ggplot() +
        geom_point(data=zscores,
             mapping=aes(x=arm, y=z, fill=label, color=label)) +
    scale_fill_manual(values = c(
        "del" = "#76bebc",
        "dup" = "#ff8a47"
    )) +
    scale_color_manual(values = c(
        "del" = "#76bebc",
        "dup" = "#ff8a47"
    )) +
	geom_hline(yintercept=0, linetype="dashed", color="gray40") +
        #scale_y_continuous(breaks=c(-100, -25, 0, 25, 100), labels=c('-10', '-5', '0', '5', '10')) +
            facet_grid(arm ~ subtype, scale="free", switch="y")  +
            ylab('z-score') +
            xlab('Chromosome Arm') +
            coord_flip()
g.zscores <- g.zscores  + theme_classic() + mytheme
ggsave("~/fig/LEMA-zscores.pdf", g.zscores, height=9.5, width=7, units='in')

###### FLDS
flds <- fread("../data/derived_data/lema_flds.csv")
setnames(flds, "sample_id", "id")

setkey(flds, id)
setkey(ids, id)
flds <- ids[flds, nomatch=NULL]
flds[,subtype:=factor(subtype, levels=c("Adenocarcinoma", "Squamous cell carcinoma", "SCLC"), labels=c("Adenocarcinoma", "Squamous cell\ncarcinoma", "SCLC"))]
median.flds <- flds[,.(freq=mean(freq)), by=.(width, subtype)]
g.flds <- ggplot() +
	  geom_line(data=flds[width <= 230], aes(width, freq, group=id), color="gray80") +
	  geom_line(data=median.flds[width <= 230], aes(width, freq), color="skyblue4", linewidth=2) +
	  facet_wrap(subtype ~ ., ncol=1) +
	  scale_x_continuous(expand = c(0, 0)) +
	  scale_y_continuous(breaks=c(0.01, 0.02, 0.03, 0.04), expand = c(0, 0)) +
	  ylab("Frequency") + xlab("Fragment length") + theme_classic() + mytheme +
	 theme(strip.text.x = element_text(hjust=0, size=15), 
		   axis.text.y=element_text(size=15),
		   axis.text.x=element_text(size=15),
		   axis.ticks.y=element_line())
ggsave("~/fig/LEMA-flds.pdf", g.flds, height=9.5, width=5, units='in')
