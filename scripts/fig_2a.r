library(data.table)
library(readxl)
library(ggplot2)

armlevels <- c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
               "7p","7q","8p","8q", "9p", "9q","10p","10q","11p","11q","12p",
               "12q","13q","14q","15q","16p","16q","17p","17q","18p","18q",
               "19p", "19q","20p","20q","21q","22q")

### CLCGP bins
clcgp.dt <- fread("../data/reference/clcgp_1mb_segments.csv")
binfilter <- fread("../data/reference/bin-cov-filter_1mb.csv")
setkey(binfilter, chr, start, end, bin)
setkey(clcgp.dt, chr, start, end, bin)
clcgp.dt <- clcgp.dt[binfilter, nomatch=NULL]
clcgp.dt[,signal:=scale(signal)]
clcgp.dt[,subtype:=factor(histology, levels=c("AD", "SQ", "SCLC"),
                          labels=c("LUAD", "LUSC", "SCLC"))]

#clcgp.mean <- clcgp.dt[,.(signal=mean(signal)), by=.(chr, start, end, arm, bin, subtype)]

clcgp.dt[,signal2:=runmed(signal, 13), by=.(Sample, arm)]
clcgp.mean <- clcgp.dt[,.(prop_dup=sum(signal2 > 0.5)/.N, prop_del=sum(signal2 < -0.5)/.N), by=.(chr, start, end, arm, bin, subtype)]

### load plasma bins
lema.dt <- fread("../data/derived_data/lema-denoised-1mb.csv")
lema.dt[,signal:=scale(svd.foreground), by=id]

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
ids <- ids[monitoring_score > 0.003]
#ids <- ids[monitoring_score > 0.002]
setkey(ids, id)
setkey(lema.dt, id)
lema.dt <- lema.dt[ids, nomatch=NULL]
lema.dt[,signal2:=runmed(signal, 13), by=.(id, arm)]

#lema.mean <- lema.dt[,.(signal=mean(signal)), by=.(chr, start, end, arm, bin, subtype)]
lema.mean <- lema.dt[,.(prop_dup=sum(signal2 > 0.5)/.N, prop_del=sum(signal2 < -0.5)/.N), by=.(chr, start, end, arm, bin, subtype)]
lema.mean[,subtype:=factor(subtype, levels=c("Adenocarcinoma", "Squamous cell carcinoma", "SCLC"),
                           labels=c("LUAD", "LUSC", "SCLC"))]


clcgp.mean[,source:="CLCGP"]
lema.mean[,source:="LEMA"]

binned.means <- rbind(clcgp.mean, lema.mean)
binned.means[subtype=="LUAD" & bin==15]
#binned.means[,signal.centered:=signal-mean(signal), by=.(source, subtype)]
binned.means[,arm:=factor(arm, levels=armlevels)]
binned.means[order(arm), bin2:=1:.N, by=.(source, subtype)]
binned.means[,source:=factor(source, levels=c("CLCGP", "LEMA"))]

#binned.means[,.(mean(signal)), by=.(arm, source, subtype)][order(subtype, arm)][subtype=="SCLC"]

# Calculate sample sizes for each source
lema.sample.size <- lema.dt[bin==2, .N, by=subtype][, source:="LEMA"]
lema.sample.size[,subtype:=factor(subtype, levels=c("Adenocarcinoma", "Squamous cell carcinoma", "SCLC"),
                           labels=c("LUAD", "LUSC", "SCLC"))]

clcgp.sample.size <- clcgp.dt[bin==2, .N, by=subtype][, source:="CLCGP"]

source.sample.size <- rbind(lema.sample.size, clcgp.sample.size)
# Create named vector for source labels
source.sample.size[,label:=paste0(source, "\nn=",N)]

luad.sample.size <- source.sample.size[subtype=="LUAD", .(source, label)]
luad.labels <- setNames(
  luad.sample.size$label,
  luad.sample.size$source
)

lusc.sample.size <- source.sample.size[subtype=="LUSC", .(source, label)]
lusc.labels <- setNames(
  lusc.sample.size$label,
  lusc.sample.size$source
)

sclc.sample.size <- source.sample.size[subtype=="SCLC", .(source, label)]
sclc.labels <- setNames(
  sclc.sample.size$label,
  sclc.sample.size$source
)

mytheme <-  theme(
  panel.grid=element_blank(),
  axis.text=element_text(size=18),
  axis.title=element_text(size=20),
  axis.text.y=element_blank(),
  axis.ticks.y=element_blank(),
  plot.caption = element_text(size = 30, face = "bold", hjust = 0.5),
  legend.position='none',
  panel.spacing = unit(0, "cm"),
  strip.placement = "outside",
  strip.text.y=element_text(color="black", face="bold", size=15),
  strip.text.y.left = element_text(angle = 0),
  strip.text.x=element_text(color="black", face="bold", size=20),
  strip.background = element_rect(fill="white", color="white"))

g.luad <- ggplot() +
        geom_bar(data=binned.means[subtype=="LUAD"],
             mapping=aes(x=bin2, y=prop_dup, fill=paste0(source, "_dup"), color=paste0(source, "_dup") ),
             stat='identity', width=1) +
    geom_bar(data=binned.means[subtype=="LUAD"],
             mapping=aes(x=bin2, y=-prop_del, fill=paste0(source, "_del"), color=paste0(source, "_del")),
             stat='identity', width=1) +
    scale_fill_manual(values = c(
        "LEMA_del" = "#76bebc",
        "LEMA_dup" = "#ff8a47",
        "CLCGP_del" = "#769cbe",
        "CLCGP_dup" = "#c75252"
    )) +
    scale_color_manual(values = c(
        "LEMA_del" = "#76bebc",
        "LEMA_dup" = "#ff8a47",
        "CLCGP_del" = "#769cbe",
        "CLCGP_dup" = "#c75252"
    )) +
        scale_y_continuous(breaks=c(-.8, -.4, 0, .4, .8), labels=c('80%', '40%', '0%', '40%', '80%'), limits=c(-1, 1)) +
            facet_grid(arm ~ source, scale="free_y", switch="y", labeller = labeller(source = luad.labels))  +
            xlab('Chromosome Arm') +
            ylab('Percentage of Cases\nWith CNA') +
            coord_flip()
g.luad <- g.luad  + labs(caption="Adenocarcinoma") + theme_classic() + mytheme
ggsave("~/fig/lema-luad.pdf", g.luad, height=15, width=8, units='in')

g.lusc <- ggplot() +
        geom_bar(data=binned.means[subtype=="LUSC"],
             mapping=aes(x=bin2, y=prop_dup, fill=paste0(source, "_dup"), color=paste0(source, "_dup") ),
             stat='identity', width=1) +
    geom_bar(data=binned.means[subtype=="LUSC"],
             mapping=aes(x=bin2, y=-prop_del, fill=paste0(source, "_del"), color=paste0(source, "_del")),
             stat='identity', width=1) +
    scale_fill_manual(values = c(
        "LEMA_del" = "#76bebc",
        "LEMA_dup" = "#ff8a47",
        "CLCGP_del" = "#769cbe",
        "CLCGP_dup" = "#c75252"
    )) +
    scale_color_manual(values = c(
        "LEMA_del" = "#76bebc",
        "LEMA_dup" = "#ff8a47",
        "CLCGP_del" = "#769cbe",
        "CLCGP_dup" = "#c75252"
    )) +
        scale_y_continuous(breaks=c(-.8, -.4, 0, .4, .8), labels=c('80%', '40%', '0%', '40%', '80%'), limits=c(-1, 1)) +
            facet_grid(arm ~ source, scale="free_y", switch="y", labeller = labeller(source = lusc.labels))  +
            xlab('Chromosome Arm') +
            ylab('Percentage of Cases\nWith CNA') +
            coord_flip()
g.lusc <- g.lusc  + labs(caption="Squamous Cell Carcinoma") + theme_classic() + mytheme
ggsave("~/fig/lema-lusc.pdf", g.lusc, height=15, width=8, units='in')

#####
g.sclc <- ggplot() +
        geom_bar(data=binned.means[subtype=="SCLC"],
             mapping=aes(x=bin2, y=prop_dup, fill=paste0(source, "_dup"), color=paste0(source, "_dup") ),
             stat='identity', width=1) +
    geom_bar(data=binned.means[subtype=="SCLC"],
             mapping=aes(x=bin2, y=-prop_del, fill=paste0(source, "_del"), color=paste0(source, "_del")),
             stat='identity', width=1) +
    scale_fill_manual(values = c(
        "LEMA_del" = "#76bebc",
        "LEMA_dup" = "#ff8a47",
        "CLCGP_del" = "#769cbe",
        "CLCGP_dup" = "#c75252"
    )) +
    scale_color_manual(values = c(
        "LEMA_del" = "#76bebc",
        "LEMA_dup" = "#ff8a47",
        "CLCGP_del" = "#769cbe",
        "CLCGP_dup" = "#c75252"
    )) +
        scale_y_continuous(breaks=c(-.8, -.4, 0, .4, .8), labels=c('80%', '40%', '0%', '40%', '80%'), limits=c(-1, 1)) +
            facet_grid(arm ~ source, scale="free_y", switch="y", labeller = labeller(source = sclc.labels))  +
            xlab('Chromosome Arm') +
            ylab('Percentage of Cases\nWith CNA') +
            coord_flip()
g.sclc <- g.sclc  + labs(caption="Small Cell Carcinoma") + theme_classic() + mytheme
ggsave("~/fig/lema-sclc.pdf", g.sclc, height=15, width=8, units='in')


##### R^2s (for results)
types <- c(rep(1, 2252), rep(-1, 2252))
cn.sclc.lema <- binned.means[subtype=="SCLC" & source=="LEMA"][order(bin2)][, c(prop_del, prop_dup)]
cn.sclc.clcgp <- binned.means[subtype=="SCLC" & source=="CLCGP"][order(bin2)][, c(prop_del, prop_dup)]
cor(cn.sclc.lema, cn.sclc.clcgp)

summary(lm(cn.sclc.lema ~ cn.sclc.clcgp + types))
summary(lm(cn.sclc.lema ~ cn.luad.clcgp + types))
summary(lm(cn.sclc.lema ~ cn.lusc.clcgp + types))


cn.luad.lema <- binned.means[subtype=="LUAD" & source=="LEMA"][order(bin2)][, c(prop_del, prop_dup)]
cn.luad.clcgp <- binned.means[subtype=="LUAD" & source=="CLCGP"][order(bin2)][, c(prop_del, prop_dup)]
cor(cn.luad.lema, cn.luad.clcgp)

summary(lm(cn.luad.lema ~ cn.luad.clcgp + types))
summary(lm(cn.luad.lema ~ cn.sclc.clcgp + types))
summary(lm(cn.luad.lema ~ cn.lusc.clcgp + types))

cn.lusc.lema <- binned.means[subtype=="LUSC" & source=="LEMA"][order(bin2)][, c(prop_del, prop_dup)]
cn.lusc.clcgp <- binned.means[subtype=="LUSC" & source=="CLCGP"][order(bin2)][, c(prop_del, prop_dup)]
cor(cn.lusc.lema, cn.lusc.clcgp)

summary(lm(cn.lusc.lema ~ cn.lusc.clcgp + types))
summary(lm(cn.lusc.lema ~ cn.sclc.clcgp + types))
summary(lm(cn.lusc.lema ~ cn.luad.clcgp + types))
