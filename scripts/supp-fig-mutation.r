library(data.table)
library(ggplot2)
library(readxl)
master <- read_excel("../data/reference/LEMA_table_edited_github.xlsx")
setDT(master)
ids <- master[,c("seqID4x", "CNCLBL", "LCNSTAGE", "LCNHIST", "LCNCAT",
                 "EFGRBAS", "KRASBAS", "ROSBAS", "BRAFBAS", "ALKBAS", "ERBBBAS",
                 "RETBAS", "METBAS", "OMBAS")]
#ids <- master[,c("oc", "seqID4x", "CNCLBL", "LCNSTAGE", "LCNHIST", "NTTBAS")]

setnames(ids, c("id", "type", "stage", "subtype", "category", "EGFR",
                "KRAS", "ROS1", "BRAF", "ALK", "ERBB", "RET", "MET", "Other"))

### Monitoring score
dtf <- fread("../data/reference/monitoring_score.4x.v3.1.tsv")
setnames(dtf, "entity_id", "id")

setkey(dtf, id)
setkey(ids, id)
ids <- dtf[ids]

##### Mutants EGFR, KRAS, ALK, BRAF, ROS1, ERBB, RET, MET
ids[,egfr_mut:=fcase(grepl("wt", EGFR), "EGFR-",
                      !grepl("wt|NT|nt", EGFR), "EGFR+",
                      default = "NT")]
ids[,kras_mut:=fcase(grepl("wt", KRAS), "KRAS-",
                      !grepl("wt|NT|nt", KRAS), "KRAS+",
                      default = "NT")]
ids[,alk_mut:=fcase(grepl("wt", ALK), "ALK-",
                      !grepl("wt|NT|nt", ALK), "ALK+",
                      default = "NT")]
ids[,braf_mut:=fcase(grepl("wt", BRAF), "BRAF-",
                      !grepl("wt|NT|nt", BRAF), "BRAF+",
                      default = "NT")]
ids[,erbb_mut:=fcase(grepl("wt", ERBB), "ERBB-",
                      !grepl("wt|NT|nt", ERBB), "ERBB+",
                      default = "NT")]
ids[,ret_mut:=fcase(grepl("wt", RET), "RET-",
                      !grepl("wt|NT|nt", RET), "RET+",
                      default = "NT")]
ids[,met_mut:=fcase(grepl("wt", MET), "MET-",
                      !grepl("wt|NT|nt", MET), "MET+",
                      default = "NT")]

ids[,no_mt:=fcase(grepl("wt", EGFR) & grepl("wt", KRAS) &
                   grepl("wt", ROS1) & grepl("wt", BRAF) &
                   grepl("wt", ALK) & grepl("wt", ERBB) &
                   grepl("wt", RET) & grepl("wt", MET) & grepl("none", Other), "NM",
              !(grepl("wt|NT|nt", EGFR) & grepl("wt|NT|nt", KRAS) &
                   grepl("wt|NT|nt", ROS1) & grepl("wt|NT|nt", BRAF) &
                    grepl("wt|NT|nt", ALK) & grepl("wt|NT|nt", ERBB) & 
                    grepl("wt|NT|nt", RET) & grepl("wt|NT|nt", MET)),
                  "Mutant",
               default = "NT")]
#ids <- ids[no_mt != "NT"]


#subtype_colors <- c("LUAD" = "#E69F00", "LUSC" = "#56B4E9", "SCLC" = "#A83253")
ids <- ids[subtype=="Adenocarcinoma"]

ids[,egfr_mut:=factor(egfr_mut, levels=c("EGFR-", "EGFR+", "NT"))]
sample_sizes_egfr <- ids[order(egfr_mut), .(n = .N), by = egfr_mut][, label := paste0(egfr_mut, "\n(n=", n, ")")]
ids[,kras_mut:=factor(kras_mut, levels=c("KRAS-", "KRAS+", "NT"))]
sample_sizes_kras <- ids[order(kras_mut), .(n = .N), by = kras_mut][, label := paste0(kras_mut, "\n(n=", n, ")")]
ids[,erbb_mut:=factor(erbb_mut, levels=c("ERBB-", "ERBB+", "NT"))]
sample_sizes_erbb <- ids[order(erbb_mut), .(n = .N), by = erbb_mut][, label := paste0(erbb_mut, "\n(n=", n, ")")]
ids[,alk_mut:=factor(alk_mut, levels=c("ALK-", "ALK+", "NT"))]
sample_sizes_alk <- ids[order(alk_mut), .(n = .N), by = alk_mut][, label := paste0(alk_mut, "\n(n=", n, ")")]
ids[,braf_mut:=factor(braf_mut, levels=c("BRAF-", "BRAF+", "NT"))]
sample_sizes_braf <- ids[order(braf_mut), .(n = .N), by = braf_mut][, label := paste0(braf_mut, "\n(n=", n, ")")]
ids[,ret_mut:=factor(ret_mut, levels=c("RET-", "RET+", "NT"))]
sample_sizes_ret <- ids[order(ret_mut), .(n = .N), by = ret_mut][, label := paste0(ret_mut, "\n(n=", n, ")")]
ids[,met_mut:=factor(met_mut, levels=c("MET-", "MET+", "NT"))]
sample_sizes_met <- ids[order(met_mut), .(n = .N), by = met_mut][, label := paste0(met_mut, "\n(n=", n, ")")]

ids[,no_mt:=factor(no_mt, levels=c("NM", "Mutant", "NT"), labels=c("No mutation", "Mutant", "NT"))]
sample_sizes_nomt <- ids[order(no_mt), .(n = .N), by = no_mt][, label := paste0(no_mt, "\n(n=", n, ")")]

#### any mutations
g.mut <- ggplot(ids[no_mt!="NT"], aes(x = no_mt, y = monitoring_score)) +
  geom_boxplot(
    width = 0.6,
    fill = "gray90",
    color = "gray40",
    outlier.shape=NA,
    position = position_dodge(width = 0.75)
  ) +
  geom_jitter(
    color="steelblue3",
    size=1.5,
    width=0.1
  ) +
  scale_y_log10(
  breaks = c(0.001, 0.003, 0.01, 0.03,  0.1, 0.3),
  labels = c("0.1%", "0.3%", "1%", "3%", "10%", "30%")
) +
  labs(
    title = "Any mutation",
    y = "DELFI-TF"
  ) +
  scale_x_discrete(labels = setNames(sample_sizes_nomt$label, sample_sizes_nomt$no_mt)) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size=16),
    axis.title.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
    legend.position = "none"
  )
ggsave("~/fig/tmp.pdf", g.mut, height=5.5, width=3.5)

#### EGFR
g.egfr <- ggplot(ids[egfr_mut != "NT"], aes(x = egfr_mut, y = monitoring_score)) +
  geom_boxplot(
    width = 0.6,
    fill = "gray90",
    color = "gray40",
    outlier.shape=NA,
    position = position_dodge(width = 0.75)
  ) +
  geom_jitter(
    color="steelblue3",
    size=1.5,
    width=0.1
  ) +
  scale_y_log10(
  breaks = c(0.001, 0.003, 0.01, 0.03,  0.1, 0.3),
  labels = c("0.1%", "0.3%", "1%", "3%", "10%", "30%")
) +
  labs(
    title = "EGFR status",
    y = "DELFI-TF"
  ) +
  scale_x_discrete(labels = setNames(sample_sizes_egfr$label, sample_sizes_egfr$egfr)) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size=14),
    axis.title.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    legend.position = "none"
  )
ggsave("~/fig/tmp.pdf", g.egfr, height=3.5, width=2.5)


#### kras
g.kras <- ggplot(ids[kras_mut != "NT"], aes(x = kras_mut, y = monitoring_score)) +
  geom_boxplot(
    width = 0.6,
    fill = "gray90",
    color = "gray40",
    outlier.shape=NA,
    position = position_dodge(width = 0.75)
  ) +
  geom_jitter(
    color="steelblue3",
    size=1.5,
    width=0.1
  ) +
  scale_y_log10(
  breaks = c(0.001, 0.003, 0.01, 0.03,  0.1, 0.3),
  labels = c("0.1%", "0.3%", "1%", "3%", "10%", "30%")
) +
  labs(
    title = "KRAS status",
    y = "DELFI-TF"
  ) +
  scale_x_discrete(labels = setNames(sample_sizes_kras$label, sample_sizes_kras$kras)) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size=14),
    axis.title.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    legend.position = "none"
  )
ggsave("~/fig/tmp.pdf", g.kras, height=3.5, width=2.5)


#### alk
g.alk <- ggplot(ids[alk_mut != "NT"], aes(x = alk_mut, y = monitoring_score)) +
  geom_boxplot(
    width = 0.6,
    fill = "gray90",
    color = "gray40",
    outlier.shape=NA,
    position = position_dodge(width = 0.75)
  ) +
  geom_jitter(
    color="steelblue3",
    size=1.5,
    width=0.1
  ) +
  scale_y_log10(
  breaks = c(0.001, 0.003, 0.01, 0.03,  0.1, 0.3),
  labels = c("0.1%", "0.3%", "1%", "3%", "10%", "30%")
) +
  labs(
    title = "ALK status",
    y = "DELFI-TF"
  ) +
  scale_x_discrete(labels = setNames(sample_sizes_alk$label, sample_sizes_alk$alk)) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size=14),
    axis.title.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    legend.position = "none"
  )
ggsave("~/fig/tmp.pdf", g.alk, height=3.5, width=2.5)


#### braf
g.braf <- ggplot(ids[braf_mut != "NT"], aes(x = braf_mut, y = monitoring_score)) +
  geom_boxplot(
    width = 0.6,
    fill = "gray90",
    color = "gray40",
    outlier.shape=NA,
    position = position_dodge(width = 0.75)
  ) +
  geom_jitter(
    color="steelblue3",
    size=1.5,
    width=0.1
  ) +
  scale_y_log10(
  breaks = c(0.001, 0.003, 0.01, 0.03,  0.1, 0.3),
  labels = c("0.1%", "0.3%", "1%", "3%", "10%", "30%")
) +
  labs(
    title = "BRAF status",
    y = "DELFI-TF"
  ) +
  scale_x_discrete(labels = setNames(sample_sizes_braf$label, sample_sizes_braf$braf)) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size=14),
    axis.title.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    legend.position = "none"
  )
ggsave("~/fig/tmp.pdf", g.braf, height=3.5, width=2.5)


#### erbb
g.erbb <- ggplot(ids[erbb_mut != "NT"], aes(x = erbb_mut, y = monitoring_score)) +
  geom_boxplot(
    width = 0.6,
    fill = "gray90",
    color = "gray40",
    outlier.shape=NA,
    position = position_dodge(width = 0.75)
  ) +
  geom_jitter(
    color="steelblue3",
    size=1.5,
    width=0.1
  ) +
  scale_y_log10(
  breaks = c(0.001, 0.003, 0.01, 0.03,  0.1, 0.3),
  labels = c("0.1%", "0.3%", "1%", "3%", "10%", "30%")
) +
  labs(
    title = "ERBB status",
    y = "DELFI-TF"
  ) +
  scale_x_discrete(labels = setNames(sample_sizes_erbb$label, sample_sizes_erbb$erbb)) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size=14),
    axis.title.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    legend.position = "none"
  )
ggsave("~/fig/tmp.pdf", g.erbb, height=3.5, width=2.5)

#### ret
g.ret <- ggplot(ids[ret_mut != "NT"], aes(x = ret_mut, y = monitoring_score)) +
  geom_boxplot(
    width = 0.6,
    fill = "gray90",
    color = "gray40",
    outlier.shape=NA,
    position = position_dodge(width = 0.75)
  ) +
  geom_jitter(
    color="steelblue3",
    size=1.5,
    width=0.1
  ) +
  scale_y_log10(
  breaks = c(0.001, 0.003, 0.01, 0.03,  0.1, 0.3),
  labels = c("0.1%", "0.3%", "1%", "3%", "10%", "30%")
) +
  labs(
    title = "RET status",
    y = "DELFI-TF"
  ) +
  scale_x_discrete(labels = setNames(sample_sizes_ret$label, sample_sizes_ret$ret)) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size=14),
    axis.title.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    legend.position = "none"
  )
ggsave("~/fig/tmp.pdf", g.ret, height=3.5, width=2.5)

#### met
g.met <- ggplot(ids[met_mut != "NT"], aes(x = met_mut, y = monitoring_score)) +
  geom_boxplot(
    width = 0.6,
    fill = "gray90",
    color = "gray40",
    outlier.shape=NA,
    position = position_dodge(width = 0.75)
  ) +
  geom_jitter(
    color="steelblue3",
    size=1.5,
    width=0.1
  ) +
  scale_y_log10(
  breaks = c(0.001, 0.003, 0.01, 0.03,  0.1, 0.3),
  labels = c("0.1%", "0.3%", "1%", "3%", "10%", "30%")
) +
  labs(
    title = "MET status",
    y = "DELFI-TF"
  ) +
  scale_x_discrete(labels = setNames(sample_sizes_met$label, sample_sizes_met$met)) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size=14),
    axis.title.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    legend.position = "none"
  )
ggsave("~/fig/tmp.pdf", g.met, height=3.5, width=2.5)

####
ids[,wilcox.test(monitoring_score[egfr_mut=="EGFR+"],
                 monitoring_score[egfr_mut=="EGFR-"])]
ids[,wilcox.test(monitoring_score[kras_mut=="KRAS+"],
                 monitoring_score[kras_mut=="KRAS-"])]
ids[,wilcox.test(monitoring_score[erbb_mut=="ERBB+"],
                 monitoring_score[erbb_mut=="ERBB-"])]
ids[,wilcox.test(monitoring_score[braf_mut=="BRAF+"],
                 monitoring_score[braf_mut=="BRAF-"])]
ids[,wilcox.test(monitoring_score[alk_mut=="ALK+"],
                 monitoring_score[alk_mut=="ALK-"])]
ids[,wilcox.test(monitoring_score[ret_mut=="RET+"],
                 monitoring_score[ret_mut=="RET-"])]
ids[,wilcox.test(monitoring_score[met_mut=="MET+"],
                 monitoring_score[met_mut=="MET-"])]

####
ids[,wilcox.test(monitoring_score[no_mt=="Mutant"],
                 monitoring_score[no_mt=="No mutation"])]
