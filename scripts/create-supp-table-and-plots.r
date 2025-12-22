library(data.table)
library(readxl)
library(ggplot2)

master <- read_excel("~/projects/lema/data/reference/LEMA_table_edited_github.xlsx")
setDT(master)

### samples that failed QC
filterids <- c("DL004703NSP0", "DL004710NSP0", "DL004885NSP0", "DL005317NSP0", "DL004755NSP0",
               "DL004956NSP0", "DL005257NSP0", "DL004678NSP0", "DL005296NSP0")
filterstr <- paste0(filterids, collapse="|")

#### Parsing storage time column
master[,Date_stored_freezerNKI := as.Date(Date_stored_freezerNKI, format='%Y-%m-%d')]
master[,fixed_Date_gathered := as.Date(fixed_Date_gathered, format='%Y-%m-%d')]
master[,storageTime := difftime(Date_stored_freezerNKI, fixed_Date_gathered, units='days')]#                                  grepl("^SCLC$", category))]
# there are some obvious errors , na these
master[storageTime <= -1, storageTime := NA]

# bin the storage times into 1, 2, 3+ days
master[storageTime == 1, storageTimeBin := '1 Day']
master[storageTime == 2, storageTimeBin := '2 Day']
master[storageTime >= 3, storageTimeBin := '3+ Days']

#ids <- master[,c("seqID4x", "CNCLBL", "LCNSTAGE", "NTTBAS", "LCNCAT")]
### Need diabetes, CVD, COPD, autoimmune disorders
ids <- master[,.(id=seqID4x, type=CNCLBL, stage=LCNSTAGE,
                 subtype=LCNHIST, category=LCNCAT, age=AGE, sex=SEX, race=RACE,
                 bmi=BMI, packyears=SUPACKYR, smoking_status=SUCAT,
                 prior_cancer=PRCNRHS, number_comorbidities=NUM_COMORB,
                 COPD=`COBAS#Chronic_lung_disease`, diabetes=`COBAS#Diabetes`,
                 CVD=`COBAS#Cardiovascular_disease`, autoimmune=`COBAS#Autoimmune_disease`,
                 antipyretic=`MEDBAS#Antipyretic`, opioid=`MEDBAS#Opioids`,
                 steroid=`MEDBAS#Systemic_steroids`, nsaid=`MEDBAS#NSAIDS`,
                 anticoagulant=`MEDBAS#Anticoagulants`, no_medication=`MEDBAS#None`,
                 T_stage=TNMBAS_cTNM_T, N_stage=TNMBAS_cTNM_N,
				 storage_time=storageTimeBin, streck=PROCESSED_STRECK)]

ids[, qc_fail := ifelse(is.na(id) | grepl(filterstr, id), "FAIL", "PASS")]
#ids <- ids[type=="Lung Cancer" & grepl("^II$|^III$|^IV$", stage) &
#           grepl("Adenocarcinoma|Squamous", subtype)]
#ids <- ids[type=="Lung Cancer" & (grepl("Adenocarcinoma|Squamous", subtype) |
# calculate time difference

ids[category=="SCLC", subtype:="SCLC"]
### Monitoring score
dtf <- fread("~/projects/lema/data/reference/monitoring_score.4x.v3.1.tsv")
setnames(dtf, "entity_id", "id")

setkey(dtf, id)
setkey(ids, id)
dt <- dtf[ids]
fwrite(dt, "LEMA_table.csv")
ids <- ids[type=="Lung Cancer"]
ids <- ids[!is.na(subtype)]

dt[,table(stage)]

dt2 <- dt[grepl("Adenocarcinoma|Squamous cell carcinoma|SCLC", subtype)]
dt2[,median(monitoring_score, na.rm=TRUE), by=subtype]

library(broom)
model <- dt2[,lm(log(monitoring_score) ~  subtype + stage)]
tidy_output <- tidy(model)
write.csv(tidy_output, "regression_results.csv", row.names = FALSE)

dt2[stage > "I",cor(monitoring_score, monitoring_score, use="complete.obs", method="spearman"), by=subtype]

dt[,wilcox.test(monitoring_score[smoking_status=="Former"], monitoring_score[smoking_status=="Current"])]
dt[,wilcox.test(monitoring_score[smoking_status=="Never"], monitoring_score[smoking_status=="Current"])]
dt[,wilcox.test(monitoring_score[smoking_status=="Never"], monitoring_score[smoking_status=="Former"])]
dt[,wilcox.test(monitoring_score[CVD==0], monitoring_score[CVD==1])]
dt[,wilcox.test(monitoring_score[autoimmune==0], monitoring_score[autoimmune==1])]
dt[,wilcox.test(monitoring_score[COPD==0], monitoring_score[COPD==1])]
dt[,wilcox.test(monitoring_score[diabetes==0], monitoring_score[diabetes==1])]
dt[,wilcox.test(monitoring_score[no_medication==0], monitoring_score[no_medication==1])]
dt[,wilcox.test(monsteroidng_score[no_medication==0], monitoring_score[no_medication==1])]
dt[,wilcox.test(monitoring_score[antipyretic==1], monitoring_score[no_medication==1])]
dt[,wilcox.test(monitoring_score[anticoagulant==1], monitoring_score[no_medication==1])]
dt[,wilcox.test(monitoring_score[nsaid==1], monitoring_score[no_medication==1])]
dt[,wilcox.test(monitoring_score[steroid==1], monitoring_score[no_medication==1])]
dt[,summary(lm(log(monitoring_score) ~ smoking_status))]




#### TUMOR STAGING PLOTS
dt[is.na(N_stage), N_stage:="Not reported"]
dt[is.na(T_stage), T_stage:="Not reported"]
dt[,N_stage:=factor(N_stage, levels=c("N0", "N1", "N2", "N3", "Nx", "Not reported"))]

sample_sizes_nstage <- dt[, .(n = .N), by = N_stage][, label := paste0(N_stage, "\n(n=", n, ")")]
#### N-stage plot
g.nstage <- ggplot(dt, aes(x = N_stage, y = monitoring_score)) +
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
    title = "Node Stage",
    y = "DELFI-TF"
  ) +
  scale_x_discrete(labels = setNames(sample_sizes_nstage$label, sample_sizes_nstage$N_stage)) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size=14),
    axis.title.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    legend.position = "none"
  )
ggsave("~/fig/tmp.pdf", g.nstage, height=3.5, width=6)

dt[,T_stage:=gsub("a|b|c|mi", "", T_stage)]
dt[,T_stage:=factor(T_stage, levels=c("T0", "T1", "T2", "T3", "T4", "Tx", "Not reported"))]

sample_sizes_tstage <- dt[order(T_stage), .(n = .N), by = T_stage][, label := paste0(T_stage, "\n(n=", n, ")")]

#### T-stage plot
g.tstage <- ggplot(dt, aes(x = T_stage, y = monitoring_score)) +
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
    title = "Tumor Stage",
    y = "DELFI-TF"
  ) +
  scale_x_discrete(labels = setNames(sample_sizes_tstage$label, sample_sizes_tstage$T_stage)) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size=14),
    axis.title.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    legend.position = "none"
  )
ggsave("~/fig/tmp.pdf", g.tstage, height=3.5, width=6)

####
dt[,TNM_stage:=paste0(T_stage, ":", N_stage)]
tnm_levels <- dt[order(T_stage, N_stage)][T_stage != "Tx" & T_stage != "T0" & N_stage!="Nx" & T_stage!="Not reported", unique(TNM_stage)]
dt.tnm <- dt[T_stage != "Tx" & T_stage != "T0" & N_stage!="Nx" & T_stage!="Not reported"]
dt.tnm[,TNM_stage:=factor(TNM_stage, levels=tnm_levels)]
sample_sizes_tnmstage <- dt.tnm[order(TNM_stage), .(n = .N), by = TNM_stage][, label := paste0(TNM_stage, "\n(n=", n, ")")]

colors <- c(
  "T1" = "#E63946",
  "T2" = "#F77F00",
  "T3" = "#06AED5",
  "T4" = "#8AC926"
)

#### TNM-stage plot
g.tnmstage <- ggplot(dt.tnm, aes(x = TNM_stage, y = monitoring_score)) +
  geom_boxplot(
    width = 0.6,
    fill = "gray90",
    color = "gray40",
    outlier.shape=NA,
    position = position_dodge(width = 0.75)
  ) +
  geom_jitter(
    aes(color=T_stage),
    size=1.5,
    width=0.1
  ) +
  scale_y_log10(
  breaks = c(0.001, 0.003, 0.01, 0.03,  0.1, 0.3),
  labels = c("0.1%", "0.3%", "1%", "3%", "10%", "30%")
  ) +
  labs(
    title = "TNM Stage",
    y = "DELFI-TF"
  ) +
  scale_x_discrete(labels = setNames(sample_sizes_tnmstage$label, sample_sizes_tnmstage$TNM_stage)) +
  scale_color_manual(values = colors) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size=14),
    axis.title.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    legend.position = "none"
  )
ggsave("~/fig/tmp.pdf", g.tnmstage, height=4, width=12)


#### STRECK VS EDTA & storage time
dt[,streck:=factor(streck, levels=c(TRUE, FALSE), labels=c("Streck", "EDTA"))]
sample_sizes_streck <- dt[order(streck), .(n = .N), by = streck][, label := paste0(streck, "\n(n=", n, ")")]

dt[,storage_time:=factor(storage_time, levels=c("1 Day", "2 Day", "3+ Days"))]
sample_sizes_storage <- dt[order(storage_time), .(n = .N), by = storage_time][, label := paste0(storage_time, "\n(n=", n, ")")]

g.streck <- ggplot(dt, aes(x = streck, y = monitoring_score)) +
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
    title = "Blood collection tube",
    y = "DELFI-TF"
  ) +
  scale_x_discrete(labels = setNames(sample_sizes_streck$label, sample_sizes_streck$streck)) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size=14),
    axis.title.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    legend.position = "none"
  )
ggsave("~/fig/tmp.pdf", g.streck, height=3.5, width=4)

g.storage <- ggplot(dt[!is.na(storage_time)], aes(x = storage_time, y = monitoring_score)) +
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
    title = "Time to Storage",
    y = "DELFI-TF"
  ) +
  scale_x_discrete(labels = setNames(sample_sizes_storage$label, sample_sizes_storage$streck)) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size=14),
    axis.title.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    legend.position = "none"
  )
ggsave("~/fig/tmp.pdf", g.storage, height=3.5, width=5)

#### Smoking
dt[,smoking_status:=factor(smoking_status, levels=c("Never", "Former", "Current"))]
sample_sizes_smoking <- dt[order(smoking_status), .(n = .N), by = smoking_status][, label := paste0(smoking_status, "\n(n=", n, ")")]

g.smoking <- ggplot(dt, aes(x = smoking_status, y = monitoring_score)) +
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
    title = "Smoking Status",
    y = "DELFI-TF"
  ) +
  scale_x_discrete(labels = setNames(sample_sizes_smoking$label, sample_sizes_smoking$n)) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size=14),
    axis.title.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    legend.position = "none"
  )
ggsave("~/fig/tmp.pdf", g.smoking, height=3.5, width=6)

#### BMI
dt[,bmi_cat:= fcase(bmi < 18.5, "Underweight",
          bmi >= 18.5 & bmi < 25, "Normal weight",
          bmi >= 25 & bmi < 30, "Overweight",
          bmi >= 30, "Obese")]

dt[,bmi_cat:=factor(bmi_cat, levels=c("Underweight", "Normal weight", "Overweight", "Obese"))]
sample_sizes_bmi <- dt[order(bmi_cat), .(n = .N), by = bmi_cat][, label := paste0(bmi_cat, "\n(n=", n, ")")]

g.bmi <- ggplot(dt[!is.na(bmi_cat)], aes(x = bmi_cat, y = monitoring_score)) +
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
    title = "BMI Category",
    y = "DELFI-TF"
  ) +
  scale_x_discrete(labels = setNames(sample_sizes_bmi$label, sample_sizes_bmi$bmi_cat)) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size=14),
    axis.title.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    legend.position = "none"
  )
ggsave("~/fig/tmp.pdf", g.bmi, height=3.5, width=6)

#### Medications
dt[,medication:=fcase(no_medication == 1, "None",
                      anticoagulant == 1, "Anticoagulants",
                      antipyretic == 1, "Antipyretics",
                      nsaid == 1, "NSAIDs",
                      steroid==1, "Systemic steroids")]


dt[,medication:=factor(medication, levels=c("None", "Anticoagulants",
                                            "Antipyretics", "NSAIDs", "Systemic steroids"))]
sample_sizes_med <- dt[order(medication), .(n = .N), by = medication][, label := paste0(medication, "\n(n=", n, ")")]

#### medication-stage plot
g.med <- ggplot(dt[!is.na(medication)], aes(x = medication, y = monitoring_score)) +
  geom_boxplot(
    width = 0.6,
    fill = "gray90",
    color = "gray40",
    outlier.shape=NA,
    position = position_dodge(width = 0.75)
  ) +
  geom_jitter(
    aes(color=factor(no_medication)),
    size=1.5,
    width=0.1
  ) +
  scale_y_log10(
  breaks = c(0.001, 0.003, 0.01, 0.03,  0.1, 0.3),
  labels = c("0.1%", "0.3%", "1%", "3%", "10%", "30%")
  ) +
  labs(
    title = "Medication",
    y = "DELFI-TF"
  ) +
  scale_x_discrete(labels = setNames(sample_sizes_med$label, sample_sizes_med$medication)) +
  scale_color_manual(values = c("0" = "#00A0C6", "1" = "#E67E22")) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size=14),
    axis.title.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    legend.position = "none"
  )
ggsave("~/fig/tmp.pdf", g.med, height=4, width=7)

######
dt[,number_comorbidities:=factor(number_comorbidities, levels=c(0, 1, 2, 3))]
sample_sizes_comorb <- dt[order(number_comorbidities), .(n = .N), by = number_comorbidities][, label := paste0(number_comorbidities, "\n(n=", n, ")")]

#### medication-stage plot
g.comorb <- ggplot(dt[!is.na(number_comorbidities)], aes(x = number_comorbidities, y = monitoring_score)) +
  geom_boxplot(
    width = 0.6,
    fill = "gray90",
    color = "gray40",
    outlier.shape=NA,
    position = position_dodge(width = 0.75)
  ) +
  geom_jitter(
    aes(color=factor(number_comorbidities)),
    size=1.5,
    width=0.1
  ) +
  scale_y_log10(
  breaks = c(0.001, 0.003, 0.01, 0.03,  0.1, 0.3),
  labels = c("0.1%", "0.3%", "1%", "3%", "10%", "30%")
  ) +
  labs(
    title = "Number of comorbidities",
    y = "DELFI-TF"
  ) +
  scale_x_discrete(labels = setNames(sample_sizes_comorb$label, sample_sizes_comorb$number_comorbidities)) +
  scale_color_manual(values = c("0" = "#E67E22", "1" = "#00A0C6", "2" = "#00A0C6", "3" = "#00A0C6")) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size=14),
    axis.title.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    legend.position = "none"
  )
ggsave("~/fig/tmp.pdf", g.comorb, height=4, width=7)


#### Comorbidities
## COPD, CVD, autoimmune, diabetes
dt[,COPD:=factor(COPD, levels=c(0, 1), labels=c("No chronic lung\ndisease", "Chronic lung\ndisease"))]
sample_sizes_copd <- dt[order(COPD), .(n = .N), by = COPD][, label := paste0(COPD, "\n(n=", n, ")")]

dt[,CVD:=factor(CVD, levels=c(0, 1), labels=c("No CVD", "CVD"))]
sample_sizes_cvd <- dt[order(CVD), .(n = .N), by = CVD][, label := paste0(CVD, "\n(n=", n, ")")]

dt[,diabetes:=factor(diabetes, levels=c(0, 1), labels=c("No diabetes", "Diabetes"))]
sample_sizes_diab <- dt[order(diabetes), .(n = .N), by = diabetes][, label := paste0(diabetes, "\n(n=", n, ")")]

dt[,autoimmune:=factor(autoimmune, levels=c(0, 1), labels=c("No autoimmune\ndisease", "Autoimmune\ndisease"))]
sample_sizes_autoimmune <- dt[order(autoimmune), .(n = .N), by = autoimmune][, label := paste0(autoimmune, "\n(n=", n, ")")]

g.copd <- ggplot(dt[!is.na(COPD)], aes(x = COPD, y = monitoring_score)) +
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
    title = "Chronic lung disease",
    y = "DELFI-TF"
  ) +
  scale_x_discrete(labels = setNames(sample_sizes_copd$label, sample_sizes_copd$copd)) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size=14),
    axis.title.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    legend.position = "none"
  )
ggsave("~/fig/tmp.pdf", g.copd, height=3.5, width=3)

g.cvd <- ggplot(dt[!is.na(CVD)], aes(x = CVD, y = monitoring_score)) +
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
    title = "Cardiovascular disease",
    y = "DELFI-TF"
  ) +
  scale_x_discrete(labels = setNames(sample_sizes_cvd$label, sample_sizes_cvd$cvd)) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size=14),
    axis.title.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    legend.position = "none"
  )
ggsave("~/fig/tmp.pdf", g.cvd, height=3.5, width=3)

g.diab <- ggplot(dt[!is.na(diabetes)], aes(x = diabetes, y = monitoring_score)) +
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
    title = "Diabetes",
    y = "DELFI-TF"
  ) +
  scale_x_discrete(labels = setNames(sample_sizes_diab$label, sample_sizes_diab$diabetes)) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size=14),
    axis.title.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    legend.position = "none"
  )
ggsave("~/fig/tmp.pdf", g.diab, height=3.5, width=3)

g.autoimmune <- ggplot(dt[!is.na(autoimmune)], aes(x = autoimmune, y = monitoring_score)) +
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
    title = "Autoimmune",
    y = "DELFI-TF"
  ) +
  scale_x_discrete(labels = setNames(sample_sizes_autoimmune$label, sample_sizes_autoimmune$autoimmune)) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size=14),
    axis.title.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    legend.position = "none"
  )
ggsave("~/fig/tmp.pdf", g.autoimmune, height=3.5, width=3)

#### Sbtype
dt2[,subtype:=factor(subtype, levels=c("Adenocarcinoma", "Squamous cell carcinoma", "SCLC"), labels=c("LUAD", "LUSC", "SCLC"))]

subtype_colors <- c("LUAD" = "#E69F00", "LUSC" = "#56B4E9", "SCLC" = "#A83253")

dt2 <- dt2[!grepl("Unknown", stage)]
dt2[,stage_cat:=ifelse(grepl("^I$|^II$", stage), "Stage I & II", "Stage III & IV")]
sample_sizes_subtype <- dt2[order(subtype), .(n = .N), by = .(subtype, stage_cat)][, label := paste0(subtype, "\n(n=", n, ")"), by=stage_cat]

setkey(dt2, subtype, stage_cat)
setkey(sample_sizes_subtype, subtype, stage_cat)
dt2 <- dt2[sample_sizes_subtype, nomatch=NULL]

#### subtype
g.subtype <- ggplot(dt2, aes(x = label, y = monitoring_score)) +
  geom_boxplot(
    width = 0.6,
    fill = "gray90",
    color = "gray40",
    outlier.shape=NA,
    position = position_dodge(width = 0.75)
  ) +
  geom_jitter(
    aes(color=subtype),
    size=1.5,
    width=0.1
  ) +
  scale_y_log10(
  	breaks = c(0.001, 0.003, 0.01, 0.03,  0.1, 0.3),
  	labels = c("0.1%", "0.3%", "1%", "3%", "10%", "30%")
  ) +
  labs(
    title = NULL,
    y = "DELFI-TF"
  ) +
  #scale_x_discrete(labels = setNames(sample_sizes_subtype$label, sample_sizes_subtype$subtype)) +
  facet_grid(.~stage_cat, scale="free_x", space="free_x") +
  scale_color_manual(values = subtype_colors) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size=14),
    axis.title.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    legend.position = "none"
  )
ggsave("~/fig/tmp.pdf", g.subtype, height=3.5, width=7)

