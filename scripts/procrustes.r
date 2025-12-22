library(vegan)
library(tidyverse)
lema.pcs <- read_csv("lema_pcs.csv")
lema.pcs <- lema.pcs %>% filter(monitoring_score > 0.003)
clcgp.pcs <- read_csv("clcgp_pcs.csv")

# Assume:
# lema.pcs  = data frame of PC scores, with column "class"
# clcgp.pcs = data frame of PC scores, with column "class"
# And PCs are in columns like PC1, PC2, ..., PC8

# Step 1: How many samples per class are possible?
table(lema.pcs$subtype)
table(clcgp.pcs$subtype)

min_by_class <- pmin(
  table(lema.pcs$subtype),
  table(clcgp.pcs$subtype)
)

min_by_class
# LUAD = 170, LUSC = 70, SCLC = 40

# Step 2: Stratified random subsampling
set.seed(1234)

subsample <- function(df, nvec) {
  do.call(rbind, lapply(names(nvec), function(cl){
    df %>% filter(subtype == cl) %>% sample_n(nvec[cl])
  }))
}

lema_sub  <- subsample(lema.pcs,  min_by_class)
clcgp_sub <- subsample(clcgp.pcs, min_by_class)

# Remove the class column so only PC scores remain
PC1 <- lema_sub %>% select(c("pc1", "pc2", "pc3"))
PC2 <- clcgp_sub %>% select(c("pc1", "pc2", "pc3"))

# Step 3: Procrustes alignment
proc <- procrustes(PC1, PC2, scale = TRUE)

# Step 4: Permutation test (p-value)
set.seed(123)
proc_test <- protest(PC1, PC2, permutations = 9999)

proc_test
