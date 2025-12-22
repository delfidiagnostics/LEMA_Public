library(tidyverse)
library(readxl)
library(tidymodels)
library(recipes)
library(data.table)
library(pROC)
library(ggpubr)


#### Use all samples (including below .3% TF) but build LUAD/LUSC model that
#### prefers fragmentation above LOB and proteins below?

proteins <- readRDS("../data/reference/proteins_zscores_lema.rds")
proteins <- proteins |> select(id=`DELFI-ID`, cea, cyfra, he4, ca125, ca153)
meta <- read_excel("../data/reference/lema_master_sheet.xlsx")
### LEMA metadata
setDT(meta)
ids <- meta[,c("oc", "seqID4x", "CNCLBL", "LCNSTAGE", "NTTBAS", "LCNCAT")]

setnames(ids, c("patient", "id", "type", "stage", "subtype", "category"))
#ids <- ids[type=="Lung Cancer" & grepl("^II$|^III$|^IV$", stage) &
#           grepl("Adenocarcinoma|Squamous", subtype)]
ids <- ids[type=="Lung Cancer" & (grepl("Adenocarcinoma|Squamous", subtype) |
                                  grepl("^SCLC$", category))]

ids[category=="SCLC", subtype:="SCLC"]

ids2 <- ids[,.(patient, id)]
proteins2 <- as.data.table(proteins)
setkey(ids2, id)
setkey(proteins2, id)
ids2[proteins2]

### Monitoring score
dtf <- fread("../data/reference/monitoring_score.4x.v3.1.tsv")
setnames(dtf, "entity_id", "id")

#dtf <- fread("~/projects/lemaAnalysis/monitoring_score.10x.v3.1.tsv")
#setnames(dtf, c("id", "monitoring_score"))
#dtf[,id:=gsub("AMSeq2", "AMSeq1", id)]

setkey(dtf, id)
setkey(ids, id)
ids <- ids[dtf, nomatch=NULL]


####
lema_pcs <- read_csv("../data/derived_data/lema_pcs.csv") ## Cutoff -> -2.5 for SCLC, 0 for ADC
#lema_preds <- lema_preds %>% mutate(first_pred=ifelse(predictions1 > -2.5, "SCLC", "NSCLC"),
#                      second_pred=ifelse(predictions2 > 0, "LUSC", "LUAD"))
#
#lema_preds <- lema_preds %>% mutate(final_pred=case_when(first_pred=="NSCLC" & second_pred=="LUAD" ~ "LUAD",
#                                           first_pred=="NSCLC" & second_pred=="LUSC" ~ "LUSC",
#                                           first_pred=="SCLC" ~ "SCLC"))
#
features <- lema_pcs %>% select(patient, id, starts_with("pc"), stage, subtype, category, dtf=monitoring_score)

features.nsclc <- features %>% mutate(category = factor(category, levels=c("NSCLC", "SCLC")))

features.luad <- features %>% filter(category != "SCLC") %>% select(-category) %>% 
    mutate(subtype = factor(subtype, levels=c("LUAD", "LUSC")))
features.luad <- inner_join(features.luad, proteins, by="id")

protein.table <- features.luad %>% select(patient, cea, cyfra, he4, ca125, ca153)
write_csv(protein.table, "../data/derived_data/lema_proteins_supp.csv")

summary(glm(subtype ~ log10(cyfra) + log10(he4) + log10(cea) + log10(ca125), data=features.luad , family="binomial"))
summary(glm(subtype ~ pc1 + log10(cyfra) + log10(he4) + log10(cea) + log10(ca125), data=features.luad , family="binomial"))

# Define feature groups
aneuploidy_features <- c("pc1", "pc2", "pc3", "pc4", "pc5", "pc6", "pc7", "pc8")
protein_features <- c("he4", "cyfra", "cea", "ca125", "ca153")

## model
#logit(P) = β₀ + β₁(TF > 0.3%) + β₂(aneuploidy features) + β₃(protein features) +
#           β₄(TF > 0.3%) × (aneuploidy features)
# Create recipe with interaction terms matching your model specification
library(hardhat)
table(features.luad$subtype)
features.luad2 <- features.luad %>% select(-patient)
features.luad2 <- features.luad2 %>%
  mutate(case_wts = ifelse(subtype == "LUSC", 261/100, 1)) %>%
  mutate(case_wts = importance_weights(case_wts))

features.luad2 <- features.luad2 %>% mutate(dtf2 = ifelse(dtf <= 0.003, 0, dtf))
#### CN signature recipe
recipe_copynumber <- recipe(subtype ~ ., data = features.luad2) %>%
    update_role(id, new_role = "ID") %>%
    update_role(stage, new_role = "STAGE") %>%
    update_role(dtf, new_role = "TF") %>%
    step_rm(he4, cyfra, cea, ca125, ca153, dtf)

#### Protein recipe
recipe_protein <- recipe(subtype ~ ., data = features.luad2) %>%
    update_role(id, new_role = "ID") %>%
    update_role(stage, new_role = "STAGE") %>%
    update_role(dtf, new_role = "TF") %>%
    step_log(he4, cyfra, cea, ca125, ca153, base=10) %>%
    step_rm(pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8, dtf, dtf2)


recipe_multimodal <- recipe(subtype ~ ., data = features.luad2) %>%
    update_role(id, new_role = "ID") %>%
    update_role(stage, new_role = "STAGE") %>%
    step_log(he4, cyfra, cea, ca125, ca153, base = 10) %>%
    step_interact(~ dtf2:(pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8)) %>%
    step_rm(dtf)

# Define regularized logistic regression model with glmnet
logistic_spec <- logistic_reg(penalty = tune(), mixture = tune()) %>%  # mixture=1 is lasso, 0 is ridge
  set_engine("glmnet") %>%
  set_mode("classification")

# Create workflow (multimodal)
workflow_multimodal <- workflow() %>%
    add_recipe(recipe_multimodal) %>%
    add_model(logistic_spec) %>% 
    add_case_weights(case_wts)

# Create workflow (copy number)
workflow_copynumber <- workflow() %>%
    add_recipe(recipe_copynumber) %>%
    add_model(logistic_spec) %>% 
    add_case_weights(case_wts)

# Create workflow (protein)
workflow_protein<- workflow() %>%
    add_recipe(recipe_protein) %>%
    add_model(logistic_spec) %>% 
    add_case_weights(case_wts)

grid_vals <- grid_regular(
  penalty(range = c(-3, -1)),     # 10^-5 to 10^-2 (small λ)
  mixture(range = c(0, 1)),
  levels = 5
)
# Set up cross-validation
set.seed(123)
cv_folds <- vfold_cv(features.luad2, v = 10, strata=subtype)

# Perform cross-validation with final model
cv_results_multimodal <- workflow_multimodal %>%
      tune_grid(
    resamples = cv_folds,
    grid = grid_vals,  # your defined grid of penalty + mixture
    metrics = metric_set(bal_accuracy, roc_auc, sensitivity, specificity),
    control = control_grid(save_pred = TRUE)
  )

cv_results_copynumber <- workflow_copynumber %>%
      tune_grid(
    resamples = cv_folds,
    grid = grid_vals,  # your defined grid of penalty + mixture
    metrics = metric_set(bal_accuracy, roc_auc, sensitivity, specificity),
    control = control_grid(save_pred = TRUE)
    )

cv_results_protein <- workflow_protein %>%
      tune_grid(
    resamples = cv_folds,
    grid = grid_vals,  # your defined grid of penalty + mixture
    metrics = metric_set(bal_accuracy, roc_auc, sensitivity, specificity),
    control = control_grid(save_pred = TRUE)
  )


# View cross-validation results
cv_metrics_copynumber <- collect_metrics(cv_results_copynumber)
best_params_copynumber <- select_best(cv_results_copynumber, metric="roc_auc")

cv_metrics_protein <- collect_metrics(cv_results_protein)
best_params_protein <- select_best(cv_results_protein, metric="roc_auc")

cv_metrics_multimodal <- collect_metrics(cv_results_multimodal)
best_params_multimodal <- select_best(cv_results_multimodal, metric="roc_auc")

# Collect predictions from each CV fold for each model
cv_preds_copynumber <- collect_predictions(cv_results_copynumber) %>% filter(.config==best_params_copynumber$.config)
cv_preds_protein <- collect_predictions(cv_results_protein) %>% filter(.config==best_params_protein$.config)
cv_preds_multimodal <- collect_predictions(cv_results_multimodal) %>% filter(.config==best_params_multimodal$.config)

# Calculate ROC curves with CI using pROC
roc_copynumber <- roc(cv_preds_copynumber$subtype,
                      cv_preds_copynumber$.pred_LUSC,
                      ci = TRUE, levels = c("LUSC", "LUAD"))

roc_protein <- roc(cv_preds_protein$subtype,
                   cv_preds_protein$.pred_LUSC,
                   ci = TRUE, levels = c("LUSC", "LUAD"))

roc_multimodal <- roc(cv_preds_multimodal$subtype,
                      cv_preds_multimodal$.pred_LUSC,
                      ci = TRUE, levels = c("LUSC", "LUAD"))

roc.test(roc_copynumber, roc_multimodal)
roc.test(roc_protein, roc_multimodal)

# Get ROC curves using tidymodels

roc_copynumber <- cv_preds_copynumber %>%
  roc_curve(truth = subtype, .pred_LUSC)

roc_protein <- cv_preds_protein %>%
  roc_curve(truth = subtype, .pred_LUSC)

roc_multimodal <- cv_preds_multimodal %>%
  roc_curve(truth = subtype, .pred_LUSC)

# Combine for plotting
roc_data <- bind_rows(
  roc_copynumber %>% mutate(model = "Copy Number"),
  roc_protein %>% mutate(model = "Protein"),
  roc_multimodal %>% mutate(model = "Multimodal")
)

# Extract AUC and CI
#auc_copynumber <- sprintf("%.3f (%.3f-%.3f)",
#                          roc_copynumber$auc,
#                          roc_copynumber$ci[1],
#                          roc_copynumber$ci[3])

auc_copynumber <- sprintf("%.2f (%.2f-%.2f)",
                          cv_metrics_copynumber$mean[2],
                          cv_metrics_copynumber$mean[2] - 1.96 * cv_metrics_copynumber$std_err[2],
                          cv_metrics_copynumber$mean[2] + 1.96 * cv_metrics_copynumber$std_err[2]
                          )

auc_protein <- sprintf("%.2f (%.2f-%.2f)",
                          cv_metrics_protein$mean[2],
                          cv_metrics_protein$mean[2] - 1.96 * cv_metrics_protein$std_err[2],
                          cv_metrics_protein$mean[2] + 1.96 * cv_metrics_protein$std_err[2]
                          )

auc_multimodal <- sprintf("%.2f (%.2f-%.2f)",
                          cv_metrics_multimodal$mean[2],
                          cv_metrics_multimodal$mean[2] - 1.96 * cv_metrics_multimodal$std_err[2],
                          cv_metrics_multimodal$mean[2] + 1.96 * cv_metrics_multimodal$std_err[2]
                          )


lab <- data.table(model=c("Multimodal", "cfDNA", "Protein"),
				  AUC=c(auc_multimodal, auc_copynumber, auc_protein),
				  x = 0.6, y=seq(0.2, 0.1, length.out=3))
lab <- lab[, text := paste0(model, " AUC: ", AUC)]


roc_copynumber <- cv_preds_copynumber %>% mutate(subtype = factor(subtype, levels=c("LUSC", "LUAD"))) %>%
  roc_curve(truth = subtype, .pred_LUAD)

roc_protein <- cv_preds_protein %>% mutate(subtype = factor(subtype, levels=c("LUSC", "LUAD"))) %>%
  roc_curve(truth = subtype, .pred_LUAD)

roc_multimodal <- cv_preds_multimodal %>% mutate(subtype = factor(subtype, levels=c("LUSC", "LUAD"))) %>%
  roc_curve(truth = subtype, .pred_LUAD)

# Combine for plotting
roc_data <- bind_rows(
  roc_copynumber %>% mutate(model = "cfDNA"),
  roc_protein %>% mutate(model = "Protein"),
  roc_multimodal %>% mutate(model = "Multimodal")
)

g <- ggplot(roc_data, aes(x = 1 - specificity, y = 1 - sensitivity, color = model)) +
  #geom_vline(xintercept=c(0.8),
  #           color="gray80", size=0.5, linetype="dashed") +
  geom_line(size=1.4) +
  scale_x_reverse(expand=c(0, 0.01),
				  breaks=c(0, 0.25, 0.5, 0.80, 1),
				  labels=as.character(
									  #c("1.0", ".75", ".50", ".20", "0"))) +
						  c("0", ".25", ".50", ".80", "1.0"))) +
						  scale_y_continuous(expand=c(0, 0.01),
							labels=as.character( c("0", ".20", ".50", ".75", "1.0"))) +

  #      scale_color_manual(values=c(rgb(171, 164, 222, maxColorValue=255),
  #                         rgb(117, 187, 149, maxColorValue=255))) +
  theme_classic(base_size=13) +
  theme(panel.grid=element_blank(),
        #           legend.position=c(0.6, 0.2),
        axis.title = element_text(size = 25),
        axis.text = element_text(size = 18),
        legend.position="none",
        aspect.ratio=0.8,
        legend.text.align=1,
        #           legend.title=element_text(size=10),
        legend.text=element_blank()) +
  xlab("Specificity") + ylab("Sensitivity") +
  geom_text(data=lab, aes(x, y, label=text), size=5.5, hjust=0, fontface="bold") +
  #geom_label(data = lab, aes(x = x, y = y, label = text, fill = model),
  #           hjust = 0, color = "white", fontface = "bold", size = 5,
  #           label.padding = unit(0.3, "lines")) +
  geom_abline(intercept = 1, slope = 1, color="gray") +
  scale_color_manual(values = c("cfDNA" = "#A6761D",
                                  "Protein" = "#666666",
                                  "Multimodal" = "#1B9E77")) +
  guides(color=guide_legend(title="Analyte"),
         fill = guide_legend(override.aes = list(colour = NULL)))
ggsave("~/fig/tmp.pdf", g)

# Print AUCs for verification
cat("AUC Values with 95% CI:\n")
cat("Copy Number:", auc_copynumber, "\n")
cat("Protein:", auc_protein, "\n")
cat("Multimodal:", auc_multimodal, "\n")

# Get predictions averaged across repeats
cv_preds_multimodal <- collect_predictions(cv_results_multimodal) %>%
  group_by(.row) %>%
  summarize(
    subtype = first(subtype),
    .pred_LUAD = mean(.pred_LUAD),
    .pred_LUSC = mean(.pred_LUSC)
  ) %>%
  ungroup()

# Find threshold that gives 80% sensitivity for LUAD
luad_samples <- cv_preds_multimodal %>% filter(subtype == "LUAD")

thresholds <- seq(0, 1, by = 0.001)
sensitivities <- sapply(thresholds, function(thresh) {
  preds <- if_else(luad_samples$.pred_LUAD > thresh, "LUAD", "LUSC")
  mean(preds == "LUAD")
})

# Find threshold closest to 80% sensitivity
target_threshold <- thresholds[which.min(abs(sensitivities - 0.80))]
achieved_sensitivity <- sensitivities[which.min(abs(sensitivities - 0.80))]

cat("Threshold for 80% LUAD sensitivity:", target_threshold, "\n")
cat("Achieved LUAD sensitivity:", achieved_sensitivity, "\n")

# Apply threshold
cv_preds_multimodal$.pred_class <- if_else(cv_preds_multimodal$.pred_LUAD > target_threshold,
                                           "LUAD", "LUSC")
cv_preds_multimodal$.pred_class <- factor(cv_preds_multimodal$.pred_class,
                                          levels = levels(cv_preds_multimodal$subtype))



#### 
# Step 1: Reshape and format
features_long <- features.luad2 %>%
  select(id, subtype, cea, cyfra, he4, ca125, ca153) %>%
  pivot_longer(
    cols = -c(id, subtype),
    names_to = "protein",
    values_to = "expression"
  ) %>%
  mutate(
    protein = case_when(
      protein == "cea"    ~ "CEA",
      protein == "cyfra"  ~ "CYFRA21-1",
      protein == "he4"    ~ "HE4",
      protein == "ca125"  ~ "CA-125",
      protein == "ca153"  ~ "CA-153"
    ),
    protein = factor(protein, levels = c("CEA", "CYFRA21-1", "HE4", "CA-125", "CA-153")),
    subtype = factor(subtype, levels = c("LUAD", "LUSC"))
  )

# Step 2: Custom colors
subtype_colors <- c("LUAD" = "#E69F00", "LUSC" = "#56B4E9")

label_positions <- features_long %>%
  group_by(protein) %>%
  summarise(label.y = max(expression, na.rm = TRUE) * 1.5)  # push above max point

# Step 2: Create comparisons manually (LUAD vs LUSC)
comparisons <- list(c("LUAD", "LUSC"))

my_p_label <- function(p) {
  case_when(
    p < 0.001 ~ "p < 0.001",
    p < 0.01 ~ "p < 0.01",
    p < 0.05 ~ "p < 0.05",
    TRUE     ~ "ns"
  )
}

pvals <- features_long %>%
  group_by(protein) %>%
  summarise(
    p = wilcox.test(expression ~ subtype)$p.value,
    .groups = "drop"
  ) %>%
  mutate(p_label = my_p_label(p))

pval_labels <- left_join(pvals, label_positions, by = "protein")

# Step 3: Plot with manual label.y
g.proteins <- ggplot(features_long, aes(x = subtype, y = expression)) +
  geom_boxplot(
    width = 0.6,
    fill = "gray90",
    color = "gray40",
    outlier.shape=NA,
    position = position_dodge(width = 0.75)
  ) +
  geom_jitter(
    aes(color = subtype),
    size = 1.2,
    position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75)
  ) +
  facet_wrap(~ protein, nrow = 1, scales = "free_y") +
  scale_y_log10() +
  scale_color_manual(values = subtype_colors) +
  geom_text(
            data = pval_labels,
            aes(x = 1, y = label.y, label = p_label),
            inherit.aes = FALSE,
            color = "red",
            size = 3.5,
            hjust=0
            ) +
  labs(
    title = NULL,
    y = "Concentration",
    color = "Subtype"
  ) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size=14.5, hjust=0),
    axis.title.x = element_blank(),
    legend.position = "none",
    axis.title.y = element_text(size=19),
    axis.text.y = element_text(size=14),
    axis.text.x = element_text(angle = 45, hjust = 1, size=14)
  )
ggsave("~/fig/tmp.pdf", g.proteins, height=5, width=10)
