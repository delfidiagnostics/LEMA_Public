# Base data
library(ggplot2)
library(ggalluvial)
library(tidyverse)
library(pROC)

lema_preds <- read_csv("lema_preds.csv") ## Cutoff -> -2.5 for SCLC, 0 for ADC
#lema_preds <- read_csv("lema_preds-0.2percentmaf.csv") ## Cutoff -> -2.5 for SCLC, 0 for ADC
#lema_preds <- read_csv("lema_preds-0.15percentmaf.csv") ## Cutoff -> -2.5 for SCLC, 0 for ADC


#### Threshold
# Create binary outcome (1 = LUSC, 0 = LUAD)
lema_preds2 <- lema_preds %>% filter(subtype!="SCLC")
true_label <- ifelse(lema_preds2$subtype == "LUSC", 1, 0)

# Find optimal threshold maximizing balanced accuracy
roc_obj <- roc(true_label, lema_preds2$predictions2)
coords_all <- coords(roc_obj, "all", ret = c("threshold", "sensitivity", "specificity"))
coords_all$balanced_acc <- (coords_all$sensitivity + coords_all$specificity) / 2

optimal_threshold <- coords_all$threshold[which.max(coords_all$balanced_acc)]
optimal_threshold <- coords(roc_obj, "best", best.method = "youden")$threshold


# Create binary outcome (1 = LUSC, 0 = LUAD)
true_label <- ifelse(lema_preds2$subtype == "LUSC", 1, 0)

# Get all thresholds with sensitivity and specificity
roc_obj <- roc(true_label, lema_preds2$predictions2)
coords_all <- coords(roc_obj, "all", ret = c("threshold", "sensitivity", "specificity"))

# Filter for thresholds meeting your constraint
viable_thresholds <- coords_all %>%
  filter(sensitivity >= 0.80, specificity >= 0.80)

if(nrow(viable_thresholds) > 0) {
  # Among viable thresholds, pick one that maximizes overall accuracy
  # or balanced accuracy
  viable_thresholds$balanced_acc <- (viable_thresholds$sensitivity + viable_thresholds$specificity) / 2
  optimal <- viable_thresholds[which.max(viable_thresholds$balanced_acc), ]

  cat("Optimal threshold:", optimal$threshold, "\n")
  cat("Sensitivity (LUSC):", optimal$sensitivity, "\n")
  cat("Specificity (LUAD):", optimal$specificity, "\n")
} else {
  cat("No threshold achieves ≥80% in both classes\n")
  # Find best compromise
  coords_all$min_acc <- pmin(coords_all$sensitivity, coords_all$specificity)
  best_compromise <- coords_all[which.max(coords_all$min_acc), ]
  cat("Best achievable minimum accuracy:", best_compromise$min_acc, "\n")
}
####

lema_preds <- lema_preds %>% mutate(first_pred=ifelse(predictions1 > -2.5, "SCLC", "NSCLC"),
#                      second_pred=ifelse(predictions2 > optimal_threshold , "LUSC", "LUAD"))
#                      second_pred=ifelse(predictions2 > -0.25 , "LUSC", "LUAD")) ## for 0.002 MAF
                      second_pred=ifelse(predictions2 > -0.15 , "LUSC", "LUAD")) ## for 0.003 MAF

lema_preds <- lema_preds %>% mutate(final_pred=case_when(first_pred=="NSCLC" & second_pred=="LUAD" ~ "LUAD",
                                           first_pred=="NSCLC" & second_pred=="LUSC" ~ "LUSC",
                                           first_pred=="SCLC" ~ "SCLC"))


table(lema_preds$subtype, lema_preds$final_pred)
total <- sum(table(lema_preds$subtype, lema_preds$final_pred))
sum(diag(table(lema_preds$subtype, lema_preds$final_pred)/total))

# Base data
data <- lema_preds %>% 
  select(original = subtype, predicted = final_pred) %>%
  count(original, predicted, name = "freq")

# Compute sample sizes
original_totals <- data %>%
  group_by(original) %>%
  summarise(n = sum(freq), .groups = "drop") %>%
  mutate(label = paste0(original, " (n=", n, ")"))

# Build label lookup
label_lookup <- setNames(original_totals$label, original_totals$original)

# Apply labels to original and predicted
data <- data %>%
  mutate(
    original_label = label_lookup[original],
    predicted_label = label_lookup[predicted] 
  ) %>%
  group_by(original_label) %>%
  mutate(percentage = paste0(round(freq / sum(freq), 2), "%")) %>%
  ungroup()

# Define stratum order with shared spacers
stratum_levels <- c(
  label_lookup["LUAD"],
  "SP1",
  label_lookup["LUSC"],
  "SP2",
  label_lookup["SCLC"]
)
predicted_levels <- c("LUAD", "SP1", "LUSC", "SP2", "SCLC")

# Add spacer rows with matching labels
spacers <- data.frame(
  original = paste0("SP", 1:2),
  predicted = paste0("SP", 1:2),
  freq = rep(20, 2),
  original_label = paste0("SP", 1:2),
  predicted_label = paste0("SP", 1:2),
  percentage = ""
)

# Combine full dataset
plot_data <- bind_rows(data, spacers)

# Set factor levels to align both sides perfectly
plot_data$original_label <- factor(plot_data$original_label, levels = stratum_levels)
plot_data$predicted_label <- factor(plot_data$predicted_label, levels = stratum_levels)

# Color palette (spacers transparent)
fill_colors <- c(
  setNames(c("#E69F00", "#56B4E9", "#A83253"), label_lookup),
  "SP1" = "transparent", "SP2" = "transparent"
)

color_colors <- c(
  setNames(c("black", "black", "black"), label_lookup),
  "SP1" = "transparent", "SP2" = "transparent"
)

# Plot
g <- ggplot(plot_data,
       aes(axis1 = original_label, axis2 = predicted_label, y = freq)) +
  geom_alluvium(aes(fill = original_label, color = original_label),
                width = 0.25, knot.pos = 0.5,
                curve_type = "quintic", alpha = 0.9) +
  geom_stratum(width = 0.25, fill = "white", color = NA) +
  geom_text(stat = "stratum",
            aes(label = ifelse(grepl("^SP", after_stat(stratum)), "", after_stat(stratum))),
            size = 4.2, fontface = "bold") +
  #geom_text(aes(label = percentage),
  #          stat = "alluvium", size = 3.5, color = "black", nudge_x = 0.15) +
  scale_x_discrete(limits = c("Histology original", "Prediction"),
                   expand = c(0.15, 0.15)) +
  scale_fill_manual(values = fill_colors) +
  scale_color_manual(values = color_colors) +
  theme_void(base_size = 13) +
  theme(
    legend.position = "none",
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    strip.text.x = element_text(size = 14, face = "bold")
  ) +
  ggtitle("Histology Classification Prediction Flow")


ggsave("~/fig/tmp.pdf", g, width=10)


###### ROC plots
library(data.table)
rocstats <- function(obs, score, levels) {
    roc <- pROC::roc
    roc <- suppressMessages(roc(obs, score,
                                levels=levels,
                                ci=TRUE))
    list(sens = rev(roc[["sensitivities"]]),
         spec = rev(roc[["specificities"]]),
         auc = roc$auc,
         lower = roc$ci[1],
         upper = roc$ci[3])
}

plotrocs <- function(data, mytheme=NULL, textsize=6.5) {
    data <- as.data.table(data)
    lab <- data[,.(lower=format(round(unique(lower),2), nsmall=2),
                   auc=format(round(unique(auc),2), nsmall=2),
                   upper=format(round(unique(upper), 2), nsmall=2))]
    lab <- lab[,
               `:=`(text=paste0("AUC: ", auc, " (",
                                lower, " - ", upper, ")"),
                 x=0.5, y=0.1)]#seq(0.4, 0.15, length.out=n))]

    A <- ggplot() +
#        geom_vline(xintercept=c(0.8),
#                   color="gray80", size=0.5, linetype="dashed") +
        geom_abline(intercept=1, slope=1, color="gray", size=1) +
        geom_line(data=data, aes(spec, sens), size=1.1,
                  color="#3b6ba5")

    A <- A + scale_x_reverse(expand=c(0, 0.01),
                             breaks=c(0, 0.25, 0.5, 0.80, 1),
                             labels=as.character(
                                 #c("1.0", ".75", ".50", ".20", "0"))) +
                                 c("0", ".25", ".50", ".80", "1.0"))) +
        scale_y_continuous(expand=c(0, 0.01),
                           labels=as.character(
                               c("0", ".25", ".50", ".75", "1.0"))) +
        ## scale_color_manual(values=palette) +
        mytheme +
        geom_text(data=lab, aes(x, y,  label=text),
                  size=textsize, hjust=0) +
        xlab("Specificity") +
        ylab("Sensitivity") +
        ##ggtitle(overall.title) +
        guides(color=guide_legend(title="Approach"))
    A
}

##### NSCLC VS SCLC ROC
nsclc_scores <- lema_preds %>% filter(category %in% c("NSCLC", "SCLC")) %>% select(category, predictions1)
rocs_nsclc <- rocstats(nsclc_scores$category, nsclc_scores$predictions1, levels=c("SCLC", "NSCLC"))

g.nsclc <- plotrocs(rocs_nsclc) +
    ggtitle("Predicting NSCLC vs. SCLC\nin cfDNA") +
    theme_classic(base_size = 13) + theme(#legend.position = "none",
                                          plot.title = element_text(size =  32),
                                          axis.title = element_text(size = 22),
                                          axis.text = element_text(size = 15),
                                          plot.margin = margin(r = 15, t = 5, b = 5, l = 5,unit = "pt")
    )
    ggsave("~/fig/tmp.pdf", g.nsclc)

##### LUAD VS LUSC ROC
luad_scores <- lema_preds %>% filter(subtype %in% c("LUAD", "LUSC")) %>% select(subtype, predictions2)
rocs_luad <- rocstats(luad_scores$subtype, luad_scores$predictions2, levels=c("LUAD", "LUSC"))

g.luad <- plotrocs(rocs_luad) +
    ggtitle("Predicting LUAD vs. LUSC\nin cfDNA") + theme_classic(base_size = 13) +
    theme_classic(base_size = 13) + theme(#legend.position = "none",
                                          plot.title = element_text(size =  32),
                                          axis.title = element_text(size = 22),
                                          axis.text = element_text(size = 15),
                                          plot.margin = margin(r = 15, t = 5, b = 5, l = 5,unit = "pt")
    )
ggsave("~/fig/tmp.pdf", g.luad)

##### Score distributions
mycolors <- c("#E69F00", "#56B4E9", "#A83253")
library(ggplot2)
g <- ggplot(lema_preds, aes(monitoring_score, predictions1+2.5, color=subtype)) 
#g <- g + geom_rect(aes(xmin=0.003, xmax=0.01, ymin =-Inf, ymax=Inf), fill="gray90" , color="gray90")
g <- g + geom_point(size=2)
g <- g + geom_hline(yintercept=0, linetype="dashed", color="blue", linewidth=1.1)
g <- g + theme_classic(base_size = 13) + theme(#legend.position = "none",
                                               plot.title = element_text(size =  32),
                                               axis.title = element_text(size = 22),
                                               axis.text = element_text(size = 15),
                                               plot.margin = margin(r = 15, t = 5, b = 5, l = 5,unit = "pt")
                                               ) +
    ylab("SCLC score") + xlab("DELFI-TF") + ggtitle("NSCLC vs. SCLC") +
    scale_x_log10(limits = c(0.003, 0.6), expand = c(0, 0), 
                 breaks = c(0.003, 0.01, 0.03, 0.1, 0.2, 0.35, 0.6),
                 labels = c("0.3%", "1%", "3%", "10%", "20%", "35%","60%")) +
    scale_color_manual(values = c(
        "LUAD" = "#E69F00",
        "LUSC" = "#56B4E9",
        "SCLC" = "#A83253"
    )) + ylim(c(-5.2, 6))

ggsave("~/fig/tmp.pdf", g, width=6)

h <- ggplot(lema_preds %>% filter(category=="NSCLC"), aes(monitoring_score, predictions2 + 0.15, color=subtype)) 
#h <- ggplot(lema_preds %>% filter(category=="NSCLC"), aes(monitoring_score, predictions2 + optimal_threshold, color=subtype)) 
#h <- h + geom_rect(aes(xmin=0.002, xmax=0.003, ymin =-Inf, ymax=Inf), fill="gray90" , color="gray90")
h <- h + geom_point(size=2)
h <- h + geom_hline(yintercept=0, linetype="dashed", color="blue", linewidth=1.1)
h <- h + ylab("LUSC score") + xlab("DELFI-TF") + ggtitle("LUAD vs. LUSC")
h <- h + theme_classic(base_size = 13) + theme(legend.position = "none",
                                               plot.title = element_text(size =  32),
                                               axis.title = element_text(size = 22),
                                               axis.text = element_text(size = 15),
                                               plot.margin = margin(r = 15, t = 5, b = 5, l = 5,unit = "pt")
                                               ) +
    scale_x_log10(limits = c(0.003, 0.6), expand = c(0, 0), 
                 breaks = c(0.003, 0.01, 0.03, 0.1, 0.2, 0.35, 0.6),
                 labels = c("0.3%", "1%", "3%", "10%", "20%", "35%","60%")) +
    scale_color_manual(values = c(
        "LUAD" = "#E69F00",
        "LUSC" = "#56B4E9",
        "SCLC" = "#A83253"
    )) + ylim(c(-5.2, 6))
ggsave("~/fig/tmp.pdf", h, width=5)
