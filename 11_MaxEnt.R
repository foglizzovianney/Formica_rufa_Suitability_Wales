# Script: 11_MaxEnt_.R
# Description: MaxEnt model 


setwd("C:/Users/2210701/OneDrive - Swansea University/Desktop/Formica_rufa_wales")

library(terra)
library(tidyverse)
library(ENMeval)
library(dismo)
library(rJava)

# LOAD DATA 
predictors <- rast("spatial/vom/predictor_stack_wales.tif")
bg_bias <- read.csv("data/processed/bg_bias_n10000.csv")
occs <- read.csv("data/processed/occs_maxent.csv")

# EXTRACT PREDICTOR VALUES 
occs_pred <- terra::extract(predictors, occs[, c("x", "y")], ID = FALSE)
occs_full <- cbind(occs, occs_pred) %>% drop_na()
cat("Occurrence points with complete data:", nrow(occs_full), "\n")

bg_pred <- terra::extract(predictors, bg_bias[, c("x", "y")], ID = FALSE)
bg_full <- cbind(bg_bias, bg_pred) %>% drop_na()
cat("Background points with complete data:", nrow(bg_full), "\n")

# SAVE EXTRACTED DATA 
# Save so we don't need to redo if it crashes
write.csv(occs_full, "data/processed/occs_full_extracted.csv", row.names = FALSE)
write.csv(bg_full, "data/processed/bg_full_extracted.csv", row.names = FALSE)


# RUN FULL ENMeval
# Following Robinson et al. 2025 exactly
# 6 feature classes x 19 regularization values = 114 combinations
# 10-fold cross validation

e.mx.full <- ENMevaluate(
  occs = occs_full[, c("x", "y", names(predictors))],
  bg = bg_full[, c("x", "y", names(predictors))],
  algorithm = 'maxent.jar',
  partitions = 'randomkfold',
  partition.settings = list(kfolds = 10),
  tune.args = list(fc = c("L", "Q", "LQ", "H", "QH", "LQH"),
                   rm = seq(1, 10, 0.5)),
  parallel = FALSE)

saveRDS(e.mx.full, "outputs/ENMeval_results_full.RDS")


# MODEL SELECTION 
res <- eval.results(e.mx.full)
write.csv(res, "outputs/tables/ENMeval_results_full.csv", row.names = FALSE)

opt.step <- res %>%
  filter(!is.na(AICc)) %>%
  filter(AICc < 2000) %>%
  filter(or.10p.avg == min(or.10p.avg, na.rm = TRUE)) %>%
  filter(auc.val.avg == max(auc.val.avg, na.rm = TRUE)) %>%
  filter(AICc == min(AICc, na.rm = TRUE))

if(nrow(opt.step) > 1) opt.step <- opt.step[1, ]

best_model_full <- eval.models(e.mx.full)[[opt.step$tune.args]]

# SAVE BEST MODEL 
saveRDS(best_model_full, "outputs/best_model_full.RDS")

# VARIABLE IMPORTANCE 
var_importance <- data.frame(
  best_model_full@results[grep("permutation.importance",
                               rownames(best_model_full@results)),])
names(var_importance) <- "perm_imp"
var_importance$variable <- gsub(".permutation.importance", 
                                rownames(var_importance), 
                                replacement = "")
var_importance <- var_importance %>% arrange(desc(perm_imp))

print(var_importance)
write.csv(var_importance,
          "outputs/tables/variable_importance_full.csv",
          row.names = FALSE)


