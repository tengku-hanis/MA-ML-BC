# Meta-analysis of ML model on mammography

# Set wd
setwd("C:/Tengku/PhD/publication/thesis-SR_MA/analysis")

# Load packages ----
library(readxl)
library(meta)
library(mada)
library(magrittr)
library(dmetatools)
library(tidyverse)

# Read data
ma_data <- read_excel("C:/Tengku/PhD/publication/thesis-SR_MA/MA data.xlsx", 
                      sheet = "meta-analysis")

# Remove data with no conf matrix
glimpse(ma_data)

ma_data2 <- 
  ma_data %>% 
  mutate(across(c(year, n, starts_with("size"), n_test_validation, accuracy:FN), 
                as.numeric),
         across(c(author_year, year, source, origin_data, type, classifier), 
                as_factor)) %>% 
  filter(!TP < 0) 

# Explore ----
## Factorise classifier
ma_data2 %>% 
  mutate(across(c(author_year, year, source, origin_data, type, classifier), fct_drop)) %>%
  select(author_year, year, source, origin_data, type, classifier) %>% 
  map(fct_count) # 36 studies

ma_data2$classifier %>% 
  fct_drop() %>% 
  fct_count() %>% 
  ggplot(aes(f,n)) +
  geom_col()

ma_data2 %<>% 
  mutate(classifier = fct_collapse(classifier, 
                                   NN = c("ANN", "BNN", "FNN", "NN"), 
                                   DL = c("CNN", "MLP", "SCNN", "TCNN", "ResNet-50", "IceptionResNet-V2", "MLP",
                                          "Inception-v3", "DCNN"),
                                   Tree_based = c("CT", "RF", "random_forest"),
                                   SVM = c("SVM", "FSVM", "CVM"),
                                   Bayes_based = c("NB", "Bayes-class")
  ),
  classifier2 = fct_collapse(classifier,
                             NonDL = c("NN", "Tree_based", "SVM", "Bayes_based", "GMM", "KNN", 
                                       "LDA", "Logistic")
  ),
  type2 = fct_collapse(type, 
                       image = c("mammo", "mammo-CEDM", "tomo"), 
                       tabular = c("mammo_data", "radiomic-mammo", "mammo-features")
  ),
  source2 = fct_collapse(source,
                         MIAS = c("MIAS", "mini-MIAS"),
  ),
  source3 = fct_collapse(source2,
                         database = c("DDSM", "IN-BREAST", "MIAS", "na", "MMD", "DDSM+MIAS")
  ),
  country = fct_collapse(origin_data,
                         others = c("iran", "portugal", "jordan", "china", "na", "korea", 
                                    "serbia"), 
  ),
  across(c(author_year, year, source, origin_data, type, classifier, classifier2, type2, 
           source2, source3, country), 
         fct_drop)
  ) 
glimpse(ma_data2)

ma_data2 %>% 
  select(classifier, classifier2, type2, source2, source3, country) %>% 
  map(fct_count) 

## Check for relatively similar models in a study 
# make unique id
ma_data2$sid <- 1:77

# group by author_year and classifier
same_model <- 
  ma_data2 %>% 
  group_by(author_year, classifier) %>% 
  summarise(n = n()) %>% 
  filter(n > 1)
same_model
sum(same_model$n) #41 model

# filter to check one-by-one
same_model2 <- 
  ma_data2 %>% 
  filter(author_year %in% same_model$author_year & classifier %in% same_model$classifier) %>% 
  mutate(accuracy2 = (TP+TN)/(TP+TN+FP+FN) *100) %>% 
  select(sid, author_year, classifier, accuracy, accuracy2, add_note, origin_data, source, predict_class) 

# remove relatively similar model based on accuracy
ma_data3 <- 
  ma_data2 %>% 
  filter(!(sid %in% c(4,5,7,8, #al-antari2020
                      11, #al-afifi2020
                      25,26, #danala2018
                      35,36))) #jebamony2020

# Study ID for the analysis
ma_data3$sid2 <- rep(paste0("model", 1:length(ma_data3$author_year)))

# Descriptive ----
descrip <- madad(ma_data3 %>% 
                   rename(names = sid2), digits = 2)
descrip

old.par <- par()
plot.new()

par(fig = c(0, 0.5, 0, 1), new = TRUE)
forest(descrip, type = "sens", sname = ma_data3$sid2, cex = 0.8, main = "Sensitivity (95% CI)")

par(fig = c(0.5, 1, 0, 1),  new = TRUE)
forest(descrip, type = "spec", sname = ma_data3$sid2, cex = 0.8, main = "Specificity (95% CI")

par(old.par)

# OVERALL ANALYSIS ----
# Univariate ----
## DOR ----
dor_mada <- madauni(ma_data3, method = "DSL") # MH=fixed effect, DSL=random effect (default)
summary(dor_mada)

dor_mada$descr$names <- ma_data3$sid2
forest(dor_mada, polycol = "lightgrey", cex = 0.8, 
       main = "Diagnostic odds ratio model (95% CI)") # in log scale

# Bivariate ----
## Reitsma model (to get AUC) ----
fit_reitsma <- reitsma(ma_data3)
summary(fit_reitsma)

## Bootstrap CI
set.seed(2021)
with(ma_data3, AUC_boot(TP, FP, FN, TN))
#AUC=0.8984889
#CI=0.8452703, 0.9034844 

plot(fit_reitsma, sroclwd = 2, sroclty = 3, main = "")
points(fpr(ma_data3), sens(ma_data3), pch = 2)
legend("bottomright", c("data", "summary estimate", "SROC", "conf. region", 
                        paste0("AUC=0.90 (95% CI: 0.85, 0.90)")), 
       pch = c(2,1,1000,1000,1000,1000), lty = c(-1,-1,3,1,-1,-1))
# substantial heterogeneity between studies

# Heterogeneity ----
# 4 methods to check heterogeneity (heterogeneity is suspected if):
# - Visual:
# 1. Asymmetrical SROC
# 2. Individual studies in SROC are largely scatter
# 3. between-study variation is greater than within-study variation in the forest plot
#    of sensitivity, specificity, DOR
# - Statistical:
# 4. Correlation coefficient of sensitivity and specificity is larger than zero

## Calculate sensitivity and specificity
cor_data <- data.frame(id=1:68)
cor_data$sn <- ma_data3$TP/(ma_data3$TP+ma_data3$FN)
cor_data$sp <- ma_data3$TN/(ma_data3$FP+ma_data3$TN)

## Transforms to logit function
cor_data$logitsn <- log(cor_data$sn/(1-cor_data$sn))
cor_data$logitsp <- log(cor_data$sp/(1-cor_data$sp))

## Correlation coefficient
cor_data2 <- 
  cor_data %>% 
  mutate(across(c(sn, sp), round, digits = 3)) %>% 
  filter(sn != 1 & sp != 1)

cor(cor_data2$logitsn, cor_data2$logitsp)

# Publication bias ----
## Deek's test - metafor package
library(metafor)
dat <- escalc(measure="OR", ai=TP, bi=FP, ci=FN, di=TN, data=ma_data3, add=1/2, to="all")
dat$nn <- dat$TN + dat$FP
dat$np <- dat$TP + dat$FN


dat$vi <- 1/(4*dat$np) + 1/(4*dat$nn)
tmp <- rma(yi, vi, data=dat, method="DL")
tmp

sav <- funnel(tmp, atransf=exp, xlim=c(0,7), at=log(c(1, 10, 100, 1000, 10000)), 
              ylim=c(0,0.3), steps=4, back=NA, level=NA, lwd=2, lty=1, refline=coef(tmp), 
              ylab="1/root(ess)", xlab = "Diagnostic odds ratio", label = TRUE, cex=0.8)
legend("bottomright", "Deeks test:\nP-value = 0.002", bty = "n")

reg <- regtest(tmp, model="lm")
reg

ys <- seq(0, 0.3, length=100)
lines(coef(reg$fit)[1] + coef(reg$fit)[2]*ys, ys, lwd=2, lty=3)

## Deek's test - meta package
dor_meta <- metabin(TP,TP+FP,FN,FN+TN, sm="DOR", comb.fixed=FALSE,comb.random=TRUE, 
                    method="Inverse", author_year,
                    data=ma_data3, subset = -63)
print(dor_meta, digits = 2)

metabias(dor_meta, method.bias = "Deeks", plotit = TRUE)

# SUBGROUP ANALYSIS ----
# Meta-regression ----
## Meta-regression ----
null_mod <- reitsma(ma_data3, formula = cbind(tsens, tfpr) ~ 1, method = "ml")

reg_class <- reitsma(ma_data3, formula = cbind(tsens, tfpr) ~ classifier, method = "ml")
summary(reg_class)

reg_type <- reitsma(ma_data3, formula = cbind(tsens, tfpr) ~ type2, method = "ml")
summary(reg_type)

reg_source <- reitsma(ma_data3, formula = cbind(tsens, tfpr) ~ source2, method = "ml")
summary(reg_source)

reg_country <- reitsma(ma_data3, formula = cbind(tsens, tfpr) ~ country, method = "ml")
summary(reg_country)

### compare 2 models
anova(null_mod, reg_class) #not sig classifier2
anova(null_mod, reg_country)
anova(null_mod, reg_source) #not sig source3

anova(null_mod, reg_type) #not sig

# classifier, country, and source explain some heterogeneity between primary studies

# Compare SROC ----
# By classifier ----
ma_data3$classifier %>% 
  fct_count()

ml_model <- c("GMM", "LDA", "Logistic")

rei_class <- 
  ma_data3 %>% 
  filter(!classifier %in% ml_model) %>% 
  nest(-classifier) %>% 
  mutate(rei_fit = map(data, ~reitsma(data = .)),
         results = map(rei_fit, summary)
  ) 

rei_class %>% 
  pluck(4)

### Table of AUC ----
rei_class$results %>% 
  pluck(17)

class_AUC <- vector("list", 0)
for (i in seq_along(rei_class$results)) {
  class_AUC[[i]] <- rei_class$results[[i]]$AUC
}

class_AIC <- vector("list", 0)
for (i in seq_along(rei_class$results)) {
  class_AIC[[i]] <- rei_class$results[[1]]$AIC
}

class_BIC <- vector("list", 0)
for (i in seq_along(rei_class$results)) {
  class_BIC[[i]] <- rei_class$results[[1]]$BIC
}

class_auc <- 
  class_AUC %>% 
  enframe() %>% 
  unnest(cols = value)

class_aic <- 
  class_AIC %>% 
  enframe() %>% 
  unnest(cols = value)

class_bic <- 
  class_BIC %>% 
  enframe() %>% 
  unnest(cols = value)

class_auc$classifier <- rep(rei_class$classifier, each = 2)
class_auc$name <- rep(c("AUC", "pAUC"), times = 6)

class_auc_wide <- 
  class_auc %>% 
  pivot_wider(names_from = name, values_from = value) 

class_metrics <- 
  data.frame(class_auc_wide, AIC = class_aic$value, BIC = class_bic$value) 
class_metrics


### compare SROC ----
#plot one model
plot(rei_class$rei_fit[[1]], lty = 3, pch = 1, main = "Classifier")
#color
Color <- c("black", "red", "green", "blue", "orange", "purple")
#add lines
walk2(rei_class %>% 
        pluck(3), Color, ~ lines(sroc(.x), col = .y))  
#list for ROCellipse
List <- list(fit = rei_class %>% pluck(3),
             col = Color,
             pch = 1:6)
#ROCellipse
pwalk(List, function(fit, col, pch) ROCellipse(fit, col = col, pch = pch, add = TRUE, lty = 3))
#legend
legend("bottomright", paste0(rei_class$classifier, "(AUC=", class_metrics$AUC %>% 
                               as.numeric() %>% 
                               round(digits = 3), 
                             ")"), 
       pch = 1:6, col = Color)

### Compare AUC ----
class_nn <- ma_data3 %>% filter(classifier == "NN")
class_dl <- ma_data3 %>% filter(classifier == "DL")
class_tree <- ma_data3 %>% filter(classifier == "Tree_based")
class_knn <- ma_data3 %>% filter(classifier == "KNN")
class_svm <- ma_data3 %>% filter(classifier == "SVM")
class_bayes <- ma_data3 %>% filter(classifier == "Bayes_based")

set.seed(2021)
class_comp1 <- AUC_comparison(class_nn$TP,class_nn$FP,class_nn$FN,class_nn$TN,
                              class_dl$TP,class_dl$FP,class_dl$FN,class_dl$TN,
                              B=10000) #not converged
class_comp2 <- AUC_comparison(class_nn$TP,class_nn$FP,class_nn$FN,class_nn$TN,
                              class_tree$TP,class_tree$FP,class_tree$FN,class_tree$TN)
class_comp3 <- AUC_comparison(class_nn$TP,class_nn$FP,class_nn$FN,class_nn$TN,
                              class_knn$TP,class_knn$FP,class_knn$FN,class_knn$TN)
class_comp4 <- AUC_comparison(class_nn$TP,class_nn$FP,class_nn$FN,class_nn$TN,
                              class_svm$TP,class_svm$FP,class_svm$FN,class_svm$TN)
class_comp5 <- AUC_comparison(class_nn$TP,class_nn$FP,class_nn$FN,class_nn$TN,
                              class_bayes$TP,class_bayes$FP,class_bayes$FN,class_bayes$TN)

class_comp6 <- AUC_comparison(class_dl$TP,class_dl$FP,class_dl$FN,class_dl$TN,
                              class_tree$TP,class_tree$FP,class_tree$FN,class_tree$TN)
class_comp7 <- AUC_comparison(class_dl$TP,class_dl$FP,class_dl$FN,class_dl$TN,
                              class_knn$TP,class_knn$FP,class_knn$FN,class_knn$TN,
                              B=10000)
class_comp8 <- AUC_comparison(class_dl$TP,class_dl$FP,class_dl$FN,class_dl$TN,
                              class_svm$TP,class_svm$FP,class_svm$FN,class_svm$TN,
                              B=10000)
class_comp9 <- AUC_comparison(class_dl$TP,class_dl$FP,class_dl$FN,class_dl$TN,
                              class_bayes$TP,class_bayes$FP,class_bayes$FN,class_bayes$TN,
                              B=10000)

class_comp10 <- AUC_comparison(class_tree$TP,class_tree$FP,class_tree$FN,class_tree$TN,
                               class_knn$TP,class_knn$FP,class_knn$FN,class_knn$TN)
class_comp11 <- AUC_comparison(class_tree$TP,class_tree$FP,class_tree$FN,class_tree$TN,
                               class_svm$TP,class_svm$FP,class_svm$FN,class_svm$TN)
class_comp12 <- AUC_comparison(class_tree$TP,class_tree$FP,class_tree$FN,class_tree$TN,
                               class_bayes$TP,class_bayes$FP,class_bayes$FN,class_bayes$TN)

class_comp13 <- AUC_comparison(class_knn$TP,class_knn$FP,class_knn$FN,class_knn$TN,
                               class_svm$TP,class_svm$FP,class_svm$FN,class_svm$TN)
class_comp14 <- AUC_comparison(class_knn$TP,class_knn$FP,class_knn$FN,class_knn$TN,
                               class_bayes$TP,class_bayes$FP,class_bayes$FN,class_bayes$TN)

class_comp15 <- AUC_comparison(class_svm$TP,class_svm$FP,class_svm$FN,class_svm$TN,
                               class_bayes$TP,class_bayes$FP,class_bayes$FN,class_bayes$TN)

# By country ----
ma_data3$country %>% 
  fct_count()

country_model <- c("US+UK")

rei_country <- 
  ma_data3 %>% 
  mutate(country = fct_recode(country, 
                              US = "us", 
                              UK = "uk",
                              Others = "others", 
                              "US+UK" = "us+uk")) %>% 
  filter(!country %in% country_model) %>% 
  nest(-country) %>% 
  mutate(rei_fit = map(data, ~reitsma(data = .)),
         results = map(rei_fit, summary)
  ) 

rei_country %>% 
  pluck(4)

## Table of AUC ----
rei_country$results %>% 
  pluck(17)

country_AUC <- vector("list", 0)
for (i in seq_along(rei_country$results)) {
  country_AUC[[i]] <- rei_country$results[[i]]$AUC
}

country_AIC <- vector("list", 0)
for (i in seq_along(rei_country$results)) {
  country_AIC[[i]] <- rei_country$results[[1]]$AIC
}

country_BIC <- vector("list", 0)
for (i in seq_along(rei_country$results)) {
  country_BIC[[i]] <- rei_country$results[[1]]$BIC
}

country_auc <- 
  country_AUC %>% 
  enframe() %>% 
  unnest(cols = value)

country_aic <- 
  country_AIC %>% 
  enframe() %>% 
  unnest(cols = value)

country_bic <- 
  country_BIC %>% 
  enframe() %>% 
  unnest(cols = value)

country_auc$country <- rep(rei_country$country, each = 2)
country_auc$name <- rep(c("AUC", "pAUC"), times = 3)

country_auc_wide <- 
  country_auc %>% 
  pivot_wider(names_from = name, values_from = value) 

country_metrics <- 
  data.frame(country_auc_wide, AIC = country_aic$value, BIC = country_bic$value) 
country_metrics

## compare SROC ----
#plot one model
plot(rei_country$rei_fit[[1]], lty = 3, pch = 1, main = "Country")
#color
Color <- c("black", "red", "green")
#add lines
walk2(rei_country %>% 
        pluck(3), Color, ~ lines(sroc(.x), col = .y))  
#list for ROCellipse
List <- list(fit = rei_country %>% pluck(3),
             col = Color,
             pch = 1:3)
#ROCellipse
pwalk(List, function(fit, col, pch) ROCellipse(fit, col = col, pch = pch, add = TRUE, lty = 3))
#legend
legend("bottomright", paste0(rei_country$country, "(AUC=", country_metrics$AUC %>% 
                               as.numeric() %>% 
                               round(digits = 2), 
                             ")"), 
       pch = 1:3, col = Color)

## Compare AUC ----
country_us <- ma_data3 %>% filter(country == "us")
country_uk <- ma_data3 %>% filter(country == "uk")
country_others <- ma_data3 %>% filter(country == "others")

set.seed(2021)
country_comp1 <- AUC_comparison(country_us$TP,country_us$FP,country_us$FN,country_us$TN,
                                country_uk$TP,country_uk$FP,country_uk$FN,country_uk$TN)
country_comp2 <- AUC_comparison(country_us$TP,country_us$FP,country_us$FN,country_us$TN,
                                country_others$TP,country_others$FP,country_others$FN,country_others$TN)
country_comp3 <- AUC_comparison(country_uk$TP,country_uk$FP,country_uk$FN,country_uk$TN,
                                country_others$TP,country_others$FP,country_others$FN,country_others$TN)

# By source ----
ma_data3$source2 %>% 
  fct_count()

source_model <- c("IN-BREAST", "na", "MMD", "DDSM+MIAS")

rei_source <- 
  ma_data3 %>% 
  mutate(source2 = fct_recode(source2, "Primary data" = "original_data")) %>% 
  filter(!source2 %in% source_model) %>% 
  nest(-source2) %>% 
  mutate(rei_fit = map(data, ~reitsma(data = .)),
         results = map(rei_fit, summary)
  ) 

rei_source %>% 
  pluck(4)

## Table of AUC ----
rei_source$results %>% 
  pluck(17)

source_AUC <- vector("list", 0)
for (i in seq_along(rei_source$results)) {
  source_AUC[[i]] <- rei_source$results[[i]]$AUC
}

source_AIC <- vector("list", 0)
for (i in seq_along(rei_source$results)) {
  source_AIC[[i]] <- rei_source$results[[1]]$AIC
}

source_BIC <- vector("list", 0)
for (i in seq_along(rei_source$results)) {
  source_BIC[[i]] <- rei_source$results[[1]]$BIC
}

source_auc <- 
  source_AUC %>% 
  enframe() %>% 
  unnest(cols = value)

source_aic <- 
  source_AIC %>% 
  enframe() %>% 
  unnest(cols = value)

source_bic <- 
  source_BIC %>% 
  enframe() %>% 
  unnest(cols = value)

source_auc$source <- rep(rei_source$source2, each = 2)
source_auc$name <- rep(c("AUC", "pAUC"), times = 3)

source_auc_wide <- 
  source_auc %>% 
  pivot_wider(names_from = name, values_from = value) 

source_metrics <- 
  data.frame(source_auc_wide, AIC = source_aic$value, BIC = source_bic$value) 
source_metrics

## compare SROC ----
#plot one model
plot(rei_source$rei_fit[[1]], lty = 3, pch = 1, main = "Source")
#color
Color <- c("black", "red", "green")
#add lines
walk2(rei_source %>% 
        pluck(3), Color, ~ lines(sroc(.x), col = .y))  
#list for ROCellipse
List <- list(fit = rei_source %>% pluck(3),
             col = Color,
             pch = 1:3)
#ROCellipse
pwalk(List, function(fit, col, pch) ROCellipse(fit, col = col, pch = pch, add = TRUE, lty = 3))
#legend
legend("bottomright", paste0(rei_source$source2, "(AUC=", source_metrics$AUC %>% 
                               as.numeric() %>% 
                               round(digits = 3), 
                             ")"), 
       pch = 1:3, col = Color)

## Compare AUC ----
source_ori <- ma_data3 %>% filter(source2 == "original_data")
source_ddsm <- ma_data3 %>% filter(source2 == "DDSM")
source_mias <- ma_data3 %>% filter(source2 == "MIAS")

set.seed(2021)
source_comp1 <- AUC_comparison(source_ori$TP,source_ori$FP,source_ori$FN,source_ori$TN,
                               source_ddsm$TP,source_ddsm$FP,source_ddsm$FN,source_ddsm$TN,
                               B = 10000)
source_comp2 <- AUC_comparison(source_ori$TP,source_ori$FP,source_ori$FN,source_ori$TN,
                               source_mias$TP,source_mias$FP,source_mias$FN,source_mias$TN)
source_comp3 <- AUC_comparison(source_ddsm$TP,source_ddsm$FP,source_ddsm$FN,source_ddsm$TN,
                               source_mias$TP,source_mias$FP,source_mias$FN,source_mias$TN,
                               B = 10000)

AUC_boot(source_ddsm$TP,source_ddsm$FP,source_ddsm$FN,source_ddsm$TN, B = 20000)


# Influential diagnostics ----
library(dmetatools)
library(doParallel)

# This will take 4 and 1/2 hours
# Create a cluster object and then register: 
cl <- makePSOCKcluster(4)
registerDoParallel(cl)

set.seed(2021)
infl_diag <- AUC_IF(ma_data3$TP, ma_data3$FP, ma_data3$FN, ma_data3$TN)

stopCluster(cl)

infl_diag
