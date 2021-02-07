#
## WAKEUP
# PURPOSE
# Binary logistic regression model to predict association with primary Endpoint (MRS "yes", "no")
#
# INPUT
# wakeupICH_FINAL_TABLE
#
# MJ 05.04.20

# packages
library(tidyverse)
library(mlbench)
library(MASS)
library(pROC)

# get data
# setwd("~/projects/WakeUp/neu")
# dataOrg <- read.csv("../wakeup_master_csv.csv", sep = ";", dec=",",stringsAsFactors = TRUE) # put directory 
# data <- dataOrg

# remove first column
# dataOrg <- dataOrg %>%
# select(-X)

# create new variables for age and NIHSS as binary variables (mrs0_1, Age >60, and NIHSS >10)
#data$AGEover60 <- cut(data$AGE, breaks=c(0,60,100),
#                  labels = c("no","yes") # new variable AGoverE60 ( = AGE > 60)

#data$NIHSSover10 <- cut(data$NIHSSSCORE, breaks=c(-1,10,40),
#                    labels = c("no","yes") # new variable NIHSSover10 (= NIHSS > 10)

#data$mrs0_1 <- cut(data$MRSSCOREd90, breaks= c(-1, 1, 6), labels = c("yes", "no"))

data <- sae
data$SAE <- as.factor(sae$SAE)
data$MRSSCOREd90 <- as.numeric(sae$MRSSCOREd90)
data$mrs0_1 <- sae$MRSSCOREd90 <= 1
data$lesion_location <- fct_collapse(data$lesion_location, vertebrobasilar = c("brainstem", "cerebellum", "pca"), anterior_circulation = "mcaOrAca", multiple = "multiple")

# put targets to first column
data <- data %>%
dplyr::select(USUBJID, MRSSCOREd0, MRSSCOREd90, mrs0_1, NIHSSSCORE, AGE, everything())

# 1. Perform analysis of primary endpoint for SAE

# create a new data frame of relevant modeling variables
newdata <- data[,c("SAE","AGE", "SEX", "mrs0_1","NIHSSSCORE", "codetrt", "stroke_volume", "LVO_V0", "lesion_location", "Arterial_hypertension", "Atrial_fibrillation", "Diabetes_mellitus_type_II")]
newdata <- newdata %>% 
  rename(age= AGE, treatment= codetrt, LVO = LVO_V0, NIHSS_baseline = NIHSSSCORE ) %>% 
  mutate(age= age/10)
#newdata$lesion_side <- ifelse(newdata$lesion_side == "right", 1,0)
#newdata$lesion_side <- as.factor(newdata$lesion_side)
newdata$mrs0_1 <- sae$MRSSCOREd90 <= 1
newdata$SEX <-as.factor(newdata$SEX)
newdata$Arterial_hypertension <- as.factor(newdata$Arterial_hypertension)
newdata$Atrial_fibrillation <- as.factor(newdata$Atrial_fibrillation)
newdata$Diabetes_mellitus_type_II <- as.factor(newdata$Diabetes_mellitus_type_II)

newdata$Arterial_hypertension[newdata$Arterial_hypertension ==2] <- NA
newdata$Atrial_fibrillation[newdata$Atrial_fibrillation ==2] <- NA
newdata$Diabetes_mellitus_type_II[newdata$Diabetes_mellitus_type_II ==2] <- NA

head(newdata)
str(newdata)

# Implementation of Logistic Regression to predict the binary outcome, exclude NA

logit_1 <- glm(mrs0_1 ~ SAE+ age + SEX+ NIHSS_baseline + LVO + stroke_volume + Arterial_hypertension + Atrial_fibrillation + Diabetes_mellitus_type_II, data= newdata, family= binomial, na.action=na.exclude)
summary(logit_1)
broom::tidy(logit_1, exponentiate= TRUE, conf.int = TRUE)
result_log <- summary.glm(logit_1)$coefficients

#export table
write.csv(result_log, file= "wakeup_binary_regression_R!.csv")

# Forest plot
plot_model(logit_1, sort.est = TRUE, axis.label= "", title= "Favorable Outcome (mRS 0-1)")

#association between LVO and stroke volume
logit_2 <- glm(LVO ~ stroke_volume, data= newdata, family = binomial, na.action= na.exclude)
broom::tidy(logit_2, exponentiate= TRUE, conf.int = TRUE)
plot_model (logit_2)

#model with cut NIHSS <10, Age > 60
logit_2 <- glm(mrs0_1 ~ I(AGE<=60) + I(NIHSSSCORE<=10) + codetrt + SAE + stroke_volume, data = newdata, family = binomial, na.action = na.exclude)
summary(logit_2)

#
#
#
#
#

# binary outcome for singel organ class "cardiac disorders", "nervoussystem disorders" and "infectious disorders"
newdata2$MRSSCOREd90 <- as.numeric(sae$MRSSCOREd90)
newdata2$mrs0_1 <- sae$MRSSCOREd90 <= 1

# nervous system disorders
model_Neuro <- glm(mrs0_1 ~ `Nervous system disorders`+  AGE + SEX + stroke_volume + Arterial_hypertension + Atrial_fibrillation + Diabetes_mellitus_type_II + NIHSSSCORE + LVO_V0, data= newdata2, family= binomial, na.action=na.exclude)
summary(model_Neuro)
broom::tidy(model_Neuro, conf.int=T, exponentiate=T)
plot_model(model_Neuro)

# cardiac disoders
model_Cardiac <- glm(mrs0_1 ~ `Cardiac disorders`+ AGE + SEX + stroke_volume + Arterial_hypertension + Atrial_fibrillation + Diabetes_mellitus_type_II + NIHSSSCORE + LVO_V0, data= newdata2, family= binomial, na.action=na.exclude)
summary(model_Cardiac)
broom::tidy(model_Cardiac, conf.int=T, exponentiate=T)
plot_model(model_Neuro)

# infectious disorders
model_Infect <- glm(mrs0_1 ~ `Infections and infestations`+  AGE + SEX + stroke_volume + Arterial_hypertension + Atrial_fibrillation + Diabetes_mellitus_type_II + NIHSSSCORE + LVO_V0, data= newdata2, family= binomial, na.action=na.exclude)
summary(model_Infect)
broom::tidy(model_Infect, conf.int=T, exponentiate=T)

#####################################################
#prediction of occurrence of SAE
newdata$SEX <- as.factor(newdata$SEX)
logit_3 <- glm(SAE ~ age + SEX + stroke_volume + Arterial_hypertension + Atrial_fibrillation + Diabetes_mellitus_type_II + NIHSS_baseline + LVO,   data= newdata, family = binomial)
summary(logit_3)
broom::tidy(logit_3, exponentiate= TRUE, conf.int = TRUE)
write.csv(result_log, file= "prediction_SAE.csv")
# Forest plot
plot_model(logit_3, sort.est = TRUE, axis.label= "")

#logitSAE <- glm(SAE ~ ., data= newdata, family = binomial, na.action = na.exclude)
#summary(logitSAE)

# Display the odds ratio and 95% CI
exp(coefficients(logit_1)[1:6])
exp(confint(logit_1)[1:6])

# stepwise model selection to minimize the AIC value
#logit_3 <- stepAIC(logit_1)
summary(logit_1)




