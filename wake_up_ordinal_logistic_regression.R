#
#WAKE UP
#PURPOSE
# --- Ordinal Logistic Regression
#
#Packages
library(tidyverse)
library(mlbench)
library(MASS)
library(pROC)

#INPUT
data <- sae

str(data)

# clean data
data$SAE <- as.factor(data$SAE)
data$MRSSCOREd90<- as.ordered(data$MRSSCOREd90)


# create a new data frame of relevant modeling variables, for ordinal model SAE in general = newdata
newdata <- data[,c("MRSSCOREd90", "SAE","AGE", "SEX","NIHSSSCORE", "codetrt", "stroke_volume", "LVO_V0", "lesion_location", "Arterial_hypertension", "Atrial_fibrillation", "Diabetes_mellitus_type_II")]

#newdata$lesion_side <- ifelse(newdata$lesion_side == "both", 1, 0)
newdata$NIHSSSCORE <- as.numeric(newdata$NIHSSSCORE)
newdata <- newdata %>% 
  rename(age= AGE, treatment= codetrt, LVO = LVO_V0, NIHSS_baseline = NIHSSSCORE ) %>% 
  mutate(age= age/10)

#preparing pre-inlesses for analysis in DF newdata (needed for analysis SAE in general)
newdata$Arterial_hypertension <- as.factor(newdata$Arterial_hypertension)
newdata$Atrial_fibrillation <- as.factor(newdata$Atrial_fibrillation)
newdata$Diabetes_mellitus_type_II <- as.factor(newdata$Diabetes_mellitus_type_II)

newdata$Arterial_hypertension[newdata$Arterial_hypertension ==2] <- NA
newdata$Atrial_fibrillation[newdata$Atrial_fibrillation ==2] <- NA
newdata$Diabetes_mellitus_type_II[newdata$Diabetes_mellitus_type_II ==2] <- NA

newdata$Arterial_hypertension <- ifelse(newdata$Arterial_hypertension ==1, 1, 0)
newdata$Atrial_fibrillation <- ifelse(newdata$Atrial_fibrillation  ==1, 1, 0)
newdata$Diabetes_mellitus_type_II <- ifelse(newdata$Diabetes_mellitus_type_II  ==1, 1, 0)


newdata <- read.csv('newdata.csv') %>% 
  mutate(MRSSCOREd90 = factor(MRSSCOREd90, levels = seq(0,6), ordered = TRUE))

################################ Ordinal logistic regression for SAE in general#######################################
############################################################################################
model <- polr(MRSSCOREd90 ~ SAE+ age + SEX+ NIHSS_baseline + LVO + stroke_volume + Arterial_hypertension + Atrial_fibrillation + Diabetes_mellitus_type_II, data= newdata, Hess= TRUE, na.action=na.exclude)
summary(model)
broom::tidy(model, conf.int=T, exponentiate=T)

#p-value calculation
(ctable <- coef(summary(model)))
model_2 <- pnorm(abs(ctable[, "t value"]), lower.tail=FALSE) *2
(ctable <- cbind(ctable, "p value" = model_2))
summary(ctable)

### proportional odds assumptions ####
require(VGAM)
res.PO <- vglm(MRSSCOREd90 ~ SAE + age + SEX + NIHSS_baseline + LVO + stroke_volume + Arterial_hypertension + Atrial_fibrillation + Diabetes_mellitus_type_II
            , family = propodds
            , data = newdata)
res.Multi <- vglm(MRSSCOREd90 ~ SAE + age + SEX + NIHSS_baseline + LVO + stroke_volume + Arterial_hypertension + Atrial_fibrillation + Diabetes_mellitus_type_II
               , family = multinomial
               , data = newdata %>% mutate(MRSSCOREd90 = factor(MRSSCOREd90, ordered = FALSE)))

p <- pchisq(deviance(res.PO) - deviance(res.Multi)
       , df = df.residual(res.PO)-df.residual(res.Multi)
       , lower.tail=FALSE)

pred.PO <- predictvglm(res.PO, type ='link', se.fit = TRUE)
pred.Multi <- predictvglm(res.Multi, type ='response', se.fit = FALSE)

df.PO <- pred.PO$fitted.values %>% as_tibble(rownames = 'ID') %>% 
  pivot_longer(2:7, values_to = 'fit') %>% 
  bind_cols(pred.PO$se.fit %>% as_tibble() %>% pivot_longer(1:6, values_to = 'se.fit') %>% dplyr::select(-name)) %>% 
  mutate(lwr = fit - 1.96*se.fit
         , upr = fit + 1.96*se.fit
         , fit.trans = res.PO@family@linkinv(fit)[,1]
         , lwr.trans = res.PO@family@linkinv(lwr)[,1]
         , upr.trans = res.PO@family@linkinv(upr)[,1]) 
df.PO %>% 
  group_by(name) %>% 
  summarise(mean = mean(fit.trans), sd= sd(fit.trans)) %>% 
  ggplot(aes(x = name, y = mean)) +
  geom_errorbar(aes(ymin = mean - 1.96*sd, ymax = mean + 1.96*sd)) +
  geom_point() +
  theme_minimal()



df.PO.Multi <- pred.Multi %>% as_tibble(rownames = 'ID') %>% 
  setNames(c('ID', paste0('p',0:6))) %>% 
  mutate(`logitlink(P[Y>=2])` = p0
         , `logitlink(P[Y>=3])` = p0 + p1
         , `logitlink(P[Y>=4])` = p0 + p1 + p2
         , `logitlink(P[Y>=5])` = p0 + p1 + p2 + p3
         , `logitlink(P[Y>=6])` = p0 + p1 + p2 + p3 + p4
         , `logitlink(P[Y>=7])` = p0 + p1 + p2 + p3 + p4 + p5) %>% 
  dplyr::select(starts_with('logit')) %>% 
  pivot_longer(1:6, values_to = 'prob') %>% dplyr::select(-name) %>% 
  bind_cols(df.PO) 

df.PO.Multi %>% 
  mutate(delta = prob - fit.trans) %>% 
  group_by(name) %>% 
  summarise(mean = mean(delta), sd = sd(delta)) %>% 
  ggplot(aes(x = name, y = mean)) +
  geom_errorbar(aes(ymin = mean - 1.96*sd, ymax = mean + 1.96*sd)) +
  geom_point() +
  theme_minimal()

df.PO.Multi %>% 
  ggplot(aes(x = prob, y = fit.trans)) +
  geom_point()
####

########################
### ordinal logistic regression für nervous system disorders and infectious disorders
# create new dataframe for ordinal model for SAE in different organ classes = newdata2
newdata2 <- data[,c("MRSSCOREd90", "Atrial_fibrillation", "NIHSSSCORE", "SEX", "stroke_volume", "Weight", "Serum_Glucose", "Arterial_hypertension", "LVO_V0", "lesion_location", "Infections and infestations", "Nervous system disorders", "Cardiac disorders", "AGE", "Arterial_hypertension", "Atrial_fibrillation", "Diabetes_mellitus_type_II")]
newdata2$`Nervous system disorders` <- as.factor(newdata2$`Nervous system disorders`)
levels(newdata2$`Nervous system disorders`) <- c('0','1', '1', '1')
newdata2$`Infections and infestations` <- as.factor(newdata2$`Infections and infestations`)
levels(newdata2$`Infections and infestations`) <- c('0','1', '1')
newdata2$`Cardiac disorders` <- as.factor(newdata2$`Cardiac disorders`)
levels(newdata2$`Cardiac disorders`) <- c('0','1', '1')

#preparing pre-ilesses for analysis in DF newdata2
newdata2$Arterial_hypertension <- as.factor(newdata2$Arterial_hypertension)
newdata2$Atrial_fibrillation <- as.factor(newdata2$Atrial_fibrillation)
newdata2$Diabetes_mellitus_type_II <- as.factor(newdata2$Diabetes_mellitus_type_II)

newdata2$Arterial_hypertension[newdata2$Arterial_hypertension ==2] <- NA
newdata2$Atrial_fibrillation[newdata2$Atrial_fibrillation ==2] <- NA
newdata2$Diabetes_mellitus_type_II[newdata2$Diabetes_mellitus_type_II ==2] <- NA

newdata2$Arterial_hypertension <- ifelse(newdata2$Arterial_hypertension ==1, 1, 0)
newdata2$Atrial_fibrillation <- ifelse(newdata2$Atrial_fibrillation  ==1, 1, 0)
newdata2$Diabetes_mellitus_type_II <- ifelse(newdata2$Diabetes_mellitus_type_II  ==1, 1, 0)

newdata2$SEX <- as.factor(newdata2$SEX)

newdata2 <- read.csv('newdata2.csv') %>% 
  mutate(MRSSCOREd90 = factor(MRSSCOREd90, levels = seq(0,6), ordered = TRUE))

#Neurological disorders
model_Neuro <- polr(MRSSCOREd90 ~ Nervous.system.disorders+ AGE + SEX + stroke_volume + Arterial_hypertension + Atrial_fibrillation + Diabetes_mellitus_type_II + NIHSSSCORE + LVO_V0, data= newdata2, Hess = TRUE, na.action=na.exclude)
summary(model_Neuro)
broom::tidy(model_Neuro, conf.int=T, exponentiate=T)
plot_model(model_Neuro)

#Infectious disorders
model_Infect <- polr(MRSSCOREd90 ~ Infections.and.infestations+  AGE + SEX + stroke_volume + Arterial_hypertension + Atrial_fibrillation + Diabetes_mellitus_type_II + NIHSSSCORE + LVO_V0, data= newdata2, Hess = TRUE, na.action=na.exclude)
summary(model_Infect)
broom::tidy(model_Infect, conf.int=T, exponentiate=T)
plot_model(model_Infect)

#Cardiac disorders
model_Cardiac <- polr(MRSSCOREd90 ~ Cardiac.disorders + AGE + SEX + stroke_volume + Arterial_hypertension + Atrial_fibrillation + Diabetes_mellitus_type_II + NIHSSSCORE + LVO_V0, data= newdata2, Hess = TRUE, na.action=na.exclude)
summary(model_Cardiac)
broom::tidy(model_Cardiac, conf.int=T, exponentiate=T)
plot_model(model_Cardiac)

#p-value calculation --> change the model for each organ system
(ctable <- coef(summary(model_Cardiac)))
model_2 <- pnorm(abs(ctable[, "t value"]), lower.tail=FALSE) *2
(ctable <- cbind(ctable, "p value" = model_2))
summary(ctable)


### PO assumption ####

res.PO <- vglm(MRSSCOREd90 ~ Cardiac.disorders + AGE + SEX + stroke_volume + Arterial_hypertension + Atrial_fibrillation + Diabetes_mellitus_type_II + NIHSSSCORE + LVO_V0
               , family = propodds
               , data = newdata2)
res.Multi <- vglm(MRSSCOREd90 ~ Cardiac.disorders + AGE + SEX + stroke_volume + Arterial_hypertension + Atrial_fibrillation + Diabetes_mellitus_type_II + NIHSSSCORE + LVO_V0
                  , family = multinomial
                  , data = newdata2 %>% mutate(MRSSCOREd90 = factor(MRSSCOREd90, ordered = FALSE)))

p <- pchisq(deviance(res.PO)-deviance(res.Multi)
            , df=df.residual(res.PO)-df.residual(res.Multi)
            , lower.tail=FALSE)
p

#######################

##### Combination of multiple SAEs ################

model_SAE<- polr(MRSSCOREd90 ~ Cardiac.disorders + Infections.and.infestations + AGE + SEX + stroke_volume + Arterial_hypertension + Atrial_fibrillation + Diabetes_mellitus_type_II + NIHSSSCORE + LVO_V0, data= newdata2, Hess = TRUE, na.action=na.exclude)
broom::tidy(model_SAE, conf.int = TRUE)
pchisq(deviance(model_Infect) - deviance(model_SAE)  # model comparison
       , df = df.residual(model_Infect) - df.residual(model_SAE)
       , lower.tail=FALSE)###################################################


m.cardiac.vglm <- vglm(MRSSCOREd90 ~ Cardiac.disorders + AGE + SEX + stroke_volume + Arterial_hypertension + Atrial_fibrillation + Diabetes_mellitus_type_II + NIHSSSCORE + LVO_V0
               , family = propodds
               , data = newdata2)
m.infect.vglm <- vglm(MRSSCOREd90 ~ Infections.and.infestations + AGE + SEX + stroke_volume + Arterial_hypertension + Atrial_fibrillation + Diabetes_mellitus_type_II + NIHSSSCORE + LVO_V0
                       , family = propodds
                       , data = newdata2)

m.SAE.vglm <- vglm(MRSSCOREd90 ~ Cardiac.disorders  + Infections.and.infestations + AGE + SEX + stroke_volume + Arterial_hypertension + Atrial_fibrillation + Diabetes_mellitus_type_II + NIHSSSCORE + LVO_V0
                      , family = propodds
                      , data = newdata2)

lrtest(m.SAE.vglm, m.cardiac.vglm)


#####################p
# ordinal regression model for prediction of SAE
#pred <- (predict(model, newdata[1:5,], type= "prob"))
#pred
summary(model)
broom::tidy(model, conf.int=T, exponentiate=T)

plot_model(model, sort.est = TRUE, axis.label = "")
result_ord_mrs0_1 <- summary.lm(model)$coefficients

write.csv(result_ord_mrs0_1, file= "wake_up_back_regression_mrs0_1_02.csv")



