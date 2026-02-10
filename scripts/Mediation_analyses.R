## Mediation analysis

install.packages("mediation")
library(mediation)
vignette("mediation")
names(MRS_full)
citation("mediation")


################################# DATA PREPARATION #####################
# Select raw columns
mediation_combined <- MRS_full[, c("m_m_precuneus", "m_m_acc", 
                                   "hipp_mean", "cortical_thickness_adsignature_dickson", 
                                   "hipp_mean_act", "activation_parietal_sup_l", 
                                   "sex", "age_spectro")] |> na.omit()
# Create composite "mean" variables
mediation_combined$mean_structure <- rowMeans(scale(mediation_combined[, c("hipp_mean", "cortical_thickness_adsignature_dickson")]), na.rm = TRUE)
mediation_combined$mean_activation <- rowMeans(scale(mediation_combined[, c("hipp_mean_act", "activation_parietal_sup_l")]), na.rm = TRUE)
mediation_combined$mean_glu <- rowMeans(scale(mediation_combined[, c("m_m_precuneus", "m_m_acc")]), na.rm = TRUE)



###############################  GLOBAL COMBINED MODEL ###########################
# Mean Structure -> Mean Glu -> Mean Activation
# M ~ X
model_m_comb <- lm(mean_glu ~ mean_structure, data = mediation_combined)
summary(model_m_comb)

# Y ~ X + M
model_y_comb <- lm(mean_activation ~ mean_structure + mean_glu, data = mediation_combined)
summary(model_y_comb)

# Mediate
med_comb <- mediate(model_m_comb, model_y_comb, 
                    treat = "mean_structure", mediator = "mean_glu", 
                    boot = TRUE, sims = 5000)
summary(med_comb)
plot(med_comb)


#### With robmed package ####
library(robmed)


rob_mediation <- fit_mediation(mediation_combined,
  x = "mean_structure",
  y = "mean_activation",
  m = "mean_glu")
summary(rob_mediation)

rob_test <- test_mediation(rob_mediation)
summary(rob_test)



######## With WRS2 package ##########
library(WRS2)
x <- mediation_combined$mean_structure
y <- mediation_combined$mean_activation
med <- mediation_combined$mean_glu
zmediate_results <- ZYmediate(x, y, med, nboot = 5000, alpha = 0.05, kappa = 0.05)
zmediate_results





################### DECOMPOSING ACTIVATION (Y) ######################
# Keeping Structure and Glu Combined, splitting Activation
## 2.1 Parietal Superior Activation
# M ~ X
model_m_par <- lm(mean_glu ~ mean_structure, data = mediation_combined)
summary(model_m_par)
# Y ~ X + M
model_y_par <- lm(activation_parietal_sup_l ~ mean_structure + mean_glu, data = mediation_combined)
summary(model_y_par)

# Mediate
med_par <- mediate(model_m_par, model_y_par, 
                   treat = "mean_structure", mediator = "mean_glu", 
                   boot = TRUE, sims = 5000)
summary(med_par)
plot(med_par)

## 2.2 Hippocampal Activation
# M ~ X
model_m_hipp_act <- lm(mean_glu ~ mean_structure, data = mediation_combined)
summary(model_m_hipp_act)
# Y ~ X + M
model_y_hipp_act <- lm(hipp_mean_act ~ mean_structure + mean_glu, data = mediation_combined)
summary(model_y_hipp_act)

# Mediate
med_hipp_act <- mediate(model_m_hipp_act, model_y_hipp_act, 
                        treat = "mean_structure", mediator = "mean_glu", 
                        boot = TRUE, sims = 5000)
summary(med_hipp_act)
plot(med_hipp_act)



################# 3. DECOMPOSING GLUTAMATE (M)#############
# Keeping Structure and Activation Combined, splitting Glu

## 3.1 Precuneus Glu
# M ~ X
model_m_prec <- lm(m_m_precuneus ~ mean_structure, data = mediation_combined)
summary(model_m_prec)
# Y ~ X + M
model_y_prec <- lm(mean_activation ~ mean_structure + m_m_precuneus, data = mediation_combined)
summary(model_y_prec)
# Mediate
med_glu_prec <- mediate(model_m_prec, model_y_prec, 
                        treat = "mean_structure", mediator = "m_m_precuneus", 
                        boot = TRUE, sims = 5000)
summary(med_glu_prec)
plot(med_glu_prec)

## 3.2 ACC Glu
# M ~ X
model_m_acc <- lm(m_m_acc ~ mean_structure, data = mediation_combined)
summary(model_m_acc)

# Y ~ X + M
model_y_acc <- lm(mean_activation ~ mean_structure + m_m_acc, data = mediation_combined)
summary(model_y_acc)

# Mediate
med_glu_acc <- mediate(model_m_acc, model_y_acc, 
                       treat = "mean_structure", mediator = "m_m_acc", 
                       boot = TRUE, sims = 5000)
summary(med_glu_acc)
plot(med_glu_acc)


################### DECOMPOSING STRUCTURE (X) #################
# Keeping Glu and Activation Combined, splitting Structure

## 4.1 Hippocampal Volume
# M ~ X
model_m_hipp_vol <- lm(mean_glu ~ hipp_mean, data = mediation_combined)
summary(model_m_hipp_vol)
# Y ~ X + M
model_y_hipp_vol <- lm(mean_activation ~ hipp_mean + mean_glu, data = mediation_combined)
summary(model_y_hipp_vol)

# Mediate
med_struc_hipp <- mediate(model_m_hipp_vol, model_y_hipp_vol, 
                          treat = "hipp_mean", mediator = "mean_glu", 
                          boot = TRUE, sims = 5000)
summary(med_struc_hipp)
plot(med_struc_hipp)

## 4.2 Cortical Thickness (AD Signature)
# M ~ X
model_m_cort <- lm(mean_glu ~ cortical_thickness_adsignature_dickson, data = mediation_combined)
summary(model_m_cort)
# Y ~ X + M
model_y_cort <- lm(mean_activation ~ cortical_thickness_adsignature_dickson + mean_glu, data = mediation_combined)
summary(model_y_cort)
# Mediate
med_struc_cort <- mediate(model_m_cort, model_y_cort, 
                          treat = "cortical_thickness_adsignature_dickson", mediator = "mean_glu", 
                          boot = TRUE, sims = 5000)
summary(med_struc_cort)
plot(med_struc_cort)






################# 4. DECOMPOSING STRUCTURE (X) #############
# Keeping Glu and Activation Combined, splitting Structure

## 4.1 Hippocampal Volume (X)
# M ~ X
model_m_hipp_vol <- lm(mean_glu ~ hipp_mean, data = mediation_combined)
summary(model_m_hipp_vol)

# Y ~ X + M
model_y_hipp_vol <- lm(mean_activation ~ hipp_mean + mean_glu, data = mediation_combined)
summary(model_y_hipp_vol)

# Mediate
med_struc_hipp <- mediate(model_m_hipp_vol, model_y_hipp_vol, 
                          treat = "hipp_mean", mediator = "mean_glu", 
                          boot = TRUE, sims = 5000)
summary(med_struc_hipp)
plot(med_struc_hipp)

## 4.2 Cortical Thickness (X)
# M ~ X
model_m_cort <- lm(mean_glu ~ cortical_thickness_adsignature_dickson, data = mediation_combined)
summary(model_m_cort)

# Y ~ X + M
model_y_cort <- lm(mean_activation ~ cortical_thickness_adsignature_dickson + mean_glu, data = mediation_combined)
summary(model_y_cort)

# Mediate
med_struc_cort <- mediate(model_m_cort, model_y_cort, 
                          treat = "cortical_thickness_adsignature_dickson", mediator = "mean_glu", 
                          boot = TRUE, sims = 5000)
summary(med_struc_cort)
plot(med_struc_cort)





################# 5. FULLY GRANULAR (ALL INDIVIDUAL VARIABLES) #############
# Testing specific paths: Specific Structure -> Specific Glu -> Specific Activation

# --- GROUP A: HIPPOCAMPAL VOLUME AS PREDICTOR (X = hipp_mean) ---

## 5.1 Hipp Volume -> Precuneus Glu -> Parietal Activation
model_m_g1 <- lm(m_m_precuneus ~ hipp_mean, data = mediation_combined)
summary(model_m_g1)
model_y_g1 <- lm(activation_parietal_sup_l ~ hipp_mean + m_m_precuneus, data = mediation_combined)
summary(model_y_g1)

med_g1 <- mediate(model_m_g1, model_y_g1, treat = "hipp_mean", mediator = "m_m_precuneus", boot = TRUE, sims = 5000)
summary(med_g1)
plot(med_g1)

## 5.2 Hipp Volume -> Precuneus Glu -> Hippocampal Activation
model_m_g2 <- lm(m_m_precuneus ~ hipp_mean, data = mediation_combined)
summary(model_m_g2)
model_y_g2 <- lm(hipp_mean_act ~ hipp_mean + m_m_precuneus, data = mediation_combined)
summary(model_y_g2)

med_g2 <- mediate(model_m_g2, model_y_g2, treat = "hipp_mean", mediator = "m_m_precuneus", boot = TRUE, sims = 5000)
summary(med_g2)
plot(med_g2)

## 5.3 Hipp Volume -> ACC Glu -> Parietal Activation
model_m_g3 <- lm(m_m_acc ~ hipp_mean, data = mediation_combined)
summary(model_m_g3)
model_y_g3 <- lm(activation_parietal_sup_l ~ hipp_mean + m_m_acc, data = mediation_combined)
summary(model_y_g3)

med_g3 <- mediate(model_m_g3, model_y_g3, treat = "hipp_mean", mediator = "m_m_acc", boot = TRUE, sims = 5000)
summary(med_g3)
plot(med_g3)

## 5.4 Hipp Volume -> ACC Glu -> Hippocampal Activation
model_m_g4 <- lm(m_m_acc ~ hipp_mean, data = mediation_combined)
summary(model_m_g4)
model_y_g4 <- lm(hipp_mean_act ~ hipp_mean + m_m_acc, data = mediation_combined)
summary(model_y_g4)

med_g4 <- mediate(model_m_g4, model_y_g4, treat = "hipp_mean", mediator = "m_m_acc", boot = TRUE, sims = 5000)
summary(med_g4)
plot(med_g4)


# --- GROUP B: CORTICAL THICKNESS AS PREDICTOR (X = cortical_thickness) ---

## 5.5 Cort Thickness -> Precuneus Glu -> Parietal Activation
model_m_g5 <- lm(m_m_precuneus ~ cortical_thickness_adsignature_dickson, data = mediation_combined)
summary(model_m_g5)
model_y_g5 <- lm(activation_parietal_sup_l ~ cortical_thickness_adsignature_dickson + m_m_precuneus, data = mediation_combined)
summary(model_y_g5)

med_g5 <- mediate(model_m_g5, model_y_g5, treat = "cortical_thickness_adsignature_dickson", mediator = "m_m_precuneus", boot = TRUE, sims = 5000)
summary(med_g5)
plot(med_g5)

## 5.6 Cort Thickness -> Precuneus Glu -> Hippocampal Activation
model_m_g6 <- lm(m_m_precuneus ~ cortical_thickness_adsignature_dickson, data = mediation_combined)
summary(model_m_g6)
model_y_g6 <- lm(hipp_mean_act ~ cortical_thickness_adsignature_dickson + m_m_precuneus, data = mediation_combined)
summary(model_y_g6)

med_g6 <- mediate(model_m_g6, model_y_g6, treat = "cortical_thickness_adsignature_dickson", mediator = "m_m_precuneus", boot = TRUE, sims = 5000)
summary(med_g6)
plot(med_g6)

## 5.7 Cort Thickness -> ACC Glu -> Parietal Activation
model_m_g7 <- lm(m_m_acc ~ cortical_thickness_adsignature_dickson, data = mediation_combined)
summary(model_m_g7)
model_y_g7 <- lm(activation_parietal_sup_l ~ cortical_thickness_adsignature_dickson + m_m_acc, data = mediation_combined)
summary(model_y_g7)

med_g7 <- mediate(model_m_g7, model_y_g7, treat = "cortical_thickness_adsignature_dickson", mediator = "m_m_acc", boot = TRUE, sims = 5000)
summary(med_g7)
plot(med_g7)

## 5.8 Cort Thickness -> ACC Glu -> Hippocampal Activation
model_m_g8 <- lm(m_m_acc ~ cortical_thickness_adsignature_dickson, data = mediation_combined)
summary(model_m_g8)
model_y_g8 <- lm(hipp_mean_act ~ cortical_thickness_adsignature_dickson + m_m_acc, data = mediation_combined)
summary(model_y_g8)

med_g8 <- mediate(model_m_g8, model_y_g8, treat = "cortical_thickness_adsignature_dickson", mediator = "m_m_acc", boot = TRUE, sims = 5000)
summary(med_g8)
plot(med_g8)












