library(readxl)   
library(janitor)  
library(jmv)
library(effectsize)
library(interactions)
library(emmeans)

## Glutamate Relates to Structural and Functional markers of Disease Severity in early Alzheimer’s disease  ######
## Arsenii Prozorov 
#ANALYSES PRÉLIMINAIRES :
#Création d’une banque de données

import_data_arsenii_20250609_hippR <- read_excel("C:/Users/okkam/Desktop/labo/article 1/rencontre_Sylvie_20251017/import_data_arsenii_20250609_hippR.xlsx")
MRS_full <- import_data_arsenii_20250609_hippR

# Clean the column name
MRS_full <- janitor::clean_names(MRS_full)
names(MRS_full)
# Convertir  en numérique
MRS_full$education <- as.numeric(MRS_full$education)
MRS_full$moca_corr_spectro <- as.numeric(MRS_full$moca_corr_spectro)
MRS_full$dprime_hit_fa <- as.numeric(MRS_full$dprime_hit_fa)
MRS_full$associative_memory_performance <- as.numeric(MRS_full$associative_memory_performance)
MRS_full$moca_corr_spectro <- as.numeric(MRS_full$moca_corr_spectro)
# Regrouper les groupes 'HC' et 'SCD' dans un seul groupe 'SCD'
MRS_full$diagnostic_nick <- ifelse(MRS_full$diagnostic_nick %in% c("HC", "SCD"), "HC", MRS_full$diagnostic_nick)

# Convertir diagnostic_nick en facteur avec l’ordre HC, SCD+, MCI
MRS_full$diagnostic_nick <- factor(MRS_full$diagnostic_nick, levels = c("HC", "SCD+", "MCI"))

# Convertir sex en facteur
MRS_full$sex <- as.factor(MRS_full$sex)

# Creer Hipp mean
MRS_full$hipp_mean <- (MRS_full$hip_l_nor_icv + MRS_full$righ_hip_vol)/2
MRS_full$hipp_mean_act <- (MRS_full$activation_hippocampus_l + MRS_full$activation_hippocampus_r)/2



##################  Winzorising ###################
# Define the winsorize function
winsorize_iqr <- function(x, iqr_multiplier) {
  qnt <- quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
  iqr <- qnt[2] - qnt[1]
  lower <- qnt[1] - iqr_multiplier * iqr
  upper <- qnt[2] + iqr_multiplier * iqr
  x[x < lower] <- lower
  x[x > upper] <- upper
  return(x)}

# Apply winsorization for 1.5 IQR variables
MRS_full$memoria_libre_correcte <- winsorize_iqr(MRS_full$memoria_libre_correcte, 1.5)
MRS_full$face_name_rappel_differe_spectro <- winsorize_iqr(MRS_full$face_name_rappel_differe_spectro, 1.5)
MRS_full$associative_memory_performance <- winsorize_iqr(MRS_full$associative_memory_performance, 1.5)
MRS_full$activation_temporal_inf_r <- winsorize_iqr(MRS_full$activation_temporal_inf_r, 1.5)
MRS_full$activation_parietal_sup_l <- winsorize_iqr(MRS_full$activation_parietal_sup_l, 1.5)
MRS_full$activation_hippocampus_l <- winsorize_iqr(MRS_full$activation_hippocampus_l, 1.5)
MRS_full$hipp_mean_act <- winsorize_iqr(MRS_full$hipp_mean_act, 1.5)
MRS_full$moca_corr_spectro <- winsorize_iqr(MRS_full$moca_corr_spectro, 1.5)
MRS_full$dprime_hit_fa <- winsorize_iqr(MRS_full$dprime_hit_fa, 1.5)
MRS_full$associative_memory_performance_sylvie <- winsorize_iqr(MRS_full$associative_memory_performance_sylvie, 1.5)


# Apply winsorization for 1.7 IQR variables
MRS_full$hip_l_nor_icv <- winsorize_iqr(MRS_full$hip_l_nor_icv, 1.7)
MRS_full$hipp_mean <- winsorize_iqr(MRS_full$hipp_mean, 1.7)
MRS_full$cortical_thickness_adsignature_dickson <- winsorize_iqr(MRS_full$cortical_thickness_adsignature_dickson, 1.7)
MRS_full$m_m_acc <- winsorize_iqr(MRS_full$m_m_acc, 1.7)
MRS_full$m_m_precuneus <- winsorize_iqr(MRS_full$m_m_precuneus, 1.7)


#jmv::descriptives(data = MRS_full, vars = vars(hipp_mean_act),sd = TRUE,range = TRUE,skew = TRUE,kurt = TRUE)


#### Centrer et mettre au carré
# education
MRS_full$education_c  <- scale(MRS_full$education, center = TRUE, scale = FALSE)
MRS_full$education_sq <- MRS_full$education_c^2
# age_spectro
MRS_full$age_spectro_c  <- scale(MRS_full$age_spectro, center = TRUE, scale = FALSE)
MRS_full$age_spectro_sq <- MRS_full$age_spectro_c^2
# m_m_precuneus
MRS_full$m_m_precuneus_c  <- scale(MRS_full$m_m_precuneus, center = TRUE, scale = FALSE)
MRS_full$m_m_precuneus_sq <- MRS_full$m_m_precuneus_c^2
# m_m_acc
MRS_full$m_m_acc_c  <- scale(MRS_full$m_m_acc, center = TRUE, scale = FALSE)
MRS_full$m_m_acc_sq <- MRS_full$m_m_acc_c^2
# moca_corr_spectro
MRS_full$moca_corr_spectro_c  <- scale(MRS_full$moca_corr_spectro, center = TRUE, scale = FALSE)
MRS_full$moca_corr_spectro_sq <- MRS_full$moca_corr_spectro_c^2
# memoria_libre_correcte
MRS_full$memoria_libre_correcte_c  <- scale(MRS_full$memoria_libre_correcte, center = TRUE, scale = FALSE)
MRS_full$memoria_libre_correcte_sq <- MRS_full$memoria_libre_correcte_c^2
# face_name_rappel_differe_spectro
MRS_full$face_name_rappel_differe_spectro_c  <- scale(MRS_full$face_name_rappel_differe_spectro, center = TRUE, scale = FALSE)
MRS_full$face_name_rappel_differe_spectro_sq <- MRS_full$face_name_rappel_differe_spectro_c^2
# hip_l_nor_icv
MRS_full$hip_l_nor_icv_c  <- scale(MRS_full$hip_l_nor_icv, center = TRUE, scale = FALSE)
MRS_full$hip_l_nor_icv_sq <- MRS_full$hip_l_nor_icv_c^2
# hipp_mean
MRS_full$hipp_mean_c  <- scale(MRS_full$hipp_mean, center = TRUE, scale = FALSE)
MRS_full$hipp_mean_sq <- MRS_full$hipp_mean_c^2
# cortical_thickness_adsignature_dickson
MRS_full$cortical_thickness_adsignature_dickson_c  <- scale(MRS_full$cortical_thickness_adsignature_dickson, center = TRUE, scale = FALSE)
MRS_full$cortical_thickness_adsignature_dickson_sq <- MRS_full$cortical_thickness_adsignature_dickson_c^2
#  dprime_hit_fa
MRS_full$dprime_hit_fa_c  <- scale(MRS_full$dprime_hit_fa, center = TRUE, scale = FALSE)
MRS_full$dprime_hit_fa_sq <- MRS_full$dprime_hit_fa_c^2
#  associative_memory_performance
MRS_full$associative_memory_performance_c  <- scale(MRS_full$associative_memory_performance, center = TRUE, scale = FALSE)
MRS_full$associative_memory_performance_sq <- MRS_full$associative_memory_performance_c^2
#  associative_memory_performance_sylvie
MRS_full$associative_memory_performance_sylvie_c  <- scale(MRS_full$associative_memory_performance_sylvie, center = TRUE, scale = FALSE)
MRS_full$associative_memory_performance_sylvie_sq <- MRS_full$associative_memory_performance_sylvie_c^2
#  activation_hippocampus_l
MRS_full$activation_hippocampus_l_c  <- scale(MRS_full$activation_hippocampus_l, center = TRUE, scale = FALSE)
MRS_full$activation_hippocampus_l_sq <- MRS_full$activation_hippocampus_l_c^2
#  activation_hipp mean
MRS_full$hipp_mean_act_c  <- scale(MRS_full$hipp_mean_act, center = TRUE, scale = FALSE)
MRS_full$hipp_mean_act_sq <- MRS_full$hipp_mean_act_c^2
#  activation_parietal_sup_l 
MRS_full$activation_parietal_sup_l_c  <- scale(MRS_full$activation_parietal_sup_l, center = TRUE, scale = FALSE)
MRS_full$activation_parietal_sup_l_sq <- MRS_full$activation_parietal_sup_l_c^2
#  activation_temporal_inf_r
MRS_full$activation_temporal_inf_r_c  <- scale(MRS_full$activation_temporal_inf_r, center = TRUE, scale = FALSE)
MRS_full$activation_temporal_inf_r_sq <- MRS_full$activation_temporal_inf_r_c^2



##### Creer de banques de donees pour les analyses

#### Structure ####
MRS_S_Prec <- MRS_full[, c("m_m_precuneus", "m_m_precuneus_c", "m_m_precuneus_sq","hipp_mean","hipp_mean_c","hipp_mean_sq","sex", "age_spectro_c",
                           "hip_l_nor_icv", "hip_l_nor_icv_c", "hip_l_nor_icv_sq","cortical_thickness_adsignature_dickson",
                           "cortical_thickness_adsignature_dickson_c","cortical_thickness_adsignature_dickson_sq")] |> na.omit()

MRS_S_ACC <- MRS_full[, c("m_m_acc", "m_m_acc_c", "m_m_acc_sq","m_m_precuneus_sq","hipp_mean","hipp_mean_c","hipp_mean_sq","sex", "age_spectro_c",
                          "hip_l_nor_icv", "hip_l_nor_icv_c", "hip_l_nor_icv_sq","cortical_thickness_adsignature_dickson",
                          "cortical_thickness_adsignature_dickson_c", "cortical_thickness_adsignature_dickson_sq")] |> na.omit()



# Activaiton 
MRS_A_Prec <- MRS_full[, c(
  "m_m_precuneus", "m_m_precuneus_c", "m_m_precuneus_sq","activation_hippocampus_l", "activation_hippocampus_l_c","hipp_mean_act","hipp_mean_act_c","hipp_mean_act_sq",
  "activation_hippocampus_l_sq","activation_parietal_sup_l", "activation_parietal_sup_l_c", "sex", "age_spectro_c",
  "activation_parietal_sup_l_sq")] |> na.omit()

MRS_A_ACC <- MRS_full[, c("m_m_acc", "m_m_acc_c", "m_m_acc_sq","hipp_mean_act","hipp_mean_act_c","hipp_mean_act_sq",
                          "activation_hippocampus_l", "activation_hippocampus_l_c", "activation_hippocampus_l_sq","sex", "age_spectro_c",
                          "activation_parietal_sup_l", "activation_parietal_sup_l_c", "activation_parietal_sup_l_sq")] |> na.omit()

# Memory 
MRS_M_Prec <- MRS_full[, c("m_m_precuneus", "m_m_precuneus_c", "m_m_precuneus_sq","face_name_rappel_differe_spectro","face_name_rappel_differe_spectro_c","face_name_rappel_differe_spectro_sq",
                           "memoria_libre_correcte", "memoria_libre_correcte_c", "memoria_libre_correcte_sq","sex", "age_spectro_c",
                           "associative_memory_performance","associative_memory_performance_c",
                           "associative_memory_performance_sq")] |> na.omit()

MRS_M_ACC <- MRS_full[, c("m_m_acc", "m_m_acc_c", "m_m_acc_sq","memoria_libre_correcte","face_name_rappel_differe_spectro_c","face_name_rappel_differe_spectro_sq",
                          "memoria_libre_correcte_c", "memoria_libre_correcte_sq","face_name_rappel_differe_spectro","sex", "age_spectro_c",
                          "associative_memory_performance","associative_memory_performance_c",
                          "associative_memory_performance_sq")] |> na.omit()
############################## ANALYSES PRINCIPALES ###########################


################## Participants charactertistics  ################
# Sex
# chi-square
participant_sex <- table(MRS_full$diagnostic_nick, MRS_full$sex)
participant_sex
chisq.test(participant_sex)
fisher.test(participant_sex[c("HC","SCD+"), ])
fisher.test(participant_sex[c("HC","MCI"), ])
fisher.test(participant_sex[c("SCD+","MCI"), ]) 

# Age
jmv::descriptives(data = MRS_full,vars = vars(age_spectro), splitBy = "diagnostic_nick",sd = TRUE,hist = TRUE)
anova_age <- aov(age_spectro ~ diagnostic_nick, data = MRS_full)
summary(anova_age)               
TukeyHSD(anova_age) 


# Education
jmv::descriptives(data = MRS_full,vars = vars(education), splitBy = "diagnostic_nick",sd = TRUE, hist = TRUE)
anova_edu <- aov(education ~ diagnostic_nick, data = MRS_full)
summary(anova_edu)
TukeyHSD(anova_edu)

# MoCA
jmv::descriptives(data = MRS_full,vars = vars(moca_corr_spectro), splitBy = "diagnostic_nick",sd = TRUE, hist = TRUE)
anova_moca <- aov(moca_corr_spectro ~ diagnostic_nick, data = MRS_full)
summary(anova_moca)
TukeyHSD(anova_moca)


# Memoria Free Word recall
jmv::descriptives(data = MRS_full,vars = vars(memoria_libre_correcte), splitBy = "diagnostic_nick",sd = TRUE, hist = TRUE)
anova_free <- aov(memoria_libre_correcte ~ diagnostic_nick, data = MRS_full)
summary(anova_free)
TukeyHSD(anova_free)

# Face-name recall
jmv::descriptives(data = MRS_full,vars = vars(face_name_rappel_differe_spectro), splitBy = "diagnostic_nick",sd = TRUE, hist = TRUE)
anova_face <- aov(face_name_rappel_differe_spectro ~ diagnostic_nick, data = MRS_full)
summary(anova_face)
TukeyHSD(anova_face)


# Associative memory
jmv::descriptives(data = MRS_full,vars = vars(associative_memory_performance), splitBy = "diagnostic_nick",sd = TRUE, hist = TRUE)
anova_assoc <- aov(associative_memory_performance ~ diagnostic_nick, data = MRS_full)
summary(anova_assoc)
TukeyHSD(anova_assoc)


################ ANOVA glut  ###################
jmv::descriptives(data = MRS_full,vars = vars(m_m_precuneus, m_m_acc), splitBy = "diagnostic_nick",sd = TRUE, hist = TRUE, iqr = TRUE)

#  ACC
anova_acc <- aov(m_m_acc ~ diagnostic_nick, data = MRS_full)
summary(anova_acc)
TukeyHSD(anova_acc)
eta_squared(anova_acc)          
# For Precuneus
anova_precuneus <- aov(m_m_precuneus ~ diagnostic_nick, data = MRS_full)
summary(anova_precuneus)
TukeyHSD(anova_precuneus)
eta_squared(anova_precuneus)          




######################### Polynomial analyses  ####################################

#### structure and memory ######


# Model 1 m_m_ACC ~ hipp mean
model_lin_acc_hipp_mean <- lm(m_m_acc ~ hipp_mean_c, data = MRS_S_ACC)
summary(model_lin_acc_hipp_mean)
AIC(model_lin_acc_hipp_mean)
model_quad_acc_hipp_mean <- lm(m_m_acc ~ hipp_mean_c + hipp_mean_sq, data = MRS_S_ACC)
summary(model_quad_acc_hipp_mean)
AIC(model_quad_acc_hipp_mean)


# Model 2 m_m_Precuneus ~ hipp mean 
model_lin_prec_hipp_mean <- lm(m_m_precuneus ~ hipp_mean_c, data = MRS_S_Prec)
summary(model_lin_prec_hipp_mean)
AIC(model_lin_prec_hipp_mean)
model_quad_prec_hipp_mean <- lm(m_m_precuneus ~ hipp_mean_c + hipp_mean_sq, data = MRS_S_Prec)
summary(model_quad_prec_hipp_mean)
AIC(model_quad_prec_hipp_mean)


# Model 3 m_m_ACC ~ thickness
model_lin_acc_thick <- lm(m_m_acc ~ cortical_thickness_adsignature_dickson_c, data = MRS_S_ACC)
summary(model_lin_acc_thick)
AIC(model_lin_acc_thick)
model_quad_acc_thick <- lm(m_m_acc ~ cortical_thickness_adsignature_dickson_c + cortical_thickness_adsignature_dickson_sq, data = MRS_S_ACC)
summary(model_quad_acc_thick)
AIC(model_quad_acc_thick)


# Model 4 m_m_Precuneus ~ thickness
model_lin_prec_thick <- lm(m_m_precuneus ~ cortical_thickness_adsignature_dickson_c, data = MRS_S_Prec)
summary(model_lin_prec_thick)
AIC(model_lin_prec_thick)
model_quad_prec_thick <- lm(m_m_precuneus ~ cortical_thickness_adsignature_dickson_c + cortical_thickness_adsignature_dickson_sq, data = MRS_S_Prec)
summary(model_quad_prec_thick)
AIC(model_quad_prec_thick)


#### Activaiton ######
names(MRS_A_ACC)
names(MRS_A_Prec)
## Model 5 ACC ~ activaiton parietal
model_lin_acc_sup_act_rev <- lm(m_m_acc ~ activation_parietal_sup_l, data = MRS_A_ACC)
summary(model_lin_acc_sup_act_rev)
AIC(model_lin_acc_sup_act_rev)

model_quad_acc_sup_act_rev <- lm(m_m_acc ~ activation_parietal_sup_l + activation_parietal_sup_l_sq, data = MRS_A_ACC)
summary(model_quad_acc_sup_act_rev)
AIC(model_quad_acc_sup_act_rev)

## Model 6 recuneus ~ activaiton parietal
model_lin_prec_sup_act_rev <- lm(m_m_precuneus ~ activation_parietal_sup_l, data = MRS_A_Prec)
summary(model_lin_prec_sup_act_rev)
AIC(model_lin_prec_sup_act_rev)

model_quad_prec_sup_act_rev <- lm(m_m_precuneus ~ activation_parietal_sup_l + activation_parietal_sup_l_sq, data = MRS_A_Prec)
summary(model_quad_prec_sup_act_rev)
AIC(model_quad_prec_sup_act_rev)


# Model 7 acc ~ hipp activation
model_lin_acc_hipp_act <- lm(m_m_acc ~ activation_hippocampus_l_c, data = MRS_A_ACC)
summary(model_lin_acc_hipp_act)
AIC(model_lin_acc_hipp_act)
model_quad_acc_hipp_act <- lm(m_m_acc ~ activation_hippocampus_l_c + activation_hippocampus_l_sq, data = MRS_A_ACC)
summary(model_quad_acc_hipp_act)
AIC(model_quad_acc_hipp_act)

# Model 8 precuneus  ~ hipp activation
model_lin_prec_hipp_act <- lm(m_m_precuneus ~ activation_hippocampus_l_c, data = MRS_A_Prec)
summary(model_lin_prec_hipp_act)
AIC(model_lin_prec_hipp_act)
model_quad_prec_hipp_act <- lm(m_m_precuneus ~ activation_hippocampus_l_c + activation_hippocampus_l_sq, data = MRS_A_Prec)
summary(model_quad_prec_hipp_act)
AIC(model_quad_prec_hipp_act)

# Model 7 acc ~ hipp activation mean
model_lin_acc_hipp_act <- lm(m_m_acc ~ hipp_mean_act_c, data = MRS_A_ACC)
summary(model_lin_acc_hipp_act)
AIC(model_lin_acc_hipp_act)
model_quad_acc_hipp_act <- lm(m_m_acc ~ hipp_mean_act_c + hipp_mean_act_sq, data = MRS_A_ACC)
summary(model_quad_acc_hipp_act)
AIC(model_quad_acc_hipp_act)

# Model 8 precuneus  ~ hipp activation mean
model_lin_prec_hipp_act <- lm(m_m_precuneus ~ hipp_mean_act_c, data = MRS_A_Prec)
summary(model_lin_prec_hipp_act)
AIC(model_lin_prec_hipp_act)
model_quad_prec_hipp_act <- lm(m_m_precuneus ~ hipp_mean_act_c + hipp_mean_act_sq, data = MRS_A_Prec)
summary(model_quad_prec_hipp_act)
AIC(model_quad_prec_hipp_act)

### Memory ###

# Model 9 memoria ~ ACC Glu
model_lin_acc_memor <- lm(memoria_libre_correcte ~ m_m_acc_c, data = MRS_M_ACC)
summary(model_lin_acc_memor)
AIC(model_lin_acc_memor)
model_quad_acc_memor <- lm(memoria_libre_correcte ~ m_m_acc_c + m_m_acc_sq, data = MRS_M_ACC)
summary(model_quad_acc_memor)
AIC(model_quad_acc_memor)

# Model 10 memoria ~ Precuneus Glu
model_lin_prec_memor <- lm(memoria_libre_correcte ~ m_m_precuneus_c, data = MRS_M_Prec)
summary(model_lin_prec_memor)
AIC(model_lin_prec_memor)
model_quad_prec_memor <- lm(memoria_libre_correcte ~ m_m_precuneus_c + m_m_precuneus_sq, data = MRS_M_Prec)
summary(model_quad_prec_memor)
AIC(model_quad_prec_memor)

# Model 11 face name ~  ACC glu
model_lin_acc_memor <- lm(face_name_rappel_differe_spectro ~ m_m_acc_c, data = MRS_M_ACC)
summary(model_lin_acc_memor)
AIC(model_lin_acc_memor)
model_quad_acc_memor <- lm(face_name_rappel_differe_spectro ~ m_m_acc_c + m_m_acc_sq, data = MRS_M_ACC)
summary(model_quad_acc_memor)
AIC(model_quad_acc_memor)

# Model 12 face name ~ Precuneus glu
model_lin_prec_memor <- lm(face_name_rappel_differe_spectro ~ m_m_precuneus_c, data = MRS_M_Prec)
summary(model_lin_prec_memor)
AIC(model_lin_prec_memor)
model_quad_prec_memor <- lm(face_name_rappel_differe_spectro ~ m_m_precuneus_c + m_m_precuneus_sq, data = MRS_M_Prec)
summary(model_quad_prec_memor)
AIC(model_quad_prec_memor)





############## Mediation analysis ##################
library(mediation)
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




