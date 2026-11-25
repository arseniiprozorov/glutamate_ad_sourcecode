library(readxl)   
library(janitor)  
library(jmv)
library(effectsize)
library(interactions)
library(emmeans)
library(mediation)

## Glutamate Relates to Structural and Functional markers of Disease Severity in early Alzheimer’s disease  ######
## Arsenii Prozorov 
#ANALYSES PRÉLIMINAIRES :
#Création d’une banque de données
import_data_arsenii_20250609_hippR <- read_excel("C:/Users/okkam/Desktop/labo/article 1/import_data_arsenii_20261017.xlsx")
MRS_full <- import_data_arsenii_20250609_hippR

# Clean the column name
MRS_full <- janitor::clean_names(MRS_full)
names(MRS_full)
# Convertir  en numérique
MRS_full$education <- as.numeric(MRS_full$education)
MRS_full$moca_corr_spectro <- as.numeric(MRS_full$moca_corr_spectro)
MRS_full$associative_memory_performance <- as.numeric(MRS_full$associative_memory_performance)
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
MRS_full$activation_parietal_sup_l <- winsorize_iqr(MRS_full$activation_parietal_sup_l, 1.5)
MRS_full$hipp_mean_act <- winsorize_iqr(MRS_full$hipp_mean_act, 1.5)
MRS_full$moca_corr_spectro <- winsorize_iqr(MRS_full$moca_corr_spectro, 1.5)
MRS_full$associative_memory_performance <- winsorize_iqr(MRS_full$associative_memory_performance, 1.5)

MRS_full$hipp_mean <- winsorize_iqr(MRS_full$hipp_mean, 1.5)
MRS_full$cortical_thickness_adsignature_dickson <- winsorize_iqr(MRS_full$cortical_thickness_adsignature_dickson, 1.5)
MRS_full$m_m_acc <- winsorize_iqr(MRS_full$m_m_acc, 1.5)
MRS_full$m_m_precuneus <- winsorize_iqr(MRS_full$m_m_precuneus, 1.5)


#jmv::descriptives(data = MRS_full, vars = vars(hipp_mean_act, activation_parietal_sup_l),sd = TRUE,range = TRUE,skew = TRUE,kurt = TRUE)


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
# memoria_libre_correcte
MRS_full$memoria_libre_correcte_c  <- scale(MRS_full$memoria_libre_correcte, center = TRUE, scale = FALSE)
MRS_full$memoria_libre_correcte_sq <- MRS_full$memoria_libre_correcte_c^2
# face_name_rappel_differe_spectro
MRS_full$face_name_rappel_differe_spectro_c  <- scale(MRS_full$face_name_rappel_differe_spectro, center = TRUE, scale = FALSE)
MRS_full$face_name_rappel_differe_spectro_sq <- MRS_full$face_name_rappel_differe_spectro_c^2
# hipp_mean
MRS_full$hipp_mean_c  <- scale(MRS_full$hipp_mean, center = TRUE, scale = FALSE)
MRS_full$hipp_mean_sq <- MRS_full$hipp_mean_c^2
# cortical_thickness_adsignature_dickson
MRS_full$cortical_thickness_adsignature_dickson_c  <- scale(MRS_full$cortical_thickness_adsignature_dickson, center = TRUE, scale = FALSE)
MRS_full$cortical_thickness_adsignature_dickson_sq <- MRS_full$cortical_thickness_adsignature_dickson_c^2
#  activation_hipp mean
MRS_full$hipp_mean_act_c  <- scale(MRS_full$hipp_mean_act, center = TRUE, scale = FALSE)
MRS_full$hipp_mean_act_sq <- MRS_full$hipp_mean_act_c^2
#  activation_parietal_sup_l 
MRS_full$activation_parietal_sup_l_c  <- scale(MRS_full$activation_parietal_sup_l, center = TRUE, scale = FALSE)
MRS_full$activation_parietal_sup_l_sq <- MRS_full$activation_parietal_sup_l_c^2


############################## ANALYSES PRINCIPALES ###########################


################## Participants charactertistics  ################
#sink("Table 1. participant_characteristics_results.txt") 
# Sex
# chi-square
participant_sex <- table(MRS_full$diagnostic_nick, MRS_full$sex)
participant_sex
chisq.test(participant_sex)
fisher.test(participant_sex[c("HC","SCD+"), ])
fisher.test(participant_sex[c("HC","MCI"), ])
fisher.test(participant_sex[c("SCD+","MCI"), ]) 

# Age
jmv::descriptives(data = MRS_full,vars = vars(age_spectro), splitBy = "diagnostic_nick",sd = TRUE,hist = TRUE, range = TRUE)
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


# Face-name recall
jmv::descriptives(data = MRS_full,vars = vars(face_name_rappel_differe_spectro), splitBy = "diagnostic_nick",sd = TRUE, hist = TRUE)
anova_face <- aov(face_name_rappel_differe_spectro ~ diagnostic_nick, data = MRS_full)
summary(anova_face)
TukeyHSD(anova_face)

# Associative memory (task)
jmv::descriptives(data = MRS_full,vars = vars(associative_memory_performance), splitBy = "diagnostic_nick",sd = TRUE, hist = TRUE)
anova_asso <- aov(associative_memory_performance ~ diagnostic_nick, data = MRS_full)
summary(anova_asso)
TukeyHSD(anova_asso)

###### ANOVA glut  #
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


#### ANOVA Structure ####
jmv::descriptives(data = MRS_full,vars = vars(hipp_mean, cortical_thickness_adsignature_dickson), splitBy = "diagnostic_nick",sd = TRUE, hist = TRUE, iqr = TRUE)

#  Hipp volume
anova_hipp_vol <- aov(hipp_mean ~ diagnostic_nick, data = MRS_full)
summary(anova_hipp_vol)
TukeyHSD(anova_hipp_vol)
eta_squared(anova_hipp_vol)          
# Thickness
anova_thick <- aov(cortical_thickness_adsignature_dickson ~ diagnostic_nick, data = MRS_full)
summary(anova_thick)


## ANOVA Activation  #####
jmv::descriptives(data = MRS_full,vars = vars(hipp_mean_act, activation_parietal_sup_l), splitBy = "diagnostic_nick",sd = TRUE, hist = TRUE, iqr = TRUE)

#  Hipp activation
anova_hipp_act <- aov(hipp_mean_act ~ diagnostic_nick, data = MRS_full)
summary(anova_hipp_act)
        
# Superior parietal
anova_parietal_act <- aov(activation_parietal_sup_l ~ diagnostic_nick, data = MRS_full)
summary(anova_parietal_act)

#sink()


######################### Polynomial analyses  ####################################
#stats::lm
#stats::lm.fit
#stats:::anova.lmlist
#stats:::anova.lmlist
#### structure and memory ######
#sink("Table 2,3 linear and quadratic models between proxies of disease severity .txt") 

# Model 1 m_m_ACC ~ hipp mean
model_lin_acc_hipp_mean <- lm(m_m_acc ~ hipp_mean_c, data = MRS_full)
summary(model_lin_acc_hipp_mean)
model_quad_acc_hipp_mean <- lm(m_m_acc ~ hipp_mean_c + hipp_mean_sq, data = MRS_full)
summary(model_quad_acc_hipp_mean)

AIC(model_lin_acc_hipp_mean, model_quad_acc_hipp_mean)
anova(model_lin_acc_hipp_mean, model_quad_acc_hipp_mean)

# Model 2 m_m_Precuneus ~ hipp mean 
model_lin_prec_hipp_mean <- lm(m_m_precuneus ~ hipp_mean_c, data = MRS_full)
summary(model_lin_prec_hipp_mean)
model_quad_prec_hipp_mean <- lm(m_m_precuneus ~ hipp_mean_c + hipp_mean_sq, data = MRS_full)
summary(model_quad_prec_hipp_mean)

AIC(model_lin_prec_hipp_mean, model_quad_prec_hipp_mean)
anova(model_lin_prec_hipp_mean, model_quad_prec_hipp_mean)

# Model 3 m_m_ACC ~ thickness
model_lin_acc_thick <- lm(m_m_acc ~ cortical_thickness_adsignature_dickson_c, data = MRS_full)
summary(model_lin_acc_thick)
model_quad_acc_thick <- lm(m_m_acc ~ cortical_thickness_adsignature_dickson_c + cortical_thickness_adsignature_dickson_sq, data = MRS_full)
summary(model_quad_acc_thick)

AIC(model_lin_acc_thick)
AIC(model_quad_acc_thick)
anova(model_lin_acc_thick, model_quad_acc_thick)

# Model 4 m_m_Precuneus ~ thickness
model_lin_prec_thick <- lm(m_m_precuneus ~ cortical_thickness_adsignature_dickson_c , data = MRS_full)
summary(model_lin_prec_thick)
model_quad_prec_thick <- lm(m_m_precuneus ~ cortical_thickness_adsignature_dickson_c + cortical_thickness_adsignature_dickson_sq, data = MRS_full)
summary(model_quad_prec_thick)

AIC(model_lin_prec_thick, model_quad_prec_thick)
anova(model_lin_prec_thick, model_quad_prec_thick)

#### Activaiton ######

## Model 5 ACC ~ activaiton parietal
model_lin_acc_sup_act_rev <- lm(m_m_acc ~ activation_parietal_sup_l_c, data = MRS_full)
summary(model_lin_acc_sup_act_rev)
model_quad_acc_sup_act_rev <- lm(m_m_acc ~ activation_parietal_sup_l_c + activation_parietal_sup_l_sq, data = MRS_full)
summary(model_quad_acc_sup_act_rev)

AIC(model_lin_acc_sup_act_rev, model_quad_acc_sup_act_rev)
anova(model_lin_acc_sup_act_rev, model_quad_acc_sup_act_rev)

## Model 6 recuneus ~ activaiton parietal
model_lin_prec_sup_act_rev <- lm(m_m_precuneus ~ activation_parietal_sup_l_c, data = MRS_full)
summary(model_lin_prec_sup_act_rev)
model_quad_prec_sup_act_rev <- lm(m_m_precuneus ~ activation_parietal_sup_l_c + activation_parietal_sup_l_sq, data = MRS_full)
summary(model_quad_prec_sup_act_rev)

AIC(model_lin_prec_sup_act_rev, model_quad_prec_sup_act_rev)
anova(model_lin_prec_sup_act_rev, model_quad_prec_sup_act_rev)


# Model 7 acc ~ hipp activation mean
model_lin_acc_hipp_act <- lm(m_m_acc ~ hipp_mean_act_c , data = MRS_full)
summary(model_lin_acc_hipp_act)
model_quad_acc_hipp_act <- lm(m_m_acc ~ hipp_mean_act_c + hipp_mean_act_sq, data = MRS_full)
summary(model_quad_acc_hipp_act)

AIC(model_lin_acc_hipp_act, model_quad_acc_hipp_act)
anova(model_lin_acc_hipp_act, model_quad_acc_hipp_act)

# Model 8 precuneus  ~ hipp activation mean
model_lin_prec_hipp_act <- lm(m_m_precuneus ~ hipp_mean_act_c, data = MRS_full)
summary(model_lin_prec_hipp_act)
model_quad_prec_hipp_act <- lm(m_m_precuneus ~ hipp_mean_act_c + hipp_mean_act_sq, data = MRS_full)
summary(model_quad_prec_hipp_act)

AIC(model_lin_prec_hipp_act)
AIC(model_quad_prec_hipp_act)
anova(model_lin_prec_hipp_act, model_quad_prec_hipp_act)


##### Memory #####

# Model 9 memoria ~ ACC Glu
model_lin_acc_memor <- lm(memoria_libre_correcte ~ m_m_acc_c, data = MRS_full)
summary(model_lin_acc_memor)
model_quad_acc_memor <- lm(memoria_libre_correcte ~ m_m_acc_c + m_m_acc_sq, data = MRS_full)
summary(model_quad_acc_memor)

AIC(model_lin_acc_memor)
AIC(model_quad_acc_memor)
anova(model_lin_acc_memor, model_quad_acc_memor)

# Model 10 memoria ~ Precuneus Glu
model_lin_prec_memor <- lm(memoria_libre_correcte ~ m_m_precuneus_c , data = MRS_full)
summary(model_lin_prec_memor)
model_quad_prec_memor <- lm(memoria_libre_correcte ~ m_m_precuneus_c + m_m_precuneus_sq, data = MRS_full)
summary(model_quad_prec_memor)

AIC(model_lin_prec_memor)
AIC(model_quad_prec_memor)
anova(model_lin_prec_memor, model_quad_prec_memor)



# Model 11 face name ~  ACC glu
model_lin_acc_facename <- lm(face_name_rappel_differe_spectro ~ m_m_acc_c, data = MRS_full)
summary(model_lin_acc_facename)
model_quad_acc_facename <- lm(face_name_rappel_differe_spectro ~ m_m_acc_c + m_m_acc_sq, data = MRS_full)
summary(model_quad_acc_facename)

AIC(model_lin_acc_facename)
AIC(model_quad_acc_facename)
anova(model_lin_acc_facename, model_quad_acc_facename)

# Model 12 face name ~ Precuneus glu
model_lin_prec_facename <- lm(face_name_rappel_differe_spectro ~ m_m_precuneus_c, data = MRS_full)
summary(model_lin_prec_facename)
model_quad_prec_facename <- lm(face_name_rappel_differe_spectro ~ m_m_precuneus_c + m_m_precuneus_sq, data = MRS_full)
summary(model_quad_prec_facename)

AIC(model_lin_prec_facename)
AIC(model_quad_prec_facename)
anova(model_lin_prec_facename, model_quad_prec_facename)


#sink() 



names(MRS_full)
############## Mediation analysis ##################
set.seed(1)
# Select raw columns
mediation_combined <- MRS_full[, c("m_m_precuneus", "m_m_acc", 
                                   "hipp_mean", "cortical_thickness_adsignature_dickson", 
                                   "hipp_mean_act", "activation_parietal_sup_l", "education", "diagnostic_nick",
                                   "sex", "age_spectro_c")] |> na.omit()


# Z-score individual variables
mediation_combined$prec_glu_z     <- as.numeric(scale(mediation_combined$m_m_precuneus))
mediation_combined$acc_glu_z      <- as.numeric(scale(mediation_combined$m_m_acc))
mediation_combined$hipp_vol_z     <- as.numeric(scale(mediation_combined$hipp_mean))
mediation_combined$thick_z        <- as.numeric(scale(mediation_combined$cortical_thickness_adsignature_dickson))
mediation_combined$hipp_act_z     <- as.numeric(scale(mediation_combined$hipp_mean_act))
mediation_combined$parietal_act_z <- as.numeric(scale(mediation_combined$activation_parietal_sup_l))

# Create composite "mean" variables
mediation_combined$mean_structure  <- rowMeans(mediation_combined[, c("hipp_vol_z", "thick_z")])
mediation_combined$mean_activation <- rowMeans(mediation_combined[, c("hipp_act_z", "parietal_act_z")])
mediation_combined$mean_glu        <- rowMeans(mediation_combined[, c("prec_glu_z", "acc_glu_z")])




###############################  GLOBAL COMBINED MODEL ###########################
#sink("mediation_main.txt") 
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


#sink() 





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
zmediate_results <- ZYmediate(x, y, med, nboot = 1000, alpha = 0.05, kappa = 0.05)
zmediate_results

