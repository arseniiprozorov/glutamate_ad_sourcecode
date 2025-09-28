library(readxl)   
library(janitor)  
library(jmv)
library(effectsize)
      

## Glutamate Relates to Structural and Functional markers of Disease Severity in early Alzheimer’s disease  ######
## Arsenii Prozorov 


#ANALYSES PRÉLIMINAIRES :
#Création d’une banque de données

import_data_arsenii_20250609 <- read_excel("C:/Users/okkam/Desktop/labo/article 1/supplementary/import_data_arsenii_20250609.xlsx")
MRS_full <- import_data_arsenii_20250609

# Clean the column names
MRS_full <- janitor::clean_names(MRS_full)

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


#Vérifier les varibales
lapply(MRS_full,class)   
levels(MRS_full$diagnostic_nick)
table(MRS_full$diagnostic_nick)


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
MRS_full$moca_corr_spectro <- winsorize_iqr(MRS_full$moca_corr_spectro, 1.5)
MRS_full$dprime_hit_fa <- winsorize_iqr(MRS_full$dprime_hit_fa, 1.5)
MRS_full$associative_memory_performance_sylvie <- winsorize_iqr(MRS_full$associative_memory_performance_sylvie, 1.5)


# Apply winsorization for 1.7 IQR variables
MRS_full$hip_l_nor_icv <- winsorize_iqr(MRS_full$hip_l_nor_icv, 1.7)
MRS_full$cortical_thickness_adsignature_dickson <- winsorize_iqr(MRS_full$cortical_thickness_adsignature_dickson, 1.7)
MRS_full$m_m_acc <- winsorize_iqr(MRS_full$m_m_acc, 1.7)
MRS_full$m_m_precuneus <- winsorize_iqr(MRS_full$m_m_precuneus, 1.7)


#jmv::descriptives(data = MRS_full, vars = vars(memoria_libre_correcte),sd = TRUE,range = TRUE,skew = TRUE,kurt = TRUE)


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
#  activation_parietal_sup_l 
MRS_full$activation_parietal_sup_l_c  <- scale(MRS_full$activation_parietal_sup_l, center = TRUE, scale = FALSE)
MRS_full$activation_parietal_sup_l_sq <- MRS_full$activation_parietal_sup_l_c^2
#  activation_temporal_inf_r
MRS_full$activation_temporal_inf_r_c  <- scale(MRS_full$activation_temporal_inf_r, center = TRUE, scale = FALSE)
MRS_full$activation_temporal_inf_r_sq <- MRS_full$activation_temporal_inf_r_c^2



##### Creer de banques de donees pour les analyses

#### Structure and Memory ####
MRS_S_Prec <- MRS_full[, c("m_m_precuneus", "m_m_precuneus_c", "m_m_precuneus_sq",
                           "hip_l_nor_icv", "hip_l_nor_icv_c", "hip_l_nor_icv_sq","cortical_thickness_adsignature_dickson",
                           "cortical_thickness_adsignature_dickson_c","cortical_thickness_adsignature_dickson_sq")] |> na.omit()

MRS_S_ACC <- MRS_full[, c("m_m_acc", "m_m_acc_c", "m_m_acc_sq",
                          "hip_l_nor_icv", "hip_l_nor_icv_c", "hip_l_nor_icv_sq","cortical_thickness_adsignature_dickson",
                          "cortical_thickness_adsignature_dickson_c", "cortical_thickness_adsignature_dickson_sq")] |> na.omit()

MRS_M_Prec <- MRS_full[, c("m_m_precuneus", "m_m_precuneus_c", "m_m_precuneus_sq",
                           "memoria_libre_correcte", "memoria_libre_correcte_c", "memoria_libre_correcte_sq",
                           "associative_memory_performance","associative_memory_performance_c",
                           "associative_memory_performance_sq")] |> na.omit()

MRS_M_ACC <- MRS_full[, c("m_m_acc", "m_m_acc_c", "m_m_acc_sq","memoria_libre_correcte",
                          "memoria_libre_correcte_c", "memoria_libre_correcte_sq",
                          "associative_memory_performance","associative_memory_performance_c",
                          "associative_memory_performance_sq")] |> na.omit()

# Activaiton 
MRS_A_Prec <- MRS_full[, c(
  "m_m_precuneus", "m_m_precuneus_c", "m_m_precuneus_sq","activation_hippocampus_l", "activation_hippocampus_l_c", 
  "activation_hippocampus_l_sq","activation_parietal_sup_l", "activation_parietal_sup_l_c", 
  "activation_parietal_sup_l_sq")] |> na.omit()

MRS_A_ACC <- MRS_full[, c("m_m_acc", "m_m_acc_c", "m_m_acc_sq",
                          "activation_hippocampus_l", "activation_hippocampus_l_c", "activation_hippocampus_l_sq",
                          "activation_parietal_sup_l", "activation_parietal_sup_l_c", "activation_parietal_sup_l_sq")] |> na.omit()


# Moderation
#### Structure and Memory ####
MRS_SM_Prec <- MRS_full[, c("m_m_precuneus", "m_m_precuneus_c", "m_m_precuneus_sq",
                            "hip_l_nor_icv", "hip_l_nor_icv_c", "hip_l_nor_icv_sq","cortical_thickness_adsignature_dickson",
                            "cortical_thickness_adsignature_dickson_c","cortical_thickness_adsignature_dickson_sq",
                            "memoria_libre_correcte", "memoria_libre_correcte_c", "memoria_libre_correcte_sq")] |> na.omit()

MRS_SM_ACC <- MRS_full[, c("m_m_acc", "m_m_acc_c", "m_m_acc_sq",
                           "hip_l_nor_icv", "hip_l_nor_icv_c", "hip_l_nor_icv_sq","cortical_thickness_adsignature_dickson",
                           "cortical_thickness_adsignature_dickson_c", "cortical_thickness_adsignature_dickson_sq","memoria_libre_correcte",
                           "memoria_libre_correcte_c", "memoria_libre_correcte_sq")] |> na.omit()


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


# d prime
jmv::descriptives(data = MRS_full,vars = vars(dprime_hit_fa), splitBy = "diagnostic_nick",sd = TRUE, hist = TRUE)
anova_dprime <- aov(dprime_hit_fa ~ diagnostic_nick, data = MRS_full)
summary(anova_dprime)
TukeyHSD(anova_dprime)

# Associative memory
jmv::descriptives(data = MRS_full,vars = vars(associative_memory_performance), splitBy = "diagnostic_nick",sd = TRUE, hist = TRUE)
anova_assoc <- aov(associative_memory_performance ~ diagnostic_nick, data = MRS_full)
summary(anova_assoc)
TukeyHSD(anova_assoc)




########################### ANOVA   ##########################

jmv::descriptives(data = MRS_full,vars = vars(m_m_precuneus, m_m_acc), splitBy = "diagnostic_nick",sd = TRUE, hist = TRUE)
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

## t test 
jmv::descriptives(data = MRS_full, vars = vars(m_m_precuneus,m_m_acc),splitBy = "sex", sd = TRUE, hist = TRUE)
t.test(m_m_precuneus ~ sex, data = MRS_full, var.equal = TRUE)



######################### Polynomial analyses  #####################################
names(MRS_full)
#### structure and memory ######

# Model 1 m_m_ACC ~ hipp left
model_lin_acc_hipp <- lm(m_m_acc ~ hip_l_nor_icv_c, data = MRS_S_ACC)
summary(model_lin_acc_hipp)
AIC(model_lin_acc_hipp)
model_quad_acc_hipp <- lm(m_m_acc ~ hip_l_nor_icv_c + hip_l_nor_icv_sq, data = MRS_S_ACC)
summary(model_quad_acc_hipp)
AIC(model_quad_acc_hipp)

# Model 2 m_m_Precuneus ~ hipp left 
model_lin_prec_hipp <- lm(m_m_precuneus ~ hip_l_nor_icv_c, data = MRS_S_Prec)
summary(model_lin_prec_hipp)
AIC(model_lin_prec_hipp)
model_quad_prec_hipp <- lm(m_m_precuneus ~ hip_l_nor_icv_c + hip_l_nor_icv_sq, data = MRS_S_Prec)
summary(model_quad_prec_hipp)
AIC(model_quad_prec_hipp)

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


### Memory ###

# Model 5 m_m_ACC ~ memoria
model_lin_acc_memor <- lm(m_m_acc ~ memoria_libre_correcte_c, data = MRS_M_ACC)
summary(model_lin_acc_memor)
AIC(model_lin_acc_memor)
model_quad_acc_memor <- lm(m_m_acc ~ memoria_libre_correcte_c + memoria_libre_correcte_sq, data = MRS_M_ACC)
summary(model_quad_acc_memor)
AIC(model_quad_acc_memor)

# Model 6 m_m_Precuneus ~ memoria
model_lin_prec_memor <- lm(m_m_precuneus ~ memoria_libre_correcte_c, data = MRS_M_Prec)
summary(model_lin_prec_memor)
AIC(model_lin_prec_memor)
model_quad_prec_memor <- lm(m_m_precuneus ~ memoria_libre_correcte_c + memoria_libre_correcte_sq, data = MRS_M_Prec)
summary(model_quad_prec_memor)
AIC(model_quad_prec_memor)


#### Activaiton ######

## Model 7ACC ~ activaiton parietal
model_lin_acc_sup_act_rev <- lm(m_m_acc_c ~ activation_parietal_sup_l, data = MRS_A_ACC)
summary(model_lin_acc_sup_act_rev)
AIC(model_lin_acc_sup_act_rev)

model_quad_acc_sup_act_rev <- lm(m_m_acc_c ~ activation_parietal_sup_l + activation_parietal_sup_l_sq, data = MRS_A_ACC)
summary(model_quad_acc_sup_act_rev)
AIC(model_quad_acc_sup_act_rev)

## Model 8 recuneus ~ activaiton parietal
model_lin_prec_sup_act_rev <- lm(m_m_precuneus_c ~ activation_parietal_sup_l, data = MRS_A_Prec)
summary(model_lin_prec_sup_act_rev)
AIC(model_lin_prec_sup_act_rev)

model_quad_prec_sup_act_rev <- lm(m_m_precuneus_c ~ activation_parietal_sup_l + activation_parietal_sup_l_sq, data = MRS_A_Prec)
summary(model_quad_prec_sup_act_rev)
AIC(model_quad_prec_sup_act_rev)


# Model 9acc ~ hipp activation
model_lin_acc_hipp_act <- lm(m_m_acc ~ activation_hippocampus_l_c, data = MRS_A_ACC)
summary(model_lin_acc_hipp_act)
AIC(model_lin_acc_hipp_act)
model_quad_acc_hipp_act <- lm(m_m_acc ~ activation_hippocampus_l_c + activation_hippocampus_l_sq, data = MRS_A_ACC)
summary(model_quad_acc_hipp_act)
AIC(model_quad_acc_hipp_act)

# Model 10precuneus  ~ hipp activation
model_lin_prec_hipp_act <- lm(m_m_precuneus ~ activation_hippocampus_l_c, data = MRS_A_Prec)
summary(model_lin_prec_hipp_act)
AIC(model_lin_prec_hipp_act)
model_quad_prec_hipp_act <- lm(m_m_precuneus ~ activation_hippocampus_l_c + activation_hippocampus_l_sq, data = MRS_A_Prec)
summary(model_quad_prec_hipp_act)
AIC(model_quad_prec_hipp_act)



### moderation ####
##  Hipp x glut
model_mod_hipp_acc  <- lm(memoria_libre_correcte ~ m_m_acc_c + hip_l_nor_icv_c:m_m_acc_c, data = MRS_SM_ACC)
summary(model_mod_hipp_acc)
model_mod_hipp_prec  <- lm(memoria_libre_correcte ~ m_m_precuneus_c + hip_l_nor_icv_c:m_m_precuneus_c, data = MRS_SM_Prec)
summary(model_mod_hipp_prec)


## thick x glut
model_mod_thick_acc <- lm(memoria_libre_correcte ~ m_m_acc_c + cortical_thickness_adsignature_dickson_c:m_m_acc_c, data = MRS_SM_ACC)
summary(model_mod_thick_acc)
model_mod_thick_prec <- lm(memoria_libre_correcte ~ m_m_precuneus_c + cortical_thickness_adsignature_dickson_c:m_m_precuneus_c, data = MRS_SM_Prec)
summary(model_mod_thick_prec)


