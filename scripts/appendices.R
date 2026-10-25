###   Controling for sex,  age,  education, diagnostic

#### structure  ######

# Model 1 m_m_ACC ~ hipp mean
model_lin_acc_hipp_mean <- lm(m_m_acc ~ hipp_mean_c + age_spectro_c + sex + education + diagnostic_nick, data = MRS_S_ACC)
summary(model_lin_acc_hipp_mean)
AIC(model_lin_acc_hipp_mean)
model_quad_acc_hipp_mean <- lm(m_m_acc ~ hipp_mean_c + hipp_mean_sq + age_spectro_c + sex + education + diagnostic_nick, data = MRS_S_ACC)
summary(model_quad_acc_hipp_mean)
AIC(model_quad_acc_hipp_mean)

anova(model_lin_acc_hipp_mean, model_quad_acc_hipp_mean)

# Model 2 m_m_Precuneus ~ hipp mean 
model_lin_prec_hipp_mean <- lm(m_m_precuneus ~ hipp_mean_c + age_spectro_c + sex + education + diagnostic_nick, data = MRS_S_Prec)
summary(model_lin_prec_hipp_mean)
AIC(model_lin_prec_hipp_mean)
model_quad_prec_hipp_mean <- lm(m_m_precuneus ~ hipp_mean_c + hipp_mean_sq + age_spectro_c + sex + education + diagnostic_nick, data = MRS_S_Prec)
summary(model_quad_prec_hipp_mean)
AIC(model_quad_prec_hipp_mean)

anova(model_lin_prec_hipp_mean, model_quad_prec_hipp_mean)

# Model 3 m_m_ACC ~ thickness
model_lin_acc_thick <- lm(m_m_acc ~ cortical_thickness_adsignature_dickson_c + age_spectro_c + sex + education + diagnostic_nick, data = MRS_S_ACC)
summary(model_lin_acc_thick)
AIC(model_lin_acc_thick)
model_quad_acc_thick <- lm(m_m_acc ~ cortical_thickness_adsignature_dickson_c + cortical_thickness_adsignature_dickson_sq + age_spectro_c + sex + education + diagnostic_nick, data = MRS_S_ACC)
summary(model_quad_acc_thick)
AIC(model_quad_acc_thick)

anova(model_lin_acc_thick, model_quad_acc_thick)

# Model 4 m_m_Precuneus ~ thickness
model_lin_prec_thick <- lm(m_m_precuneus ~ cortical_thickness_adsignature_dickson_c  + age_spectro_c + sex + education + diagnostic_nick, data = MRS_S_Prec)
summary(model_lin_prec_thick)
AIC(model_lin_prec_thick)
model_quad_prec_thick <- lm(m_m_precuneus ~ cortical_thickness_adsignature_dickson_c + cortical_thickness_adsignature_dickson_sq + age_spectro_c + sex + education + diagnostic_nick, data = MRS_S_Prec)
summary(model_quad_prec_thick)
AIC(model_quad_prec_thick)

anova(model_lin_prec_thick, model_quad_prec_thick)


#### Activaiton ######
names(MRS_A_ACC)
names(MRS_A_Prec)
## Model 5 ACC ~ activaiton parietal
model_lin_acc_sup_act_rev <- lm(m_m_acc ~ activation_parietal_sup_l + age_spectro_c + sex + education + diagnostic_nick, data = MRS_A_ACC)
summary(model_lin_acc_sup_act_rev)
AIC(model_lin_acc_sup_act_rev)
model_quad_acc_sup_act_rev <- lm(m_m_acc ~ activation_parietal_sup_l + activation_parietal_sup_l_sq + age_spectro_c + sex + education + diagnostic_nick, data = MRS_A_ACC)
summary(model_quad_acc_sup_act_rev)
AIC(model_quad_acc_sup_act_rev)

anova(model_lin_acc_sup_act_rev, model_quad_acc_sup_act_rev)

## Model 6 Precuneus ~ activaiton parietal
model_lin_prec_sup_act_rev <- lm(m_m_precuneus ~ activation_parietal_sup_l + age_spectro_c + sex + education + diagnostic_nick, data = MRS_A_Prec)
summary(model_lin_prec_sup_act_rev)
AIC(model_lin_prec_sup_act_rev)
model_quad_prec_sup_act_rev <- lm(m_m_precuneus ~ activation_parietal_sup_l + activation_parietal_sup_l_sq + age_spectro_c + sex + education + diagnostic_nick, data = MRS_A_Prec)
summary(model_quad_prec_sup_act_rev)
AIC(model_quad_prec_sup_act_rev)

anova(model_lin_prec_sup_act_rev, model_quad_prec_sup_act_rev)

# Model 7 acc ~ hipp activation mean
model_lin_acc_hipp_act <- lm(m_m_acc ~ hipp_mean_act_c + age_spectro_c + sex + education + diagnostic_nick , data = MRS_A_ACC)
summary(model_lin_acc_hipp_act)
AIC(model_lin_acc_hipp_act)
model_quad_acc_hipp_act <- lm(m_m_acc ~ hipp_mean_act_c + hipp_mean_act_sq + age_spectro_c + sex + education + diagnostic_nick, data = MRS_A_ACC)
summary(model_quad_acc_hipp_act)
AIC(model_quad_acc_hipp_act)

anova(model_lin_acc_hipp_act, model_quad_acc_hipp_act)

# Model 8 precuneus  ~ hipp activation mean
model_lin_prec_hipp_act <- lm(m_m_precuneus ~ hipp_mean_act_c + age_spectro_c + sex + education + diagnostic_nick, data = MRS_A_Prec)
summary(model_lin_prec_hipp_act)
AIC(model_lin_prec_hipp_act)
model_quad_prec_hipp_act <- lm(m_m_precuneus ~ hipp_mean_act_c + hipp_mean_act_sq + age_spectro_c + sex + education + diagnostic_nick, data = MRS_A_Prec)
summary(model_quad_prec_hipp_act)
AIC(model_quad_prec_hipp_act)

anova(model_lin_prec_hipp_act, model_quad_prec_hipp_act)


### Memory ###

# Model 9 memoria ~ ACC Glu
model_lin_acc_memor <- lm(memoria_libre_correcte ~ m_m_acc_c + age_spectro_c + sex + education + diagnostic_nick, data = MRS_M_ACC)
summary(model_lin_acc_memor)
AIC(model_lin_acc_memor)
model_quad_acc_memor <- lm(memoria_libre_correcte ~ m_m_acc_c + m_m_acc_sq + age_spectro_c + sex + education + diagnostic_nick, data = MRS_M_ACC)
summary(model_quad_acc_memor)
AIC(model_quad_acc_memor)

anova(model_lin_acc_memor, model_quad_acc_memor)

# Model 10 memoria ~ Precuneus Glu
model_lin_prec_memor <- lm(memoria_libre_correcte ~ m_m_precuneus_c  + age_spectro_c + sex + education + diagnostic_nick, data = MRS_M_Prec)
summary(model_lin_prec_memor)
AIC(model_lin_prec_memor)
model_quad_prec_memor <- lm(memoria_libre_correcte ~ m_m_precuneus_c + m_m_precuneus_sq   + age_spectro_c + sex + education + diagnostic_nick, data = MRS_M_Prec)
summary(model_quad_prec_memor)
AIC(model_quad_prec_memor)

anova(model_lin_prec_memor, model_quad_prec_memor)

# Model 11 face name ~  ACC glu
model_lin_acc_facename <- lm(face_name_rappel_differe_spectro ~ m_m_acc_c + age_spectro_c + sex + education + diagnostic_nick, data = MRS_M_ACC)
summary(model_lin_acc_facename)
AIC(model_lin_acc_facename)
model_quad_acc_facename <- lm(face_name_rappel_differe_spectro ~ m_m_acc_c + m_m_acc_sq + age_spectro_c + sex + education + diagnostic_nick, data = MRS_M_ACC)
summary(model_quad_acc_facename)
AIC(model_quad_acc_facename)

# Model 12 face name ~ Precuneus glu
model_lin_prec_memor <- lm(face_name_rappel_differe_spectro ~ m_m_precuneus_c + age_spectro_c + sex + education + diagnostic_nick, data = MRS_M_Prec)
summary(model_lin_prec_memor)
AIC(model_lin_prec_memor)
model_quad_prec_memor <- lm(face_name_rappel_differe_spectro ~ m_m_precuneus_c + m_m_precuneus_sq + age_spectro_c + sex + education + diagnostic_nick, data = MRS_M_Prec)
summary(model_quad_prec_memor)
AIC(model_quad_prec_memor)


############## Mediation analysis  (controlling for covariates)

# Mean Structure -> Mean Glu -> Mean Activation
# M ~ X
model_m_comb <- lm(mean_glu ~ mean_structure + age_spectro_c + sex + education + diagnostic_nick, data = mediation_combined)
summary(model_m_comb)

# Y ~ X + M
model_y_comb <- lm(mean_activation ~ mean_structure + mean_glu + age_spectro_c + sex + education + diagnostic_nick, data = mediation_combined)
summary(model_y_comb)

# Mediate
med_comb <- mediate(model_m_comb, model_y_comb, 
                    treat = "mean_structure", mediator = "mean_glu", 
                    boot = TRUE, sims = 1000)
summary(med_comb)
plot(med_comb)



## Contrloing for group (diagnostic)?
#Models by subgroup
#model_quad_prec_thick_CU <- lm(m_m_precuneus ~ cortical_thickness_adsignature_dickson_c + cortical_thickness_adsignature_dickson_sq + age_spectro_c + sex + education, data = subset(MRS_S_Prec, diagnostic_nick == "HC"))
#summary(model_quad_prec_thick_CU)
#model_quad_prec_thick_SCD <- lm(m_m_precuneus ~ cortical_thickness_adsignature_dickson_c + cortical_thickness_adsignature_dickson_sq + age_spectro_c + sex + education, data = subset(MRS_S_Prec, diagnostic_nick == "SCD+"))
#summary(model_quad_prec_thick_SCD)
#model_quad_prec_thick_MCI <- lm(m_m_precuneus ~ cortical_thickness_adsignature_dickson_c + cortical_thickness_adsignature_dickson_sq + age_spectro_c + sex + education, data = subset(MRS_S_Prec, diagnostic_nick == "MCI"))
#summary(model_quad_prec_thick_MCI)


# Interaction ?










######################################################### Table 2 ############################
#import_data_arsenii_20250609_hippR <- read_excel("C:/Users/okkam/Desktop/labo/article 1/rencontre_Sylvie_20251017/import_data_arsenii_20250609_hippR.xlsx")
MRS_outliers <- import_data_arsenii_20250609_hippR

# Clean the column name
MRS_outliers <- janitor::clean_names(MRS_outliers)
names(MRS_outliers)
# Convertir  en numérique
MRS_outliers$education <- as.numeric(MRS_outliers$education)
MRS_outliers$moca_corr_spectro <- as.numeric(MRS_outliers$moca_corr_spectro)
MRS_outliers$dprime_hit_fa <- as.numeric(MRS_outliers$dprime_hit_fa)
MRS_outliers$associative_memory_performance <- as.numeric(MRS_outliers$associative_memory_performance)
MRS_outliers$moca_corr_spectro <- as.numeric(MRS_outliers$moca_corr_spectro)
# Regrouper les groupes 'HC' et 'SCD' dans un seul groupe 'SCD'
MRS_outliers$diagnostic_nick <- ifelse(MRS_outliers$diagnostic_nick %in% c("HC", "SCD"), "HC", MRS_outliers$diagnostic_nick)

# Convertir diagnostic_nick en facteur avec l’ordre HC, SCD+, MCI
MRS_outliers$diagnostic_nick <- factor(MRS_outliers$diagnostic_nick, levels = c("HC", "SCD+", "MCI"))

# Convertir sex en facteur
MRS_outliers$sex <- as.factor(MRS_outliers$sex)

# Creer Hipp mean
MRS_outliers$hipp_mean <- (MRS_outliers$hip_l_nor_icv + MRS_outliers$righ_hip_vol)/2
MRS_outliers$hipp_mean_act <- (MRS_outliers$activation_hippocampus_l + MRS_outliers$activation_hippocampus_r)/2




#### Centrer et mettre au carré
# education
MRS_outliers$education_c  <- scale(MRS_outliers$education, center = TRUE, scale = FALSE)
MRS_outliers$education_sq <- MRS_outliers$education_c^2
# age_spectro
MRS_outliers$age_spectro_c  <- scale(MRS_outliers$age_spectro, center = TRUE, scale = FALSE)
MRS_outliers$age_spectro_sq <- MRS_outliers$age_spectro_c^2
# m_m_precuneus
MRS_outliers$m_m_precuneus_c  <- scale(MRS_outliers$m_m_precuneus, center = TRUE, scale = FALSE)
MRS_outliers$m_m_precuneus_sq <- MRS_outliers$m_m_precuneus_c^2
# m_m_acc
MRS_outliers$m_m_acc_c  <- scale(MRS_outliers$m_m_acc, center = TRUE, scale = FALSE)
MRS_outliers$m_m_acc_sq <- MRS_outliers$m_m_acc_c^2
# moca_corr_spectro
MRS_outliers$moca_corr_spectro_c  <- scale(MRS_outliers$moca_corr_spectro, center = TRUE, scale = FALSE)
MRS_outliers$moca_corr_spectro_sq <- MRS_outliers$moca_corr_spectro_c^2
# memoria_libre_correcte
MRS_outliers$memoria_libre_correcte_c  <- scale(MRS_outliers$memoria_libre_correcte, center = TRUE, scale = FALSE)
MRS_outliers$memoria_libre_correcte_sq <- MRS_outliers$memoria_libre_correcte_c^2
# face_name_rappel_differe_spectro
MRS_outliers$face_name_rappel_differe_spectro_c  <- scale(MRS_outliers$face_name_rappel_differe_spectro, center = TRUE, scale = FALSE)
MRS_outliers$face_name_rappel_differe_spectro_sq <- MRS_outliers$face_name_rappel_differe_spectro_c^2
# hip_l_nor_icv
MRS_outliers$hip_l_nor_icv_c  <- scale(MRS_outliers$hip_l_nor_icv, center = TRUE, scale = FALSE)
MRS_outliers$hip_l_nor_icv_sq <- MRS_outliers$hip_l_nor_icv_c^2
# hipp_mean
MRS_outliers$hipp_mean_c  <- scale(MRS_outliers$hipp_mean, center = TRUE, scale = FALSE)
MRS_outliers$hipp_mean_sq <- MRS_outliers$hipp_mean_c^2
# cortical_thickness_adsignature_dickson
MRS_outliers$cortical_thickness_adsignature_dickson_c  <- scale(MRS_outliers$cortical_thickness_adsignature_dickson, center = TRUE, scale = FALSE)
MRS_outliers$cortical_thickness_adsignature_dickson_sq <- MRS_outliers$cortical_thickness_adsignature_dickson_c^2
#  dprime_hit_fa
MRS_outliers$dprime_hit_fa_c  <- scale(MRS_outliers$dprime_hit_fa, center = TRUE, scale = FALSE)
MRS_outliers$dprime_hit_fa_sq <- MRS_outliers$dprime_hit_fa_c^2
#  associative_memory_performance
MRS_outliers$associative_memory_performance_c  <- scale(MRS_outliers$associative_memory_performance, center = TRUE, scale = FALSE)
MRS_outliers$associative_memory_performance_sq <- MRS_outliers$associative_memory_performance_c^2
#  associative_memory_performance_sylvie
MRS_outliers$associative_memory_performance_sylvie_c  <- scale(MRS_outliers$associative_memory_performance_sylvie, center = TRUE, scale = FALSE)
MRS_outliers$associative_memory_performance_sylvie_sq <- MRS_outliers$associative_memory_performance_sylvie_c^2
#  activation_hippocampus_l
MRS_outliers$activation_hippocampus_l_c  <- scale(MRS_outliers$activation_hippocampus_l, center = TRUE, scale = FALSE)
MRS_outliers$activation_hippocampus_l_sq <- MRS_outliers$activation_hippocampus_l_c^2
#  activation_hipp mean
MRS_outliers$hipp_mean_act_c  <- scale(MRS_outliers$hipp_mean_act, center = TRUE, scale = FALSE)
MRS_outliers$hipp_mean_act_sq <- MRS_outliers$hipp_mean_act_c^2
#  activation_parietal_sup_l 
MRS_outliers$activation_parietal_sup_l_c  <- scale(MRS_outliers$activation_parietal_sup_l, center = TRUE, scale = FALSE)
MRS_outliers$activation_parietal_sup_l_sq <- MRS_outliers$activation_parietal_sup_l_c^2
#  activation_temporal_inf_r
MRS_outliers$activation_temporal_inf_r_c  <- scale(MRS_outliers$activation_temporal_inf_r, center = TRUE, scale = FALSE)
MRS_outliers$activation_temporal_inf_r_sq <- MRS_outliers$activation_temporal_inf_r_c^2



#### Data Banks Subset Creation (Outlier / Raw Dataset) ####

# Structure
MRS_S_Prec_out <- MRS_outliers[, c("m_m_precuneus", "m_m_precuneus_c", "m_m_precuneus_sq", "hipp_mean", "hipp_mean_c", "hipp_mean_sq", "sex", "age_spectro_c",
                                   "education", "diagnostic_nick", "hip_l_nor_icv", "hip_l_nor_icv_c", "hip_l_nor_icv_sq", "cortical_thickness_adsignature_dickson",
                                   "cortical_thickness_adsignature_dickson_c", "cortical_thickness_adsignature_dickson_sq")] |> na.omit()

MRS_S_ACC_out <- MRS_outliers[, c("m_m_acc", "m_m_acc_c", "m_m_acc_sq", "m_m_precuneus_sq", "hipp_mean", "hipp_mean_c", "hipp_mean_sq", "sex", "age_spectro_c",
                                  "education", "diagnostic_nick", "hip_l_nor_icv", "hip_l_nor_icv_c", "hip_l_nor_icv_sq", "cortical_thickness_adsignature_dickson",
                                  "cortical_thickness_adsignature_dickson_c", "cortical_thickness_adsignature_dickson_sq")] |> na.omit()

# Activation
MRS_A_Prec_out <- MRS_outliers[, c("education", "diagnostic_nick", "m_m_precuneus", "m_m_precuneus_c", "m_m_precuneus_sq", "activation_hippocampus_l", "activation_hippocampus_l_c", "hipp_mean_act", "hipp_mean_act_c", "hipp_mean_act_sq",
                                   "activation_hippocampus_l_sq", "activation_parietal_sup_l", "activation_parietal_sup_l_c", "sex", "age_spectro_c",
                                   "activation_parietal_sup_l_sq")] |> na.omit()

MRS_A_ACC_out <- MRS_outliers[, c("m_m_acc", "m_m_acc_c", "m_m_acc_sq", "hipp_mean_act", "hipp_mean_act_c", "hipp_mean_act_sq",
                                  "education", "diagnostic_nick", "activation_hippocampus_l", "activation_hippocampus_l_c", "activation_hippocampus_l_sq", "sex", "age_spectro_c",
                                  "activation_parietal_sup_l", "activation_parietal_sup_l_c", "activation_parietal_sup_l_sq")] |> na.omit()

# Memory
MRS_M_Prec_out <- MRS_outliers[, c("m_m_precuneus", "m_m_precuneus_c", "m_m_precuneus_sq", "face_name_rappel_differe_spectro", "face_name_rappel_differe_spectro_c", "face_name_rappel_differe_spectro_sq",
                                   "education", "diagnostic_nick", "memoria_libre_correcte", "memoria_libre_correcte_c", "memoria_libre_correcte_sq", "sex", "age_spectro_c",
                                   "associative_memory_performance", "associative_memory_performance_c",
                                   "associative_memory_performance_sq")] |> na.omit()

MRS_M_ACC_out <- MRS_outliers[, c("m_m_acc", "m_m_acc_c", "m_m_acc_sq", "memoria_libre_correcte", "face_name_rappel_differe_spectro_c", "face_name_rappel_differe_spectro_sq",
                                  "education", "diagnostic_nick", "memoria_libre_correcte_c", "memoria_libre_correcte_sq", "face_name_rappel_differe_spectro", "sex", "age_spectro_c",
                                  "associative_memory_performance", "associative_memory_performance_c",
                                  "associative_memory_performance_sq")] |> na.omit()


library(robustbase)

#### Structure Models ####

# Model 1 m_m_ACC ~ hipp mean
model_rob_lin_acc_hipp_mean <- lmrob(m_m_acc ~ hipp_mean_c, data = MRS_S_ACC_out)
summary(model_rob_lin_acc_hipp_mean)

model_rob_quad_acc_hipp_mean <- lmrob(m_m_acc ~ hipp_mean_c + hipp_mean_sq, data = MRS_S_ACC_out)
summary(model_rob_quad_acc_hipp_mean)

anova(model_rob_lin_acc_hipp_mean, model_rob_quad_acc_hipp_mean, test = "Wald")


# Model 2 m_m_Precuneus ~ hipp mean 
model_rob_lin_prec_hipp_mean <- lmrob(m_m_precuneus ~ hipp_mean_c, data = MRS_S_Prec_out)
summary(model_rob_lin_prec_hipp_mean)

model_rob_quad_prec_hipp_mean <- lmrob(m_m_precuneus ~ hipp_mean_c + hipp_mean_sq, data = MRS_S_Prec_out)
summary(model_rob_quad_prec_hipp_mean)

anova(model_rob_lin_prec_hipp_mean, model_rob_quad_prec_hipp_mean, test = "Wald")


# Model 3 m_m_ACC ~ thickness
model_rob_lin_acc_thick <- lmrob(m_m_acc ~ cortical_thickness_adsignature_dickson_c, data = MRS_S_ACC_out)
summary(model_rob_lin_acc_thick)

model_rob_quad_acc_thick <- lmrob(m_m_acc ~ cortical_thickness_adsignature_dickson_c + cortical_thickness_adsignature_dickson_sq, data = MRS_S_ACC_out)
summary(model_rob_quad_acc_thick)

anova(model_rob_lin_acc_thick, model_rob_quad_acc_thick, test = "Wald")


# Model 4 m_m_Precuneus ~ thickness
model_rob_lin_prec_thick <- lmrob(m_m_precuneus ~ cortical_thickness_adsignature_dickson_c, data = MRS_S_Prec_out)
summary(model_rob_lin_prec_thick)

model_rob_quad_prec_thick <- lmrob(m_m_precuneus ~ cortical_thickness_adsignature_dickson_c + cortical_thickness_adsignature_dickson_sq, data = MRS_S_Prec_out)
summary(model_rob_quad_prec_thick)

anova(model_rob_lin_prec_thick, model_rob_quad_prec_thick, test = "Wald")


#### Activation Models ####

## Model 5 ACC ~ activation parietal
model_rob_lin_acc_sup_act_rev <- lmrob(m_m_acc ~ activation_parietal_sup_l, data = MRS_A_ACC_out)
summary(model_rob_lin_acc_sup_act_rev)

model_rob_quad_acc_sup_act_rev <- lmrob(m_m_acc ~ activation_parietal_sup_l + activation_parietal_sup_l_sq, data = MRS_A_ACC_out)
summary(model_rob_quad_acc_sup_act_rev)

anova(model_rob_lin_acc_sup_act_rev, model_rob_quad_acc_sup_act_rev, test = "Wald")


## Model 6 Precuneus ~ activation parietal
model_rob_lin_prec_sup_act_rev <- lmrob(m_m_precuneus ~ activation_parietal_sup_l, data = MRS_A_Prec_out)
summary(model_rob_lin_prec_sup_act_rev)

model_rob_quad_prec_sup_act_rev <- lmrob(m_m_precuneus ~ activation_parietal_sup_l + activation_parietal_sup_l_sq, data = MRS_A_Prec_out)
summary(model_rob_quad_prec_sup_act_rev)

anova(model_rob_lin_prec_sup_act_rev, model_rob_quad_prec_sup_act_rev, test = "Wald")


# Model 7 acc ~ hipp activation mean
model_rob_lin_acc_hipp_act <- lmrob(m_m_acc ~ hipp_mean_act_c, data = MRS_A_ACC_out)
summary(model_rob_lin_acc_hipp_act)

model_rob_quad_acc_hipp_act <- lmrob(m_m_acc ~ hipp_mean_act_c + hipp_mean_act_sq, data = MRS_A_ACC_out)
summary(model_rob_quad_acc_hipp_act)

anova(model_rob_lin_acc_hipp_act, model_rob_quad_acc_hipp_act, test = "Wald")


# Model 8 precuneus ~ hipp activation mean
model_rob_lin_prec_hipp_act <- lmrob(m_m_precuneus ~ hipp_mean_act_c, data = MRS_A_Prec_out)
summary(model_rob_lin_prec_hipp_act)

model_rob_quad_prec_hipp_act <- lmrob(m_m_precuneus ~ hipp_mean_act_c + hipp_mean_act_sq, data = MRS_A_Prec_out)
summary(model_rob_quad_prec_hipp_act)

anova(model_rob_lin_prec_hipp_act, model_rob_quad_prec_hipp_act, test = "Wald")


#### Memory Models ####

# Model 9 memoria ~ ACC Glu
model_rob_lin_acc_memor <- lmrob(memoria_libre_correcte ~ m_m_acc_c, data = MRS_M_ACC_out)
summary(model_rob_lin_acc_memor)

model_rob_quad_acc_memor <- lmrob(memoria_libre_correcte ~ m_m_acc_c + m_m_acc_sq, data = MRS_M_ACC_out)
summary(model_rob_quad_acc_memor)

anova(model_rob_lin_acc_memor, model_rob_quad_acc_memor, test = "Wald")


# Model 10 memoria ~ Precuneus Glu
model_rob_lin_prec_memor <- lmrob(memoria_libre_correcte ~ m_m_precuneus_c, data = MRS_M_Prec_out)
summary(model_rob_lin_prec_memor)

model_rob_quad_prec_memor <- lmrob(memoria_libre_correcte ~ m_m_precuneus_c + m_m_precuneus_sq, data = MRS_M_Prec_out)
summary(model_rob_quad_prec_memor)

anova(model_rob_lin_prec_memor, model_rob_quad_prec_memor, test = "Wald")


# Model 11 face name ~ ACC glu
model_rob_lin_acc_facename <- lmrob(face_name_rappel_differe_spectro ~ m_m_acc_c, data = MRS_M_ACC_out)
summary(model_rob_lin_acc_facename)

model_rob_quad_acc_facename <- lmrob(face_name_rappel_differe_spectro ~ m_m_acc_c + m_m_acc_sq, data = MRS_M_ACC_out)
summary(model_rob_quad_acc_facename)

anova(model_rob_lin_acc_facename, model_rob_quad_acc_facename, test = "Wald")


# Model 12 face name ~ Precuneus glu
model_rob_lin_prec_facename <- lmrob(face_name_rappel_differe_spectro ~ m_m_precuneus_c, data = MRS_M_Prec_out)
summary(model_rob_lin_prec_facename)

model_rob_quad_prec_facename <- lmrob(face_name_rappel_differe_spectro ~ m_m_precuneus_c + m_m_precuneus_sq, data = MRS_M_Prec_out)
summary(model_rob_quad_prec_facename)

anova(model_rob_lin_prec_facename, model_rob_quad_prec_facename, test = "Wald")



###### sensetivity analyses for memory quadratic: removing outliers ############

model_lin_prec_memor <- lm(memoria_libre_correcte ~ m_m_precuneus_c, data = MRS_M_Prec)
summary(model_lin_prec_memor)
AIC(model_lin_prec_memor)
model_quad_prec_memor <- lm(memoria_libre_correcte ~ m_m_precuneus_c + m_m_precuneus_sq, data = MRS_M_Prec)
summary(model_quad_prec_memor)
AIC(model_quad_prec_memor)

######### remove extreme 2 values ##########
MRS_M_Prec$prec_z <- scale(MRS_M_Prec$m_m_precuneus_c)

# order by z-score
ord <- order(MRS_M_Prec$prec_z)

# keep everything except the 2 smallest and 2 largest
keep_index <- ord[3:(length(ord) - 2)]

# trimmed dataset
MRS_M_Prec_trim4 <- MRS_M_Prec[keep_index, ]


model_lin_prec_memor_trim4  <- lm(memoria_libre_correcte ~ m_m_precuneus_c,
                                  data = MRS_M_Prec_trim4)
summary(model_lin_prec_memor_trim4)
AIC(model_lin_prec_memor_trim4)
model_quad_prec_memor_trim4 <- lm(memoria_libre_correcte ~ m_m_precuneus_c + m_m_precuneus_sq,
                                  data = MRS_M_Prec_trim4)
summary(model_quad_prec_memor_trim4)
AIC(model_quad_prec_memor_trim4)

anova(model_lin_prec_memor_trim4,model_quad_prec_memor_trim4)






############## Mediation analysis  
# Select raw columns
mediation_combined <- MRS_outliers[, c("m_m_precuneus", "m_m_acc", 
                                       "hipp_mean", "cortical_thickness_adsignature_dickson", 
                                       "hipp_mean_act", "activation_parietal_sup_l", "education", "diagnostic_nick",
                                       "sex", "age_spectro_c")] |> na.omit()
# Create composite "mean" variables
mediation_combined$mean_structure <- rowMeans(scale(mediation_combined[, c("hipp_mean", "cortical_thickness_adsignature_dickson")]), na.rm = TRUE)
mediation_combined$mean_activation <- rowMeans(scale(mediation_combined[, c("hipp_mean_act", "activation_parietal_sup_l")]), na.rm = TRUE)
mediation_combined$mean_glu <- rowMeans(scale(mediation_combined[, c("m_m_precuneus", "m_m_acc")]), na.rm = TRUE)



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
                    boot = TRUE, sims = 1000)
summary(med_comb)
plot(med_comb)



######## With WRS2 package ##########
library(WRS2)
x <- mediation_combined$mean_structure
y <- mediation_combined$mean_activation
med <- mediation_combined$mean_glu
zmediate_results <- ZYmediate(x, y, med, nboot = 1000, alpha = 0.05, kappa = 0.05)
zmediate_results










############# Table 3. Mediations sensetivity analyses ############################

library(mediation)

###############################################################################
# 1. Cortical Thickness -> ACC Glu -> Superior Parietal Activation
###############################################################################
model_m_1 <- lm(m_m_acc ~ cortical_thickness_adsignature_dickson, data = mediation_combined)
model_y_1 <- lm(activation_parietal_sup_l ~ cortical_thickness_adsignature_dickson + m_m_acc, data = mediation_combined)

med_1 <- mediate(model_m_1, model_y_1, treat = "cortical_thickness_adsignature_dickson", mediator = "m_m_acc", boot = TRUE, sims = 1000)
summary(med_1)
plot(med_1)


###############################################################################
# 2. Cortical Thickness -> ACC Glu -> Hippocampal Activation
###############################################################################
model_m_2 <- lm(m_m_acc ~ cortical_thickness_adsignature_dickson, data = mediation_combined)
model_y_2 <- lm(hipp_mean_act ~ cortical_thickness_adsignature_dickson + m_m_acc, data = mediation_combined)

med_2 <- mediate(model_m_2, model_y_2, treat = "cortical_thickness_adsignature_dickson", mediator = "m_m_acc", boot = TRUE, sims = 1000)
summary(med_2)
plot(med_2)


###############################################################################
# 3. Cortical Thickness -> Precuneus Glu -> Superior Parietal Activation
###############################################################################
model_m_3 <- lm(m_m_precuneus ~ cortical_thickness_adsignature_dickson, data = mediation_combined)
model_y_3 <- lm(activation_parietal_sup_l ~ cortical_thickness_adsignature_dickson + m_m_precuneus, data = mediation_combined)

med_3 <- mediate(model_m_3, model_y_3, treat = "cortical_thickness_adsignature_dickson", mediator = "m_m_precuneus", boot = TRUE, sims = 1000)
summary(med_3)
plot(med_3)


###############################################################################
# 4. Cortical Thickness -> Precuneus Glu -> Hippocampal Activation
###############################################################################
model_m_4 <- lm(m_m_precuneus ~ cortical_thickness_adsignature_dickson, data = mediation_combined)
model_y_4 <- lm(hipp_mean_act ~ cortical_thickness_adsignature_dickson + m_m_precuneus, data = mediation_combined)

med_4 <- mediate(model_m_4, model_y_4, treat = "cortical_thickness_adsignature_dickson", mediator = "m_m_precuneus", boot = TRUE, sims = 1000)
summary(med_4)
plot(med_4)


###############################################################################
# 5. Hippocampal Volume -> ACC Glu -> Superior Parietal Activation
###############################################################################
model_m_5 <- lm(m_m_acc ~ hipp_mean, data = mediation_combined)
model_y_5 <- lm(activation_parietal_sup_l ~ hipp_mean + m_m_acc, data = mediation_combined)

med_5 <- mediate(model_m_5, model_y_5, treat = "hipp_mean", mediator = "m_m_acc", boot = TRUE, sims = 1000)
summary(med_5)
plot(med_5)


###############################################################################
# 6. Hippocampal Volume -> ACC Glu -> Hippocampal Activation
###############################################################################
model_m_6 <- lm(m_m_acc ~ hipp_mean, data = mediation_combined)
model_y_6 <- lm(hipp_mean_act ~ hipp_mean + m_m_acc, data = mediation_combined)

med_6 <- mediate(model_m_6, model_y_6, treat = "hipp_mean", mediator = "m_m_acc", boot = TRUE, sims = 1000)
summary(med_6)
plot(med_6)


###############################################################################
# 7. Hippocampal Volume -> Precuneus Glu -> Superior Parietal Activation
###############################################################################
model_m_7 <- lm(m_m_precuneus ~ hipp_mean, data = mediation_combined)
model_y_7 <- lm(activation_parietal_sup_l ~ hipp_mean + m_m_precuneus, data = mediation_combined)

med_7 <- mediate(model_m_7, model_y_7, treat = "hipp_mean", mediator = "m_m_precuneus", boot = TRUE, sims = 1000)
summary(med_7)
plot(med_7)


###############################################################################
# 8. Hippocampal Volume -> Precuneus Glu -> Hippocampal Activation
###############################################################################
model_m_8 <- lm(m_m_precuneus ~ hipp_mean, data = mediation_combined)
model_y_8 <- lm(hipp_mean_act ~ hipp_mean + m_m_precuneus, data = mediation_combined)

med_8 <- mediate(model_m_8, model_y_8, treat = "hipp_mean", mediator = "m_m_precuneus", boot = TRUE, sims = 1000)
summary(med_8)
plot(med_8)





library(mediation)

########################   Squared     #######################################################

# Corrected column selection
mediation_sq <- MRS_full[, c("cortical_thickness_adsignature_dickson_c", 
                             "m_m_precuneus_c", 
                             "activation_parietal_sup_l_c", 
                             "hipp_mean_act_c")] |> na.omit()

# Flatten scale() matrices to 1D numeric vectors
mediation_sq$cortical_thickness_adsignature_dickson_c <- as.numeric(mediation_sq$cortical_thickness_adsignature_dickson_c)
mediation_sq$m_m_precuneus_c                          <- as.numeric(mediation_sq$m_m_precuneus_c)
mediation_sq$activation_parietal_sup_l_c                 <- as.numeric(mediation_sq$activation_parietal_sup_l_c)
mediation_sq$hipp_mean_act_c                          <- as.numeric(mediation_sq$hipp_mean_act_c)


sd_x <- sd(mediation_sq$cortical_thickness_adsignature_dickson_c)

###############################################################################
# Quadratic: Cortical Thickness -> Precuneus Glu -> Parietal Activation
###############################################################################
model_m_quad <- lm(m_m_precuneus_c ~ cortical_thickness_adsignature_dickson_c + I(cortical_thickness_adsignature_dickson_c^2), data = mediation_sq)
summary(model_m_quad)
model_y_quad <- lm(activation_parietal_sup_l_c ~ cortical_thickness_adsignature_dickson_c + I(cortical_thickness_adsignature_dickson_c^2) + m_m_precuneus_c, data = mediation_sq)
summary(model_y_quad)
med_quad <- mediate(model_m_quad, model_y_quad, treat = "cortical_thickness_adsignature_dickson_c", mediator = "m_m_precuneus_c", control.value = -sd_x, treat.value = sd_x, boot = TRUE, sims = 1000)
summary(med_quad)
plot(med_quad)


model_m_quad <- lm(m_m_precuneus_c ~ cortical_thickness_adsignature_dickson_c + I(cortical_thickness_adsignature_dickson_c^2), data = mediation_sq)
summary(model_m_quad)
model_y_quad <- lm(hipp_mean_act_c ~ cortical_thickness_adsignature_dickson_c + I(cortical_thickness_adsignature_dickson_c^2) + m_m_precuneus_c, data = mediation_sq)
summary(model_y_quad)
med_quad <- mediate(model_m_quad, model_y_quad, treat = "cortical_thickness_adsignature_dickson_c", mediator = "m_m_precuneus_c", control.value = -sd_x, treat.value = sd_x, boot = TRUE, sims = 5000)
summary(med_quad)
plot(med_quad)






mean_x <- mean(mediation_sq$cortical_thickness_adsignature_dickson_c)
sd_x   <- sd(mediation_sq$cortical_thickness_adsignature_dickson_c)
delta  <- 0.01  # Small step to approximate derivative

# Indirect effect at LOW Cortical Thickness (-1 SD)
med_low <- mediate(
  model_m_quad, model_y_quad, 
  treat = "cortical_thickness_adsignature_dickson_c", 
  mediator = "m_m_precuneus_c",
  control.value = (mean_x - sd_x) - (delta / 2),
  treat.value   = (mean_x - sd_x) + (delta / 2),
  boot = TRUE, sims = 1000
)

# Indirect effect at HIGH Cortical Thickness (+1 SD)
med_high <- mediate(
  model_m_quad, model_y_quad, 
  treat = "cortical_thickness_adsignature_dickson_c", 
  mediator = "m_m_precuneus_c",
  control.value = (mean_x + sd_x) - (delta / 2),
  treat.value   = (mean_x + sd_x) + (delta / 2),
  boot = TRUE, sims = 1000
)

summary(med_low)   # Shows indirect effect on left side of parabola
summary(med_high)  # Shows indirect effect on right side of parabola















