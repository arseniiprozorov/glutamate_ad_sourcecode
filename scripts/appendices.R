## Supplementary Table 1 - Sensetivity analysis with normalisation of cortical thickness by surface 
# AD-SIGNATURE CORTICAL THICKNESS DERIVATION & SENSITIVITY FORMULAS
# ## Regions from FreeSurfer Analysis (Desikan-Killiany Atlas):
#   - Medial Temporal (entorhinal and parahippocampal cortices)
#   - Inferior Temporal Cortex
#   - Temporal Pole
#   - Superior Frontal Gyrus
#   - Inferior Parietal Lobule
#   - Superior Parietal Lobule
#   - Supramarginal Gyrus
#   - Precuneus
#   - Inferior Frontal Cortex (pars opercularis, pars orbitalis, and pars triangularis)
#   - Pericalcarine Cortex
#
# EXTRACTION & DERIVATION LOGIC:
# 1. Bilateral Averaging (Step 1):
#    For every individual region (ROI), the bilateral average thickness was 
#    computed first across hemispheres:
#    Thickness_ROI = (Left_ROI_Thickness + Right_ROI_Thickness) / 2
#
# 2. Sub-regional Composite Averaging (Step 2 - for parenthetical clusters):
#    For compound structures containing multiple sub-regions:
#    - Medial Temporal = Mean(Bilateral Entorhinal, Bilateral Parahippocampal)
#    - Inferior Frontal = Mean(Bilateral Pars Opercularis, Pars Orbitalis, Pars Triangularis)
#
# 3. Primary Metric (Unweighted Mean):
#    The primary AD-signature thickness composite was computed as the simple 
#    unweighted arithmetic average across all target ROI mean thickness values (in mm):
#    AD_Signature_Thickness_Unweighted = Mean(ROI_1, ROI_2, ..., ROI_k)
#
# 4. Sensitivity Metric (Surface-Area-Weighted Mean):
#    To control for regional surface area variations, the weighted metric 
#    was computed using total surface area weighting across constituent ROIs:
#    AD_Signature_Thickness_Weighted = 
#      sum(Thickness_ROI_i * SurfaceArea_ROI_i) / sum(SurfaceArea_ROI_i)






### Supplementary table 3 - quality control MRS

final_dataset_used_with_QC_ACC_Precuneus <- read_excel("C:/Users/okkam/Desktop/labo/article 1/NeuroImage Clinical/MRS/analysis_supplementary/final_dataset_used_with_QC_ACC_Precuneus.xlsx")
MRS_QC <- final_dataset_used_with_QC_ACC_Precuneus

# Clean the column name
MRS_QC <- janitor::clean_names(MRS_QC)
names(MRS_QC)

# Regrouper les groupes 'HC' et 'SCD' dans un seul groupe 'SCD'
MRS_QC$diagnostic_nick <- ifelse(MRS_QC$diagnostic_nick %in% c("HC", "SCD"), "HC", MRS_QC$diagnostic_nick)
MRS_QC$diagnostic_nick <- factor(MRS_QC$diagnostic_nick, levels = c("HC", "SCD+", "MCI"))


jmv::descriptives(data = MRS_QC, vars = vars(m_m_acc_2, m_m_crlb_acc, percent_crlb_acc, snr_acc, fwhm_acc, gm_acc, wm_acc, csf_acc), sd = TRUE, range = TRUE, skew = TRUE, kurt = TRUE)
jmv::descriptives(data = MRS_QC, vars = vars(m_m_precuneus_2, m_m_crlb_precuneus, percent_crlb_precuneus, snr_precuneus, fwhm_precuneus, gm_precuneus, wm_precuneus, csf_precuneus), sd = TRUE, range = TRUE, skew = TRUE, kurt = TRUE)


jmv::descriptives(data = MRS_QC, vars = vars(m_m_acc_2, m_m_crlb_acc, percent_crlb_acc, snr_acc, fwhm_acc, gm_acc, wm_acc, csf_acc), sd = TRUE, range = TRUE, skew = TRUE, kurt = TRUE, splitBy = "diagnostic_nick")
jmv::descriptives(data = MRS_QC, vars = vars(m_m_precuneus_2, m_m_crlb_precuneus, percent_crlb_precuneus, snr_precuneus, fwhm_precuneus, gm_precuneus, wm_precuneus, csf_precuneus), sd = TRUE, range = TRUE, skew = TRUE, kurt = TRUE, splitBy = "diagnostic_nick")

summary(aov(m_m_crlb_acc ~ diagnostic_nick, data = MRS_QC))
summary(aov(snr_acc ~ diagnostic_nick, data = MRS_QC))
summary(aov(fwhm_acc ~ diagnostic_nick, data = MRS_QC))
summary(aov(gm_acc ~ diagnostic_nick, data = MRS_QC))
summary(aov(wm_acc ~ diagnostic_nick, data = MRS_QC))
summary(aov(csf_acc ~ diagnostic_nick, data = MRS_QC))

summary(aov(m_m_crlb_precuneus ~ diagnostic_nick, data = MRS_QC))
summary(aov(snr_precuneus ~ diagnostic_nick, data = MRS_QC))
summary(aov(fwhm_precuneus ~ diagnostic_nick, data = MRS_QC))
summary(aov(gm_precuneus ~ diagnostic_nick, data = MRS_QC))
summary(aov(wm_precuneus ~ diagnostic_nick, data = MRS_QC))
summary(aov(csf_precuneus ~ diagnostic_nick, data = MRS_QC))








######################### Q-Q Plots (Residual Normality) #########################
#pdf("QQ_Plots_Models_1_to_10.pdf", width = 10, height = 5)
# Model 1: m_m_ACC ~ hipp mean
par(mfrow = c(1, 2))
plot(model_lin_acc_hipp_mean, which = 2, main = "Model 1 Linear: ACC ~ Hipp Vol")
plot(model_quad_acc_hipp_mean, which = 2, main = "Model 1 Quad: ACC ~ Hipp Vol")

# Model 2: m_m_Precuneus ~ hipp mean
par(mfrow = c(1, 2))
plot(model_lin_prec_hipp_mean, which = 2, main = "Model 2 Linear: Prec ~ Hipp Vol")
plot(model_quad_prec_hipp_mean, which = 2, main = "Model 2 Quad: Prec ~ Hipp Vol")

# Model 3: m_m_ACC ~ thickness
par(mfrow = c(1, 2))
plot(model_lin_acc_thick, which = 2, main = "Model 3 Linear: ACC ~ Thickness")
plot(model_quad_acc_thick, which = 2, main = "Model 3 Quad: ACC ~ Thickness")

# Model 4: m_m_Precuneus ~ thickness
par(mfrow = c(1, 2))
plot(model_lin_prec_thick, which = 2, main = "Model 4 Linear: Prec ~ Thickness")
plot(model_quad_prec_thick, which = 2, main = "Model 4 Quad: Prec ~ Thickness")
max(cooks.distance(model_quad_prec_thick))
cooks.distance(model_quad_prec_thick)[92]

# Model 5: ACC ~ parietal activation
par(mfrow = c(1, 2))
plot(model_lin_acc_sup_act_rev, which = 2, main = "Model 5 Linear: ACC ~ Parietal Act")
plot(model_quad_acc_sup_act_rev, which = 2, main = "Model 5 Quad: ACC ~ Parietal Act")

# Model 6: Precuneus ~ parietal activation
par(mfrow = c(1, 2))
plot(model_lin_prec_sup_act_rev, which = 2, main = "Model 6 Linear: Prec ~ Parietal Act")
plot(model_quad_prec_sup_act_rev, which = 2, main = "Model 6 Quad: Prec ~ Parietal Act")

# Model 7: ACC ~ hipp activation mean
par(mfrow = c(1, 2))
plot(model_lin_acc_hipp_act, which = 2, main = "Model 7 Linear: ACC ~ Hipp Act")
plot(model_quad_acc_hipp_act, which = 2, main = "Model 7 Quad: ACC ~ Hipp Act")

# Model 8: Precuneus ~ hipp activation mean
par(mfrow = c(1, 2))
plot(model_lin_prec_hipp_act, which = 2, main = "Model 8 Linear: Prec ~ Hipp Act")
plot(model_quad_prec_hipp_act, which = 2, main = "Model 8 Quad: Prec ~ Hipp Act")

# Model 9: memoria ~ ACC Glu
par(mfrow = c(1, 2))
plot(model_lin_acc_memor, which = 2, main = "Model 9 Linear: Memoria ~ ACC Glu")
plot(model_quad_acc_memor, which = 2, main = "Model 9 Quad: Memoria ~ ACC Glu")

# Model 10: memoria ~ Precuneus Glu
par(mfrow = c(1, 2))
plot(model_lin_prec_memor, which = 2, main = "Model 10 Linear: Memoria ~ Prec Glu")
plot(model_quad_prec_memor, which = 2, main = "Model 10 Quad: Memoria ~ Prec Glu")
max(cooks.distance(model_quad_prec_memor))
cooks.distance(model_quad_prec_memor)[92]

#dev.off()









################### Supplementary Table 2  - sensetivity 3.29 exlusions ############################

import_data_arsenii_20250609_hippR <- read_excel("C:/Users/okkam/Desktop/labo/article 1/import_data_arsenii_20261017.xlsx")
MRS_outliers <- import_data_arsenii_20250609_hippR

# Clean the column name
MRS_outliers <- janitor::clean_names(MRS_outliers)
names(MRS_outliers)

# Creer Hipp mean
MRS_outliers$hipp_mean <- (MRS_outliers$hip_l_nor_icv + MRS_outliers$righ_hip_vol)/2
MRS_outliers$hipp_mean_act <- (MRS_outliers$activation_hippocampus_l + MRS_outliers$activation_hippocampus_r)/2



#### Define the ±3.29 SD exclusion function (p < .001, two-tailed) ####
exclude_sd <- function(x, sd_multiplier = 3.29) {
  z <- as.numeric(scale(x, center = TRUE, scale = TRUE))
  x[!is.na(z) & abs(z) > sd_multiplier] <- NA
  return(x)}

MRS_outliers[c("m_m_precuneus", "m_m_acc", "memoria_libre_correcte","face_name_rappel_differe_spectro",
               "hipp_mean", "cortical_thickness_adsignature_dickson",
               "hipp_mean_act", "activation_parietal_sup_l")] <-
  lapply(MRS_outliers[c("m_m_precuneus", "m_m_acc", "memoria_libre_correcte","face_name_rappel_differe_spectro",
                        "hipp_mean", "cortical_thickness_adsignature_dickson",
                        "hipp_mean_act", "activation_parietal_sup_l")], exclude_sd)


#### Centrer et mettre au carré
# m_m_precuneus
MRS_outliers$m_m_precuneus_c  <- scale(MRS_outliers$m_m_precuneus, center = TRUE, scale = FALSE)
MRS_outliers$m_m_precuneus_sq <- MRS_outliers$m_m_precuneus_c^2
# m_m_acc
MRS_outliers$m_m_acc_c  <- scale(MRS_outliers$m_m_acc, center = TRUE, scale = FALSE)
MRS_outliers$m_m_acc_sq <- MRS_outliers$m_m_acc_c^2
# memoria_libre_correcte
MRS_outliers$memoria_libre_correcte_c  <- scale(MRS_outliers$memoria_libre_correcte, center = TRUE, scale = FALSE)
MRS_outliers$memoria_libre_correcte_sq <- MRS_outliers$memoria_libre_correcte_c^2
# face_name_rappel_differe_spectro
MRS_outliers$face_name_rappel_differe_spectro_c  <- scale(MRS_outliers$face_name_rappel_differe_spectro, center = TRUE, scale = FALSE)
MRS_outliers$face_name_rappel_differe_spectro_sq <- MRS_outliers$face_name_rappel_differe_spectro_c^2
# hipp_mean
MRS_outliers$hipp_mean_c  <- scale(MRS_outliers$hipp_mean, center = TRUE, scale = FALSE)
MRS_outliers$hipp_mean_sq <- MRS_outliers$hipp_mean_c^2
# cortical_thickness_adsignature_dickson
MRS_outliers$cortical_thickness_adsignature_dickson_c  <- scale(MRS_outliers$cortical_thickness_adsignature_dickson, center = TRUE, scale = FALSE)
MRS_outliers$cortical_thickness_adsignature_dickson_sq <- MRS_outliers$cortical_thickness_adsignature_dickson_c^2
#  activation_hipp mean
MRS_outliers$hipp_mean_act_c  <- scale(MRS_outliers$hipp_mean_act, center = TRUE, scale = FALSE)
MRS_outliers$hipp_mean_act_sq <- MRS_outliers$hipp_mean_act_c^2
#  activation_parietal_sup_l 
MRS_outliers$activation_parietal_sup_l_c  <- scale(MRS_outliers$activation_parietal_sup_l, center = TRUE, scale = FALSE)
MRS_outliers$activation_parietal_sup_l_sq <- MRS_outliers$activation_parietal_sup_l_c^2


#### Structure Models ####

# Model 1 m_m_ACC ~ hipp mean
model_rob_lin_acc_hipp_mean <- lm(m_m_acc ~ hipp_mean_c, data = MRS_outliers)
summary(model_rob_lin_acc_hipp_mean)

model_rob_quad_acc_hipp_mean <- lm(m_m_acc ~ hipp_mean_c + hipp_mean_sq, data = MRS_outliers)
summary(model_rob_quad_acc_hipp_mean)

AIC(model_rob_lin_acc_hipp_mean)
AIC(model_rob_quad_acc_hipp_mean)
anova(model_rob_lin_acc_hipp_mean, model_rob_quad_acc_hipp_mean)


# Model 2 m_m_Precuneus ~ hipp mean 
model_rob_lin_prec_hipp_mean <- lm(m_m_precuneus ~ hipp_mean_c, data = MRS_outliers)
summary(model_rob_lin_prec_hipp_mean)

model_rob_quad_prec_hipp_mean <- lm(m_m_precuneus ~ hipp_mean_c + hipp_mean_sq, data = MRS_outliers)
summary(model_rob_quad_prec_hipp_mean)

AIC(model_rob_lin_prec_hipp_mean)
AIC(model_rob_quad_prec_hipp_mean)
anova(model_rob_lin_prec_hipp_mean, model_rob_quad_prec_hipp_mean)


# Model 3 m_m_ACC ~ thickness
model_rob_lin_acc_thick <- lm(m_m_acc ~ cortical_thickness_adsignature_dickson_c, data = MRS_outliers)
summary(model_rob_lin_acc_thick)

model_rob_quad_acc_thick <- lm(m_m_acc ~ cortical_thickness_adsignature_dickson_c + cortical_thickness_adsignature_dickson_sq, data = MRS_outliers)
summary(model_rob_quad_acc_thick)


AIC(model_rob_lin_acc_thick)
AIC(model_rob_quad_acc_thick)
anova(model_rob_lin_acc_thick, model_rob_quad_acc_thick)


# Model 4 m_m_Precuneus ~ thickness
model_rob_lin_prec_thick <- lm(m_m_precuneus ~ cortical_thickness_adsignature_dickson_c, data = MRS_outliers)
summary(model_rob_lin_prec_thick)

model_rob_quad_prec_thick <- lm(m_m_precuneus ~ cortical_thickness_adsignature_dickson_c + cortical_thickness_adsignature_dickson_sq, data = MRS_outliers)
summary(model_rob_quad_prec_thick)

AIC(model_rob_lin_prec_thick)
AIC(model_rob_quad_prec_thick)
anova(model_rob_lin_prec_thick, model_rob_quad_prec_thick)


#### Activation Models ####

## Model 5 ACC ~ activation parietal
model_rob_lin_acc_sup_act_rev <- lm(m_m_acc ~ activation_parietal_sup_l_c, data = MRS_outliers)
summary(model_rob_lin_acc_sup_act_rev)

model_rob_quad_acc_sup_act_rev <- lm(m_m_acc ~ activation_parietal_sup_l_c + activation_parietal_sup_l_sq, data = MRS_outliers)
summary(model_rob_quad_acc_sup_act_rev)

AIC(model_rob_lin_acc_sup_act_rev)
AIC(model_rob_quad_acc_sup_act_rev)
anova(model_rob_lin_acc_sup_act_rev, model_rob_quad_acc_sup_act_rev)


## Model 6 Precuneus ~ activation parietal
model_rob_lin_prec_sup_act_rev <- lm(m_m_precuneus ~ activation_parietal_sup_l_c, data = MRS_outliers)
summary(model_rob_lin_prec_sup_act_rev)

model_rob_quad_prec_sup_act_rev <- lm(m_m_precuneus ~ activation_parietal_sup_l_c + activation_parietal_sup_l_sq, data = MRS_outliers)
summary(model_rob_quad_prec_sup_act_rev)


AIC(model_rob_lin_prec_sup_act_rev)
AIC(model_rob_quad_prec_sup_act_rev)
anova(model_rob_lin_prec_sup_act_rev, model_rob_quad_prec_sup_act_rev)


# Model 7 acc ~ hipp activation mean
model_rob_lin_acc_hipp_act <- lm(m_m_acc ~ hipp_mean_act_c, data = MRS_outliers)
summary(model_rob_lin_acc_hipp_act)

model_rob_quad_acc_hipp_act <- lm(m_m_acc ~ hipp_mean_act_c + hipp_mean_act_sq, data = MRS_outliers)
summary(model_rob_quad_acc_hipp_act)

AIC(model_rob_lin_acc_hipp_act)
AIC(model_rob_quad_acc_hipp_act)
anova(model_rob_lin_acc_hipp_act, model_rob_quad_acc_hipp_act)


# Model 8 precuneus ~ hipp activation mean
model_rob_lin_prec_hipp_act <- lm(m_m_precuneus ~ hipp_mean_act_c, data = MRS_outliers)
summary(model_rob_lin_prec_hipp_act)

model_rob_quad_prec_hipp_act <- lm(m_m_precuneus ~ hipp_mean_act_c + hipp_mean_act_sq, data = MRS_outliers)
summary(model_rob_quad_prec_hipp_act)

AIC(model_rob_lin_prec_hipp_act)
AIC(model_rob_quad_prec_hipp_act)
anova(model_rob_lin_prec_hipp_act, model_rob_quad_prec_hipp_act)


#### Memory Models ####

# Model 9 memoria ~ ACC Glu
model_rob_lin_acc_memor <- lm(memoria_libre_correcte ~ m_m_acc_c, data = MRS_outliers)
summary(model_rob_lin_acc_memor)

model_rob_quad_acc_memor <- lm(memoria_libre_correcte ~ m_m_acc_c + m_m_acc_sq, data = MRS_outliers)
summary(model_rob_quad_acc_memor)

AIC(model_rob_lin_acc_memor)
AIC(model_rob_quad_acc_memor)
anova(model_rob_lin_acc_memor, model_rob_quad_acc_memor)


# Model 10 memoria ~ Precuneus Glu
model_rob_lin_prec_memor <- lm(memoria_libre_correcte ~ m_m_precuneus_c, data = MRS_outliers)
summary(model_rob_lin_prec_memor)

model_rob_quad_prec_memor <- lm(memoria_libre_correcte ~ m_m_precuneus_c + m_m_precuneus_sq, data = MRS_outliers)
summary(model_rob_quad_prec_memor)

AIC(model_rob_lin_prec_memor)
AIC(model_rob_quad_prec_memor)
anova(model_rob_lin_prec_memor, model_rob_quad_prec_memor)







###   Controling for sex,  age,  education

#### structure  ######

# Model 1 m_m_ACC ~ hipp mean
model_lin_acc_hipp_mean <- lm(m_m_acc ~ hipp_mean_c + age_spectro_c + sex + education, data = MRS_S_ACC)
summary(model_lin_acc_hipp_mean)
AIC(model_lin_acc_hipp_mean)
model_quad_acc_hipp_mean <- lm(m_m_acc ~ hipp_mean_c + hipp_mean_sq + age_spectro_c + sex + education + diagnostic_nick, data = MRS_S_ACC)
summary(model_quad_acc_hipp_mean)
AIC(model_quad_acc_hipp_mean)

anova(model_lin_acc_hipp_mean, model_quad_acc_hipp_mean)

# Model 2 m_m_Precuneus ~ hipp mean 
model_lin_prec_hipp_mean <- lm(m_m_precuneus ~ hipp_mean_c + age_spectro_c + sex + education, data = MRS_S_Prec)
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








############## Mediation analysis ############### 
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
### Global model ####

############### thick --> Glu Prec ---> activation parietal ##############
# 1. Fit models with parietal activaiton
model_m_quad <- lm(m_m_precuneus_c ~ cortical_thickness_adsignature_dickson_c + 
                     I(cortical_thickness_adsignature_dickson_c^2), 
                   data = mediation_sq)
summary(model_m_quad)
model_y_quad <- lm(activation_parietal_sup_l_c ~ m_m_precuneus_c + cortical_thickness_adsignature_dickson_c + I(cortical_thickness_adsignature_dickson_c^2), 
                   data = mediation_sq)
summary(model_y_quad)


# Mediate
med_comb <- mediate(model_m_quad, model_y_quad, 
                    treat = "cortical_thickness_adsignature_dickson_c", mediator = "m_m_precuneus_c", 
                    boot = TRUE, sims = 5000)
summary(med_comb)
plot(med_comb)
# 2. Segment 1: Ascending Limb (Low to Moderate Atrophy: -1 SD to 0 SD)
med_ascending <- mediate(model_m_quad, model_y_quad, 
                         treat = "cortical_thickness_adsignature_dickson_c", 
                         mediator = "m_m_precuneus_c",
                         control.value = -1.0, treat.value = 0.0, 
                         boot = TRUE, sims = 1000)
summary(med_ascending)

# 3. Segment 2: Apex Region (Transition around sample mean: -0.5 SD to +0.5 SD)
med_apex <- mediate(model_m_quad, model_y_quad, 
                    treat = "cortical_thickness_adsignature_dickson_c", 
                    mediator = "m_m_precuneus_c",
                    control.value = -0.5, treat.value = 0.5, 
                    boot = TRUE, sims = 1000)
summary(med_apex)

# 4. Segment 3: Descending Limb (Moderate to Severe Atrophy: 0 SD to +1 SD)
med_descending <- mediate(model_m_quad, model_y_quad, 
                          treat = "cortical_thickness_adsignature_dickson_c", 
                          mediator = "m_m_precuneus_c",
                          control.value = 0.0, treat.value = 1.0, 
                          boot = TRUE, sims = 1000)
summary(med_descending)




############### thick --> Glu Prec ---> activation hippocampus ##############
# 1. Fit models with hippocampal activaiton
model_y_quad <- lm(hipp_mean_act_c ~ cortical_thickness_adsignature_dickson_c + 
                     I(cortical_thickness_adsignature_dickson_c^2) + 
                     m_m_precuneus_c, 
                   data = mediation_sq)
summary(model_y_quad)

# 2. Segment 1: Ascending Limb (Low to Moderate Atrophy: -1 SD to 0 SD)
med_ascending <- mediate(model_m_quad, model_y_quad, 
                         treat = "cortical_thickness_adsignature_dickson_c", 
                         mediator = "m_m_precuneus_c",
                         control.value = -1.0, treat.value = 0.0, 
                         boot = TRUE, sims = 1000)
summary(med_ascending)

# 3. Segment 2: Apex Region (Transition around sample mean: -0.5 SD to +0.5 SD)
med_apex <- mediate(model_m_quad, model_y_quad, 
                    treat = "cortical_thickness_adsignature_dickson_c", 
                    mediator = "m_m_precuneus_c",
                    control.value = -0.5, treat.value = 0.5, 
                    boot = TRUE, sims = 1000)
summary(med_apex)

# 4. Segment 3: Descending Limb (Moderate to Severe Atrophy: 0 SD to +1 SD)
med_descending <- mediate(model_m_quad, model_y_quad, 
                          treat = "cortical_thickness_adsignature_dickson_c", 
                          mediator = "m_m_precuneus_c",
                          control.value = 0.0, treat.value = 1.0, 
                          boot = TRUE, sims = 1000)
summary(med_descending)







library(lavaan)

# 1. Ensure quadratic term exists in dataframe
mediation_sq$cortical_thickness_sq <- mediation_sq$cortical_thickness_adsignature_dickson_c^2

# 2. Define SEM syntax with labeled paths
sem_model <- '
  # Mediator equation (M ~ X_lin + X_quad)
  m_m_precuneus_c ~ a1*cortical_thickness_adsignature_dickson_c + a2*cortical_thickness_sq

  # Outcome equation (Y ~ M + X_lin + X_quad)
  activation_parietal_sup_l_c ~ b*m_m_precuneus_c + c1*cortical_thickness_adsignature_dickson_c + c2*cortical_thickness_sq

  # Defined Parameters: Indirect Effects
  ind_linear    := a1 * b
  ind_quadratic := a2 * b
  
  # Defined Parameters: Total Indirect & Direct Effects
  ind_total     := ind_linear + ind_quadratic
  total_linear  := c1 + ind_linear
'

# 3. Fit SEM model with nonparametric bootstrapping for CI estimation
set.seed(123)
fit_sem <- sem(
  model = sem_model, 
  data = mediation_sq, 
  se = "bootstrap", 
  bootstrap = 5000
)

# 4. View results with standardized estimates and percentile CIs
summary(fit_sem, fit.measures = TRUE, standardize = TRUE, ci = TRUE, rsquare = TRUE)



 

