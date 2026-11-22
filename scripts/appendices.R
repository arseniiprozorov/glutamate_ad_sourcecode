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


import_data_arsenii_20250609_hippR <- read_excel("C:/Users/okkam/Desktop/labo/article 1/rencontre_Sylvie_20250923/surface/import_data_arsenii_20250609_hippR_surface.xlsx")
MRS_surface <- import_data_arsenii_20250609_hippR

# Clean the column name
MRS_surface <- janitor::clean_names(MRS_surface)
names(MRS_surface)

MRS_surface$thick_normalized_surface <- winsorize_iqr(MRS_surface$thick_normalized_surface, 1.5)
MRS_surface$m_m_acc <- winsorize_iqr(MRS_surface$m_m_acc, 1.5)
MRS_surface$m_m_precuneus <- winsorize_iqr(MRS_surface$m_m_precuneus, 1.5)

#jmv::descriptives(data = MRS_surface, vars = vars(thick_normalized_surface),sd = TRUE,range = TRUE,skew = TRUE,kurt = TRUE)

#### Centrer et mettre au carré
# cortical_thickness_adsignature_dickson
MRS_surface$thick_normalized_surface_c  <- scale(MRS_surface$thick_normalized_surface, center = TRUE, scale = FALSE)
MRS_surface$thick_normalized_surface_sq <- MRS_surface$thick_normalized_surface_c^2


### Analysis #####

# Model 3 m_m_ACC ~ thickness
model_lin_acc_thick_surf <- lm(m_m_acc ~ thick_normalized_surface_c, data = MRS_surface)
summary(model_lin_acc_thick_surf)
model_quad_acc_thick_surf <- lm(m_m_acc ~ thick_normalized_surface_c + thick_normalized_surface_sq, data = MRS_surface)
summary(model_quad_acc_thick_surf)

AIC(model_lin_acc_thick_surf)
AIC(model_quad_acc_thick_surf)
anova(model_lin_acc_thick_surf, model_quad_acc_thick_surf)

# Model 4 m_m_Precuneus ~ thickness
model_lin_prec_thick_surf <- lm(m_m_precuneus ~ thick_normalized_surface_c , data = MRS_surface)
summary(model_lin_prec_thick_surf)
model_quad_prec_thick_surf <- lm(m_m_precuneus ~ thick_normalized_surface_c + thick_normalized_surface_sq, data = MRS_surface)
summary(model_quad_prec_thick_surf)

AIC(model_lin_prec_thick_surf)
AIC(model_quad_prec_thick_surf)
anova(model_lin_prec_thick_surf, model_quad_prec_thick_surf)




## Supplementary table 2 - quality control MRS

# MRS exlusions
# chi-square
acc_status <- ifelse(is.na(MRS_full$m_m_acc), "Excluded", "Included")
tab_acc    <- table(MRS_full$diagnostic_nick, acc_status, useNA = "ifany")
tab_acc
chisq.test(tab_acc)
fisher.test(tab_acc)


precuneus_status <- ifelse(is.na(MRS_full$m_m_precuneus), "Excluded", "Included")
tab_precuneus    <- table(MRS_full$diagnostic_nick, precuneus_status, useNA = "ifany")
tab_precuneus
chisq.test(tab_precuneus)
fisher.test(tab_precuneus)



# Exlusions fMRI per group
# chi-square
act_parietal_status <- ifelse(is.na(MRS_full$activation_parietal_sup_l), "Excluded", "Included")
tab_parietal_status    <- table(MRS_full$diagnostic_nick, act_parietal_status, useNA = "ifany")
tab_parietal_status
chisq.test(tab_parietal_status)
fisher.test(tab_parietal_status)


act_hipp_status <- ifelse(is.na(MRS_full$hipp_mean_act), "Excluded", "Included")
tab_act_hipp   <- table(MRS_full$diagnostic_nick, act_hipp_status, useNA = "ifany")
tab_act_hipp
chisq.test(tab_act_hipp)
fisher.test(tab_act_hipp)




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
pdf("QQ_Plots_Models_1_to_10.pdf", width = 10, height = 5)
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

dev.off()









###   Controling for sex,  age,  education

#### structure  ######

# Model 1 m_m_ACC ~ hipp mean
model_lin_acc_hipp_mean <- lm(m_m_acc ~ hipp_mean_c + age_spectro_c + sex + education, data = MRS_full)
summary(model_lin_acc_hipp_mean)
AIC(model_lin_acc_hipp_mean)
model_quad_acc_hipp_mean <- lm(m_m_acc ~ hipp_mean_c + hipp_mean_sq + age_spectro_c + sex + education, data = MRS_full)
summary(model_quad_acc_hipp_mean)
AIC(model_quad_acc_hipp_mean)

anova(model_lin_acc_hipp_mean, model_quad_acc_hipp_mean)

# Model 2 m_m_Precuneus ~ hipp mean 
model_lin_prec_hipp_mean <- lm(m_m_precuneus ~ hipp_mean_c + age_spectro_c + sex + education, data = MRS_full)
summary(model_lin_prec_hipp_mean)
AIC(model_lin_prec_hipp_mean)
model_quad_prec_hipp_mean <- lm(m_m_precuneus ~ hipp_mean_c + hipp_mean_sq + age_spectro_c + sex + education, data = MRS_full)
summary(model_quad_prec_hipp_mean)
AIC(model_quad_prec_hipp_mean)

anova(model_lin_prec_hipp_mean, model_quad_prec_hipp_mean)

# Model 3 m_m_ACC ~ thickness
model_lin_acc_thick <- lm(m_m_acc ~ cortical_thickness_adsignature_dickson_c + age_spectro_c + sex + education, data = MRS_full)
summary(model_lin_acc_thick)
AIC(model_lin_acc_thick)
model_quad_acc_thick <- lm(m_m_acc ~ cortical_thickness_adsignature_dickson_c + cortical_thickness_adsignature_dickson_sq + age_spectro_c + sex + education, data = MRS_full)
summary(model_quad_acc_thick)
AIC(model_quad_acc_thick)

anova(model_lin_acc_thick, model_quad_acc_thick)

# Model 4 m_m_Precuneus ~ thickness
model_lin_prec_thick <- lm(m_m_precuneus ~ cortical_thickness_adsignature_dickson_c  + age_spectro_c + sex + education, data = MRS_full)
summary(model_lin_prec_thick)
AIC(model_lin_prec_thick)
model_quad_prec_thick <- lm(m_m_precuneus ~ cortical_thickness_adsignature_dickson_c + cortical_thickness_adsignature_dickson_sq + age_spectro_c + sex + education, data = MRS_full)
summary(model_quad_prec_thick)
AIC(model_quad_prec_thick)

anova(model_lin_prec_thick, model_quad_prec_thick)


#### Activaiton ######

## Model 5 ACC ~ activaiton parietal
model_lin_acc_sup_act_rev <- lm(m_m_acc ~ activation_parietal_sup_l_c + age_spectro_c + sex + education , data = MRS_full)
summary(model_lin_acc_sup_act_rev)
AIC(model_lin_acc_sup_act_rev)
model_quad_acc_sup_act_rev <- lm(m_m_acc ~ activation_parietal_sup_l_c + activation_parietal_sup_l_sq + age_spectro_c + sex + education, data = MRS_full)
summary(model_quad_acc_sup_act_rev)
AIC(model_quad_acc_sup_act_rev)

anova(model_lin_acc_sup_act_rev, model_quad_acc_sup_act_rev)

## Model 6 Precuneus ~ activaiton parietal
model_lin_prec_sup_act_rev <- lm(m_m_precuneus ~ activation_parietal_sup_l_c + age_spectro_c + sex + education, data = MRS_full)
summary(model_lin_prec_sup_act_rev)
AIC(model_lin_prec_sup_act_rev)
model_quad_prec_sup_act_rev <- lm(m_m_precuneus ~ activation_parietal_sup_l_c + activation_parietal_sup_l_sq + age_spectro_c + sex + education, data = MRS_full)
summary(model_quad_prec_sup_act_rev)
AIC(model_quad_prec_sup_act_rev)

anova(model_lin_prec_sup_act_rev, model_quad_prec_sup_act_rev)

# Model 7 acc ~ hipp activation mean
model_lin_acc_hipp_act <- lm(m_m_acc ~ hipp_mean_act_c + age_spectro_c + sex + education , data = MRS_full)
summary(model_lin_acc_hipp_act)
AIC(model_lin_acc_hipp_act)
model_quad_acc_hipp_act <- lm(m_m_acc ~ hipp_mean_act_c + hipp_mean_act_sq + age_spectro_c + sex + education, data = MRS_full)
summary(model_quad_acc_hipp_act)
AIC(model_quad_acc_hipp_act)

anova(model_lin_acc_hipp_act, model_quad_acc_hipp_act)

# Model 8 precuneus  ~ hipp activation mean
model_lin_prec_hipp_act <- lm(m_m_precuneus ~ hipp_mean_act_c + age_spectro_c + sex + education, data = MRS_full)
summary(model_lin_prec_hipp_act)
AIC(model_lin_prec_hipp_act)
model_quad_prec_hipp_act <- lm(m_m_precuneus ~ hipp_mean_act_c + hipp_mean_act_sq + age_spectro_c + sex + education, data = MRS_full)
summary(model_quad_prec_hipp_act)
AIC(model_quad_prec_hipp_act)

anova(model_lin_prec_hipp_act, model_quad_prec_hipp_act)


### Memory ###

# Model 9 memoria ~ ACC Glu
model_lin_acc_memor <- lm(memoria_libre_correcte ~ m_m_acc_c + age_spectro_c + sex + education, data = MRS_full)
summary(model_lin_acc_memor)
AIC(model_lin_acc_memor)
model_quad_acc_memor <- lm(memoria_libre_correcte ~ m_m_acc_c + m_m_acc_sq + age_spectro_c + sex + education , data = MRS_full)
summary(model_quad_acc_memor)
AIC(model_quad_acc_memor)

anova(model_lin_acc_memor, model_quad_acc_memor)

# Model 10 memoria ~ Precuneus Glu
model_lin_prec_memor <- lm(memoria_libre_correcte ~ m_m_precuneus_c  + age_spectro_c + sex + education, data = MRS_full)
summary(model_lin_prec_memor)
AIC(model_lin_prec_memor)
model_quad_prec_memor <- lm(memoria_libre_correcte ~ m_m_precuneus_c + m_m_precuneus_sq   + age_spectro_c + sex + education , data = MRS_full)
summary(model_quad_prec_memor)
AIC(model_quad_prec_memor)

anova(model_lin_prec_memor, model_quad_prec_memor)

# Model 11 face name ~  ACC glu
model_lin_acc_facename <- lm(face_name_rappel_differe_spectro ~ m_m_acc_c + age_spectro_c + sex + education, data = MRS_full)
summary(model_lin_acc_facename)
AIC(model_lin_acc_facename)
model_quad_acc_facename <- lm(face_name_rappel_differe_spectro ~ m_m_acc_c + m_m_acc_sq + age_spectro_c + sex + education , data = MRS_full)
summary(model_quad_acc_facename)
AIC(model_quad_acc_facename)

anova(model_lin_acc_facename, model_quad_acc_facename)

# Model 12 face name ~ Precuneus glu
model_lin_prec_facename <- lm(face_name_rappel_differe_spectro ~ m_m_precuneus_c + age_spectro_c + sex + education , data = MRS_full)
summary(model_lin_prec_facename)
AIC(model_lin_prec_facename)
model_quad_prec_facename <- lm(face_name_rappel_differe_spectro ~ m_m_precuneus_c + m_m_precuneus_sq + age_spectro_c + sex + education, data = MRS_full)
summary(model_quad_prec_facename)
AIC(model_quad_prec_facename)

anova(model_lin_prec_facename, model_quad_prec_facename)













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











############# Table 6. Mediations sensetivity analyses ############################
set.seed(1)
#sink("mediation_supplementary.txt") 
# A. Mean Structure -> Mean Glu -> Mean Activation + covariates
model_m_comb <- lm(mean_glu ~ mean_structure + age_spectro_c + sex + education, data = mediation_combined)
summary(model_m_comb)

model_y_comb <- lm(mean_activation ~ mean_structure + mean_glu + age_spectro_c + sex + education, data = mediation_combined)
summary(model_y_comb)

med_comb <- mediate(model_m_comb, model_y_comb, 
                    treat = "mean_structure", mediator = "mean_glu", 
                    boot = TRUE, sims = 5000)
summary(med_comb)
plot(med_comb)


# B. Individual Region-Pair Linear Models
# 1. Cortical Thickness -> ACC Glu -> Superior Parietal Activation
model_m_1 <- lm(acc_glu_z ~ thick_z, data = mediation_combined)
model_y_1 <- lm(parietal_act_z ~ thick_z + acc_glu_z, data = mediation_combined)

med_1 <- mediate(model_m_1, model_y_1, treat = "thick_z", mediator = "acc_glu_z", boot = TRUE, sims = 5000)
summary(med_1)
plot(med_1)

# 2. Cortical Thickness -> ACC Glu -> Hippocampal Activation
model_m_2 <- lm(acc_glu_z ~ thick_z, data = mediation_combined)
model_y_2 <- lm(hipp_act_z ~ thick_z + acc_glu_z, data = mediation_combined)

med_2 <- mediate(model_m_2, model_y_2, treat = "thick_z", mediator = "acc_glu_z", boot = TRUE, sims = 5000)
summary(med_2)
plot(med_2)

# 3. Cortical Thickness -> Precuneus Glu -> Superior Parietal Activation
model_m_3 <- lm(prec_glu_z ~ thick_z, data = mediation_combined)
model_y_3 <- lm(parietal_act_z ~ thick_z + prec_glu_z, data = mediation_combined)

med_3 <- mediate(model_m_3, model_y_3, treat = "thick_z", mediator = "prec_glu_z", boot = TRUE, sims = 5000)
summary(med_3)
plot(med_3)

# 4. Cortical Thickness -> Precuneus Glu -> Hippocampal Activation
model_m_4 <- lm(prec_glu_z ~ thick_z, data = mediation_combined)
model_y_4 <- lm(hipp_act_z ~ thick_z + prec_glu_z, data = mediation_combined)

med_4 <- mediate(model_m_4, model_y_4, treat = "thick_z", mediator = "prec_glu_z", boot = TRUE, sims = 5000)
summary(med_4)
plot(med_4)

# 5. Hippocampal Volume -> ACC Glu -> Superior Parietal Activation
model_m_5 <- lm(acc_glu_z ~ hipp_vol_z, data = mediation_combined)
model_y_5 <- lm(parietal_act_z ~ hipp_vol_z + acc_glu_z, data = mediation_combined)

med_5 <- mediate(model_m_5, model_y_5, treat = "hipp_vol_z", mediator = "acc_glu_z", boot = TRUE, sims = 5000)
summary(med_5)
plot(med_5)

# 6. Hippocampal Volume -> ACC Glu -> Hippocampal Activation
model_m_6 <- lm(acc_glu_z ~ hipp_vol_z, data = mediation_combined)
model_y_6 <- lm(hipp_act_z ~ hipp_vol_z + acc_glu_z, data = mediation_combined)

med_6 <- mediate(model_m_6, model_y_6, treat = "hipp_vol_z", mediator = "acc_glu_z", boot = TRUE, sims = 5000)
summary(med_6)
plot(med_6)

# 7. Hippocampal Volume -> Precuneus Glu -> Superior Parietal Activation
model_m_7 <- lm(prec_glu_z ~ hipp_vol_z, data = mediation_combined)
model_y_7 <- lm(parietal_act_z ~ hipp_vol_z + prec_glu_z, data = mediation_combined)

med_7 <- mediate(model_m_7, model_y_7, treat = "hipp_vol_z", mediator = "prec_glu_z", boot = TRUE, sims = 5000)
summary(med_7)
plot(med_7)

# 8. Hippocampal Volume -> Precuneus Glu -> Hippocampal Activation
model_m_8 <- lm(prec_glu_z ~ hipp_vol_z, data = mediation_combined)
model_y_8 <- lm(hipp_act_z ~ hipp_vol_z + prec_glu_z, data = mediation_combined)

med_8 <- mediate(model_m_8, model_y_8, treat = "hipp_vol_z", mediator = "prec_glu_z", boot = TRUE, sims = 5000)
summary(med_8)
plot(med_8)


# C. Quadratic Path-a Models (N = 88 complete-case subset)
mediation_sq <- mediation_combined[, c("prec_glu_z", 
                                       "thick_z", 
                                       "hipp_act_z", 
                                       "parietal_act_z")] |> na.omit()

# 9. Quadratic: Cortical Thickness -> Precuneus Glu -> Parietal Activation
model_m_9 <- lm(prec_glu_z ~ thick_z + I(thick_z^2), data = mediation_sq)
summary(model_m_9)
model_y_9 <- lm(parietal_act_z ~ thick_z + prec_glu_z, data = mediation_sq)
summary(model_y_9)
med_9 <- mediate(model_m_9, model_y_9, treat = "thick_z", mediator = "prec_glu_z", boot = TRUE, sims = 5000)
summary(med_9)
plot(med_9)

# 10. Quadratic: Cortical Thickness -> Precuneus Glu -> Hippocampal Activation
model_m_10 <- lm(prec_glu_z ~ thick_z + I(thick_z^2), data = mediation_sq)
summary(model_m_10)
model_y_10 <- lm(hipp_act_z ~ thick_z + prec_glu_z, data = mediation_sq)
summary(model_y_10)
med_10 <- mediate(model_m_10, model_y_10, treat = "thick_z", mediator = "prec_glu_z", boot = TRUE, sims = 5000)
summary(med_10)
plot(med_10)


#sink() 







## Supplementary figure 1S Violin plots 
# --- Long format for ACC and Precuneus (unchanged) ---
long_data_mM <- MRS_full %>%
  dplyr::select(diagnostic_nick, m_m_acc, m_m_precuneus) %>%
  dplyr::rename(ACC = m_m_acc, Precuneus = m_m_precuneus) %>%
  tidyr::pivot_longer(c(ACC, Precuneus),
                      names_to = "Region",
                      values_to = "Value")


# --- Helper: pointy violins + EXACT Tukey overlays (no caps) ---
make_violin_plot_A <- function(region_name, y_label,
                               y_min = 13, y_max = 16,
                               show_points = TRUE,
                               violin_bw = 0.35,
                               violin_adjust = 1,
                               headroom = 0.45,
                               box_halfwidth = 0.12,
                               median_width   = 0.26) {
  
  axis_pad  <- 0.02; box_eps <- 1e-6
  
  dat_long <- long_data_mM %>%
    dplyr::filter(Region == region_name) %>%
    dplyr::mutate(group_f = factor(diagnostic_nick, levels = c("HC","SCD+","MCI")),
                  xi_num  = as.numeric(group_f))
  
  tallest <- suppressWarnings(max(dat_long$Value, na.rm = TRUE))
  y_upper  <- max(y_max, tallest) + headroom
  y_lower  <- y_min - axis_pad
  
  # Tukey box/whisker stats (no artificial minimum whiskers)
  box_df <- dat_long %>%
    dplyr::group_by(group_f) %>%
    dplyr::summarise(
      xi   = as.numeric(first(factor(group_f, levels = c("HC","SCD+","MCI")))),
      vals = list(sort(stats::na.omit(Value))),
      .groups = "drop"
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      q1  = as.numeric(stats::quantile(vals, 0.25, names = FALSE, type = 7)),
      med = as.numeric(stats::quantile(vals, 0.50, names = FALSE, type = 7)),
      q3  = as.numeric(stats::quantile(vals, 0.75, names = FALSE, type = 7)),
      iqr = q3 - q1,
      loF = q1 - 1.5 * iqr,
      hiF = q3 + 1.5 * iqr,
      wlow  = if (length(vals) == 0) q1 else {
        vv <- vals[vals >= loF]
        if (length(vv)) min(vv) else q1
      },
      whigh = if (length(vals) == 0) q3 else {
        vv <- vals[vals <= hiF]
        if (length(vv)) max(vv) else q3
      },
      eps = 0.015,
      draw_q1 = if (iqr == 0) med - eps else q1,
      draw_q3 = if (iqr == 0) med + eps else q3,
      draw_q3d = if ((draw_q3 - draw_q1) < 1e-8) draw_q3 + box_eps else draw_q3,
      xmin = xi - box_halfwidth, xmax = xi + box_halfwidth,
      mmin = xi - median_width/2, mmax = xi + median_width/2
    ) %>% dplyr::ungroup()
  
  ggplot(dat_long, aes(x = group_f, y = Value, fill = group_f)) +
    geom_violin(trim = FALSE, scale = "width", width = 0.8,
                bw = violin_bw, adjust = violin_adjust,
                color = "black", linewidth = 0.4, na.rm = TRUE) +
    { if (show_points)
      geom_point(position = position_jitter(width = 0.07, height = 0),
                 size = 1.0, alpha = 0.45, shape = 16, stroke = 0, na.rm = TRUE)
      else NULL } +
    geom_rect(data = box_df,
              aes(xmin = xmin, xmax = xmax, ymin = draw_q1, ymax = draw_q3d),
              inherit.aes = FALSE, fill = NA, color = "black", linewidth = 0.65) +
    geom_segment(data = box_df,
                 aes(x = mmin, xend = mmax, y = med, yend = med),
                 inherit.aes = FALSE, linewidth = 0.75, color = "black") +
    geom_segment(data = box_df,
                 aes(x = xi, xend = xi, y = wlow, yend = draw_q1),
                 inherit.aes = FALSE, linewidth = 0.7, lineend = "butt", color = "black") +
    geom_segment(data = box_df,
                 aes(x = xi, xend = xi, y = draw_q3, yend = whigh),
                 inherit.aes = FALSE, linewidth = 0.7, lineend = "butt", color = "black") +
    scale_fill_manual(values = c("HC"="#D9D9D9","SCD+"="#8F8F8F","MCI"="#4A4A4A"),
                      name = "Group") +
    scale_x_discrete(labels = c("HC","SCD+","MCI"), drop = FALSE) +
    scale_y_continuous(breaks = seq(y_min, y_max, 1),
                       limits = c(y_lower, y_upper),
                       expand = expansion(mult = c(0, 0.01))) +
    coord_cartesian(clip = "off") +
    labs(x = "", y = y_label) +
    theme_classic(base_size = 14) +
    theme(
      axis.title.y = element_text(size = 16),
      axis.text.x  = element_text(size = 14),
      axis.text.y  = element_text(size = 14),
      plot.margin  = margin(t = 26, r = 18, b = 10, l = 22),
      legend.position = "bottom")}


# --- Build panels (no significance elements) ---
p_acc  <- make_violin_plot("ACC", "ACC Glu (mM)", y_min = 12, y_max = 16, show_points = TRUE)
p_prec <- make_violin_plot("Precuneus", "Precuneus Glu (mM)", y_min = 12, y_max = 16, show_points = TRUE)

panel <- p_acc + p_prec + patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "right")

# A/B tags
panel <- panel + patchwork::plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 12),
        plot.tag.position = c(0.02, 0.98))

# --- Export ---
out_dir <- "C:/Users/okkam/Desktop/labo/article 1/Brain/methodes_resultats"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ragg::agg_tiff(file.path(out_dir, "fig_violin_ACC_Precuneus_panel_1200dpi.tiff"),
               width = 7, height = 4, units = "in", res = 1200, compression = "lzw")
print(panel)
dev.off()









#### supplementary ####
# -----------------------------
# Data preparation
# -----------------------------
DF <- MRS_full %>%
  transmute(
    memoria_libre_correcte = as.numeric(memoria_libre_correcte),
    mmprec_z               = as.numeric(scale(m_m_precuneus)),
    group                  = factor(diagnostic_nick, levels = c("HC","SCD+","MCI"))
  ) %>%
  drop_na()

# -----------------------------
# Trim: remove 2 lowest and 2 highest z-scores
# -----------------------------
ord <- order(DF$mmprec_z)                   # order by z-score
keep_index <- ord[3:(length(ord) - 2)]      # keep all except 2 smallest & 2 largest

DF_trim4 <- DF[keep_index, ]

# -----------------------------
# Quadratic fit (trimmed data)
# -----------------------------
mod_quad_trim4 <- lm(memoria_libre_correcte ~ mmprec_z + I(mmprec_z^2),
                     data = DF_trim4)
summary(mod_quad_trim4)

# -----------------------------
# Axis setup
# -----------------------------
y_min    <- floor(min(DF_trim4$memoria_libre_correcte, na.rm = TRUE))
y_max    <- ceiling(max(DF_trim4$memoria_libre_correcte, na.rm = TRUE))
y_breaks <- seq(y_min, y_max, by = 1)

# -----------------------------
# Plot (non-inverted x-axis, extremes removed)
# -----------------------------
p_quad_trim4 <- ggplot(DF_trim4, aes(x = mmprec_z, y = memoria_libre_correcte)) +
  geom_point(aes(shape = group), color = "grey30", size = 1.8, alpha = 0.85) +
  stat_smooth(
    method = "lm", formula = y ~ x + I(x^2),
    se = FALSE, color = "black", linewidth = 0.9
  ) +
  scale_shape_manual(values = c(16, 17, 15), name = "Group",
                     labels = c("CU", "SCD+", "MCI")) +  # HC → CU if you want
  scale_y_continuous(
    breaks = y_breaks,
    limits = c(y_min, y_max),
    expand = c(0.02, 0.02)
  ) +
  labs(
    x = "Precuneus Glu (z-score)",
    y = "Memoria free word recall (correct)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 11),
    plot.margin  = margin(8, 8, 8, 8)
  )

# -----------------------------
# Export (1200-DPI TIFF)
# -----------------------------
out_dir <- "C:/Users/okkam/Desktop/labo/article 1/Brain/methodes_resultats"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ragg::agg_tiff(
  filename = file.path(out_dir, "fig_memoria_vs_precuneusGluZ_quad_trim4_1200dpi.tiff"),
  width = 6.5, height = 4.3, units = "in", res = 1200, compression = "lzw"
)
print(p_quad_trim4)
dev.off()

