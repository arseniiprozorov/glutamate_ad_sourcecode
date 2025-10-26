###   Controling sex
## ANOVA
#  ACC
anova_acc <- aov(m_m_acc ~ diagnostic_nick + sex, data = MRS_full)
summary(anova_acc)
TukeyHSD(anova_acc)
eta_squared(anova_acc)          
# For Precuneus
anova_precuneus <- aov(m_m_precuneus ~ diagnostic_nick + sex, data = MRS_full)
summary(anova_precuneus)
TukeyHSD(anova_precuneus)
eta_squared(anova_precuneus)   



# Model 1 m_m_ACC ~ hipp mean + sex
model_lin_acc_hipp_mean <- lm(m_m_acc ~ hipp_mean_c + sex, data = MRS_S_ACC)
summary(model_lin_acc_hipp_mean)
AIC(model_lin_acc_hipp_mean)
model_quad_acc_hipp_mean <- lm(m_m_acc ~ hipp_mean_c + hipp_mean_sq + sex, data = MRS_S_ACC)
summary(model_quad_acc_hipp_mean)
AIC(model_quad_acc_hipp_mean)

# Model 2 m_m_Precuneus ~ hipp mean + sex
model_lin_prec_hipp_mean <- lm(m_m_precuneus ~ hipp_mean_c + sex, data = MRS_S_Prec)
summary(model_lin_prec_hipp_mean)
AIC(model_lin_prec_hipp_mean)
model_quad_prec_hipp_mean <- lm(m_m_precuneus ~ hipp_mean_c + hipp_mean_sq + sex, data = MRS_S_Prec)
summary(model_quad_prec_hipp_mean)
AIC(model_quad_prec_hipp_mean)

# Model 3 m_m_ACC ~ thickness + sex
model_lin_acc_thick <- lm(m_m_acc ~ cortical_thickness_adsignature_dickson_c + sex, data = MRS_S_ACC)
summary(model_lin_acc_thick)
AIC(model_lin_acc_thick)
model_quad_acc_thick <- lm(m_m_acc ~ cortical_thickness_adsignature_dickson_c + cortical_thickness_adsignature_dickson_sq + sex, data = MRS_S_ACC)
summary(model_quad_acc_thick)
AIC(model_quad_acc_thick)

# Model 4 m_m_Precuneus ~ thickness + sex
model_lin_prec_thick <- lm(m_m_precuneus ~ cortical_thickness_adsignature_dickson_c + sex, data = MRS_S_Prec)
summary(model_lin_prec_thick)
AIC(model_lin_prec_thick)
model_quad_prec_thick <- lm(m_m_precuneus ~ cortical_thickness_adsignature_dickson_c + cortical_thickness_adsignature_dickson_sq + sex, data = MRS_S_Prec)
summary(model_quad_prec_thick)
AIC(model_quad_prec_thick)

#### Activation ######

## Model 5 ACC ~ activation parietal + sex
model_lin_acc_sup_act_rev <- lm(m_m_acc_c ~ activation_parietal_sup_l + sex, data = MRS_A_ACC)
summary(model_lin_acc_sup_act_rev)
AIC(model_lin_acc_sup_act_rev)

model_quad_acc_sup_act_rev <- lm(m_m_acc_c ~ activation_parietal_sup_l + activation_parietal_sup_l_sq + sex, data = MRS_A_ACC)
summary(model_quad_acc_sup_act_rev)
AIC(model_quad_acc_sup_act_rev)

## Model 6 Precuneus ~ activation parietal + sex
model_lin_prec_sup_act_rev <- lm(m_m_precuneus_c ~ activation_parietal_sup_l + sex, data = MRS_A_Prec)
summary(model_lin_prec_sup_act_rev)
AIC(model_lin_prec_sup_act_rev)

model_quad_prec_sup_act_rev <- lm(m_m_precuneus_c ~ activation_parietal_sup_l + activation_parietal_sup_l_sq + sex, data = MRS_A_Prec)
summary(model_quad_prec_sup_act_rev)
AIC(model_quad_prec_sup_act_rev)

# Model 7 ACC ~ hippocampal activation + sex
model_lin_acc_hipp_act <- lm(m_m_acc ~ activation_hippocampus_l_c + sex, data = MRS_A_ACC)
summary(model_lin_acc_hipp_act)
AIC(model_lin_acc_hipp_act)
model_quad_acc_hipp_act <- lm(m_m_acc ~ activation_hippocampus_l_c + activation_hippocampus_l_sq + sex, data = MRS_A_ACC)
summary(model_quad_acc_hipp_act)
AIC(model_quad_acc_hipp_act)

# Model 8 Precuneus ~ hippocampal activation + sex
model_lin_prec_hipp_act <- lm(m_m_precuneus ~ activation_hippocampus_l_c + sex, data = MRS_A_Prec)
summary(model_lin_prec_hipp_act)
AIC(model_lin_prec_hipp_act)
model_quad_prec_hipp_act <- lm(m_m_precuneus ~ activation_hippocampus_l_c + activation_hippocampus_l_sq + sex, data = MRS_A_Prec)
summary(model_quad_prec_hipp_act)
AIC(model_quad_prec_hipp_act)

### Memory ###

# Model 9 memoria ~ ACC Glu + sex
model_lin_acc_memor <- lm(memoria_libre_correcte ~ m_m_acc_c + sex, data = MRS_M_ACC)
summary(model_lin_acc_memor)
AIC(model_lin_acc_memor)
model_quad_acc_memor <- lm(memoria_libre_correcte ~ m_m_acc_c + m_m_acc_sq + sex, data = MRS_M_ACC)
summary(model_quad_acc_memor)
AIC(model_quad_acc_memor)

# Model 10 memoria ~ Precuneus Glu + sex
model_lin_prec_memor <- lm(memoria_libre_correcte ~ m_m_precuneus_c + sex, data = MRS_M_Prec)
summary(model_lin_prec_memor)
AIC(model_lin_prec_memor)
model_quad_prec_memor <- lm(memoria_libre_correcte ~ m_m_precuneus_c + m_m_precuneus_sq + sex, data = MRS_M_Prec)
summary(model_quad_prec_memor)
AIC(model_quad_prec_memor)

# Model 11 face name ~ ACC Glu + sex
model_lin_acc_memor <- lm(face_name_rappel_differe_spectro ~ m_m_acc_c + sex, data = MRS_M_ACC)
summary(model_lin_acc_memor)
AIC(model_lin_acc_memor)
model_quad_acc_memor <- lm(face_name_rappel_differe_spectro ~ m_m_acc_c + m_m_acc_sq + sex, data = MRS_M_ACC)
summary(model_quad_acc_memor)
AIC(model_quad_acc_memor)

# Model 12 face name ~ Precuneus Glu + sex
model_lin_prec_memor <- lm(face_name_rappel_differe_spectro ~ m_m_precuneus_c + sex, data = MRS_M_Prec)
summary(model_lin_prec_memor)
AIC(model_lin_prec_memor)
model_quad_prec_memor <- lm(face_name_rappel_differe_spectro ~ m_m_precuneus_c + m_m_precuneus_sq + sex, data = MRS_M_Prec)
summary(model_quad_prec_memor)
AIC(model_quad_prec_memor)
