

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



