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







###### sensetivity analyses for memory quadratic ############

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






# Model 9  memoria ~ ACC Glu  (unchanged dataset)
# =================================================
model_lin_acc_memor <- lm(
  memoria_libre_correcte ~ m_m_acc_c,
  data = MRS_M_ACC
)
summary(model_lin_acc_memor)
AIC(model_lin_acc_memor)

model_quad_acc_memor <- lm(
  memoria_libre_correcte ~ m_m_acc_c + m_m_acc_sq,
  data = MRS_M_ACC
)
summary(model_quad_acc_memor)
AIC(model_quad_acc_memor)


# =================================================
# Model 10 memoria ~ Precuneus Glu  (TRIMMED)
# =================================================
model_lin_prec_memor_trim4 <- lm(
  memoria_libre_correcte ~ m_m_precuneus_c,
  data = MRS_M_Prec_trim4
)
summary(model_lin_prec_memor_trim4)
AIC(model_lin_prec_memor_trim4)

model_quad_prec_memor_trim4 <- lm(
  memoria_libre_correcte ~ m_m_precuneus_c + m_m_precuneus_sq,
  data = MRS_M_Prec_trim4
)
summary(model_quad_prec_memor_trim4)
AIC(model_quad_prec_memor_trim4)









######################################################    Just for SCD+ and MCI groups ################################


MRS_SCD_MCI <- MRS_full[MRS_full$diagnostic_nick == "SCD+" | MRS_full$diagnostic_nick == "MCI", ]
MRS_SCD_MCI$diagnostic_nick <- droplevels(MRS_SCD_MCI$diagnostic_nick)


MRS_HC <- MRS_full[MRS_full$diagnostic_nick == "HC",]
MRS_HC$diagnostic_nick <- droplevels(MRS_HC$diagnostic_nick)
names(MRS_HC)


model_quad_age <- lm(m_m_precuneus ~ age_spectro_c + age_spectro_sq, data = MRS_HC)
summary(model_quad_age)


model_quad_age <- lm(m_m_precuneus ~ age_spectro_c + age_spectro_sq, data = MRS_SCD_MCI)
summary(model_quad_age)


model_mod_age <- lm(memoria_libre_correcte ~ m_m_precuneus_c + I(m_m_precuneus_c^2) + age_spectro_c, data = MRS_full)
summary(model_mod_age)



#### structure and memory ######


# Model 1 m_m_ACC ~ hipp mean
model_lin_acc_hipp_mean <- lm(m_m_acc ~ hipp_mean_c, data = MRS_HC)
summary(model_lin_acc_hipp_mean)
AIC(model_lin_acc_hipp_mean)
model_quad_acc_hipp_mean <- lm(m_m_acc ~ hipp_mean_c + hipp_mean_sq, data = MRS_HC)
summary(model_quad_acc_hipp_mean)
AIC(model_quad_acc_hipp_mean)


# Model 2 m_m_Precuneus ~ hipp mean 
model_lin_prec_hipp_mean <- lm(m_m_precuneus ~ hipp_mean_c, data = MRS_SCD_MCI)
summary(model_lin_prec_hipp_mean)
AIC(model_lin_prec_hipp_mean)
model_quad_prec_hipp_mean <- lm(m_m_precuneus ~ hipp_mean_c + hipp_mean_sq, data = MRS_HC)
summary(model_quad_prec_hipp_mean)
AIC(model_quad_prec_hipp_mean)


# Model 3 m_m_ACC ~ thickness
model_lin_acc_thick <- lm(m_m_acc ~ cortical_thickness_adsignature_dickson_c, data = MRS_SCD_MCI)
summary(model_lin_acc_thick)
AIC(model_lin_acc_thick)
model_quad_acc_thick <- lm(m_m_acc ~ cortical_thickness_adsignature_dickson_c + cortical_thickness_adsignature_dickson_sq, data = MRS_HC)
summary(model_quad_acc_thick)
AIC(model_quad_acc_thick)


# Model 4 m_m_Precuneus ~ thickness
model_lin_prec_thick <- lm(m_m_precuneus ~ cortical_thickness_adsignature_dickson_c, data = MRS_SCD_MCI)
summary(model_lin_prec_thick)
AIC(model_lin_prec_thick)
model_quad_prec_thick <- lm(m_m_precuneus ~ cortical_thickness_adsignature_dickson_c + cortical_thickness_adsignature_dickson_sq, data = MRS_HC)
summary(model_quad_prec_thick)
AIC(model_quad_prec_thick)


#### Activaiton ######
names(MRS_A_ACC)
names(MRS_A_Prec)
## Model 5 ACC ~ activaiton parietal
model_lin_acc_sup_act_rev <- lm(m_m_acc_c ~ activation_parietal_sup_l, data = MRS_SCD_MCI)
summary(model_lin_acc_sup_act_rev)
AIC(model_lin_acc_sup_act_rev)

model_quad_acc_sup_act_rev <- lm(m_m_acc_c ~ activation_parietal_sup_l + activation_parietal_sup_l_sq, data = MRS_HC)
summary(model_quad_acc_sup_act_rev)
AIC(model_quad_acc_sup_act_rev)

## Model 6 recuneus ~ activaiton parietal
model_lin_prec_sup_act_rev <- lm(m_m_precuneus_c ~ activation_parietal_sup_l, data = MRS_SCD_MCI)
summary(model_lin_prec_sup_act_rev)
AIC(model_lin_prec_sup_act_rev)

model_quad_prec_sup_act_rev <- lm(m_m_precuneus_c ~ activation_parietal_sup_l + activation_parietal_sup_l_sq, data = MRS_SCD_MCI)
summary(model_quad_prec_sup_act_rev)
AIC(model_quad_prec_sup_act_rev)


# Model 7 acc ~ hipp activation
model_lin_acc_hipp_act <- lm(m_m_acc ~ activation_hippocampus_l_c, data = MRS_SCD_MCI)
summary(model_lin_acc_hipp_act)
AIC(model_lin_acc_hipp_act)
model_quad_acc_hipp_act <- lm(m_m_acc ~ activation_hippocampus_l_c + activation_hippocampus_l_sq, data = MRS_SCD_MCI)
summary(model_quad_acc_hipp_act)
AIC(model_quad_acc_hipp_act)

# Model 8 precuneus  ~ hipp activation
model_lin_prec_hipp_act <- lm(m_m_precuneus ~ activation_hippocampus_l_c, data = MRS_SCD_MCI)
summary(model_lin_prec_hipp_act)
AIC(model_lin_prec_hipp_act)
model_quad_prec_hipp_act <- lm(m_m_precuneus ~ activation_hippocampus_l_c + activation_hippocampus_l_sq, data = MRS_SCD_MCI)
summary(model_quad_prec_hipp_act)
AIC(model_quad_prec_hipp_act)

# Model 7 acc ~ hipp activation mean
model_lin_acc_hipp_act <- lm(m_m_acc ~ hipp_mean_act_c, data = MRS_SCD_MCI)
summary(model_lin_acc_hipp_act)
AIC(model_lin_acc_hipp_act)
model_quad_acc_hipp_act <- lm(m_m_acc ~ hipp_mean_act_c + hipp_mean_act_sq, data = MRS_SCD_MCI)
summary(model_quad_acc_hipp_act)
AIC(model_quad_acc_hipp_act)

# Model 8 precuneus  ~ hipp activation mean
model_lin_prec_hipp_act <- lm(m_m_precuneus ~ hipp_mean_act_c, data = MRS_SCD_MCI)
summary(model_lin_prec_hipp_act)
AIC(model_lin_prec_hipp_act)
model_quad_prec_hipp_act <- lm(m_m_precuneus ~ hipp_mean_act_c + hipp_mean_act_sq, data = MRS_SCD_MCI)
summary(model_quad_prec_hipp_act)
AIC(model_quad_prec_hipp_act)

### Memory ###

# Model 9 memoria ~ ACC Glu
model_lin_acc_memor <- lm(memoria_libre_correcte ~ m_m_acc_c, data = MRS_HC)
summary(model_lin_acc_memor)
AIC(model_lin_acc_memor)
model_quad_acc_memor <- lm(memoria_libre_correcte ~ m_m_acc_c + m_m_acc_sq, data = MRS_HC)
summary(model_quad_acc_memor)
AIC(model_quad_acc_memor)

# Model 10 memoria ~ Precuneus Glu
model_lin_prec_memor <- lm(memoria_libre_correcte ~ m_m_precuneus_c, data = MRS_HC)
summary(model_lin_prec_memor)
AIC(model_lin_prec_memor)
model_quad_prec_memor <- lm(memoria_libre_correcte ~ m_m_precuneus_c + m_m_precuneus_sq, data = MRS_HC)
summary(model_quad_prec_memor)
AIC(model_quad_prec_memor)

# Model 11 face name ~  ACC glu
model_lin_acc_memor <- lm(face_name_rappel_differe_spectro ~ m_m_acc_c, data = MRS_HC)
summary(model_lin_acc_memor)
AIC(model_lin_acc_memor)
model_quad_acc_memor <- lm(face_name_rappel_differe_spectro ~ m_m_acc_c + m_m_acc_sq, data = MRS_HC)
summary(model_quad_acc_memor)
AIC(model_quad_acc_memor)

# Model 12 face name ~ Precuneus glu
model_lin_prec_memor <- lm(face_name_rappel_differe_spectro ~ m_m_precuneus_c, data = MRS_HC)
summary(model_lin_prec_memor)
AIC(model_lin_prec_memor)
model_quad_prec_memor <- lm(face_name_rappel_differe_spectro ~ m_m_precuneus_c + m_m_precuneus_sq, data = MRS_HC)
summary(model_quad_prec_memor)
AIC(model_quad_prec_memor)


# Model 13 associative memory ~  ACC glu
model_lin_acc_memor <- lm(associative_memory_performance ~ m_m_acc_c, data = MRS_SCD_MCI)
summary(model_lin_acc_memor)
AIC(model_lin_acc_memor)
model_quad_acc_memor <- lm(associative_memory_performance ~ m_m_acc_c + m_m_acc_sq, data = MRS_SCD_MCI)
summary(model_quad_acc_memor)
AIC(model_quad_acc_memor)

# Model 14 assosiative memory ~ Precuneus glu
model_lin_prec_memor <- lm(associative_memory_performance ~ m_m_precuneus_c, data = MRS_HC)
summary(model_lin_prec_memor)
AIC(model_lin_prec_memor)
model_quad_prec_memor <- lm(associative_memory_performance ~ m_m_precuneus_c + m_m_precuneus_sq, data = MRS_HC)
summary(model_quad_prec_memor)
AIC(model_quad_prec_memor)




