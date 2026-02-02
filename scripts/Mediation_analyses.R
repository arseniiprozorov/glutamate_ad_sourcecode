## Mediation analysis

install.packages("mediation")
library(mediation)
vignette("mediation")
names(MRS_full)

# Precuneus Dataset
MRS_Prec <- MRS_full[, c("m_m_precuneus", "hipp_mean", 
  "cortical_thickness_adsignature_dickson", "hipp_mean_act", 
  "activation_parietal_sup_l", "face_name_rappel_differe_spectro", 
  "memoria_libre_correcte", "sex", "age_spectro_c")] |> na.omit()

# ACC Dataset 
MRS_ACC <- MRS_full[, c("m_m_acc", "hipp_mean", "cortical_thickness_adsignature_dickson", 
  "hipp_mean_act", "activation_parietal_sup_l", "face_name_rappel_differe_spectro", 
  "memoria_libre_correcte", "sex", "age_spectro_c")] |> na.omit()

# Check the final sample sizes
print(paste("Subjects in MRS_Prec:", nrow(MRS_Prec)))
print(paste("Subjects in MRS_ACC:", nrow(MRS_ACC)))



cor.test(MRS_Prec$cortical_thickness_adsignature_dickson,MRS_Prec$hipp_mean)
cor.test(MRS_Prec$hipp_mean_act,MRS_Prec$activation_parietal_sup_l)


# Create a temporary subset of just the two columns
cols_struc_prec <- MRS_Prec[, c("cortical_thickness_adsignature_dickson", "hipp_mean")]
cols_struc_acc  <- MRS_ACC[, c("cortical_thickness_adsignature_dickson", "hipp_mean")]

MRS_Prec$mean_structure <- rowMeans(scale(cols_struc_prec), na.rm = TRUE)
MRS_ACC$mean_structure  <- rowMeans(scale(cols_struc_acc), na.rm = TRUE)
# 1. Select just the two Activation columns
cols_act_prec <- MRS_Prec[, c("hipp_mean_act", "activation_parietal_sup_l")]
cols_act_acc  <- MRS_ACC[, c("hipp_mean_act", "activation_parietal_sup_l")]
MRS_Prec$mean_activation <- rowMeans(scale(cols_act_prec), na.rm = TRUE)
MRS_ACC$mean_activation  <- rowMeans(scale(cols_act_acc), na.rm = TRUE)

MRS_full <- MRS_Prec$mean_activation 
MRS_full <- MRS_ACC$mean_structure


# Model 1 structure --> Glu (precuneus) --> activation
#  M <- X
mediation_struc_prec <- lm(m_m_precuneus ~ mean_structure, data = MRS_Prec)
summary(mediation_struc_prec)
#M Y˜X + M
mediation_act_struc <- lm(activation_parietal_sup_l ~ mean_structure + m_m_precuneus, data = MRS_Prec)
summary(mediation_act_struc)

#  mediate
mediation_stru_prec_act <- mediate(mediation_struc_prec, mediation_act_struc, 
  treat = "mean_structure", mediator = "m_m_precuneus", 
  boot = TRUE, sims = 1000)
plot(mediation_stru_prec_act)
summary(mediation_stru_prec_act)


# Model 2 structure --> Glu (acc) --> activation
#  M <- X
mediation_struc_acc <- lm(m_m_acc ~ mean_structure, data = MRS_ACC)
summary(mediation_struc_acc)
#M Y˜X + M
mediation_act_struc_acc <- lm(activation_parietal_sup_l ~ mean_structure + m_m_acc, data = MRS_ACC)
summary(mediation_act_struc_acc)

# mediate
mediation_stru_acc_act <- mediate(mediation_struc_acc, mediation_act_struc_acc, 
                                   treat = "mean_structure", mediator = "m_m_acc", 
                                   boot = TRUE, sims = 5000)
summary(mediation_stru_acc_act)





######################### Hypothesis 2 - compensation vs excitotoxcitiy ############
# 1. Select only the raw columns that actually exist in MRS_full
mediation_combined <- MRS_full[, c("m_m_precuneus", "m_m_acc", "hipp_mean", 
                                   "cortical_thickness_adsignature_dickson", 
                                   "hipp_mean_act", "activation_parietal_sup_l", 
                                   "sex", "age_spectro")] |> na.omit()

# 2. Create the composite "mean" variables inside this dataset
mediation_combined$mean_structure <- rowMeans(scale(mediation_combined[, c("hipp_mean", "cortical_thickness_adsignature_dickson")]), na.rm = TRUE)
mediation_combined$mean_activation <- rowMeans(scale(mediation_combined[, c("hipp_mean_act", "activation_parietal_sup_l")]), na.rm = TRUE)
mediation_combined$mean_glu <- rowMeans(scale(mediation_combined[, c("m_m_precuneus", "m_m_acc")]), na.rm = TRUE)

cor.test(mediation_combined$m_m_precuneus,mediation_combined$m_m_acc)


## 1. model combined
#  M <- X
mediation_struc_glu <- lm(mean_glu ~ mean_structure, data = mediation_combined)
summary(mediation_struc_glu)
#M Y˜X + M
mediation_act_struc <- lm(mean_activation ~ mean_structure + mean_glu, data = mediation_combined)
summary(mediation_act_struc)

#  mediate
mediation_stru_glu_act <- mediate(mediation_struc_glu, mediation_act_struc, 
                                   treat = "mean_structure", mediator = "mean_glu", 
                                   boot = TRUE, sims = 1000)
plot(mediation_stru_glu_act)
summary(mediation_stru_glu_act)



## 2 model parietal superior activaiton
#  M <- X
mediation_struc_glu <- lm(mean_glu ~ mean_structure, data = mediation_combined)
summary(mediation_struc_glu)
#M Y˜X + M
mediation_act_struc <- lm(activation_parietal_sup_l ~ mean_structure + mean_glu, data = mediation_combined)
summary(mediation_act_struc)

#  mediate
mediation_stru_glu_act <- mediate(mediation_struc_glu, mediation_act_struc, 
                                  treat = "mean_structure", mediator = "mean_glu", 
                                  boot = TRUE, sims = 1000)
plot(mediation_stru_glu_act)
summary(mediation_stru_glu_act)


## 3 model hippocmapal activation
#  M <- X
mediation_struc_glu <- lm(mean_glu ~ mean_structure, data = mediation_combined)
summary(mediation_struc_glu)
#M Y˜X + M
mediation_act_struc <- lm(hipp_mean_act ~ mean_structure + mean_glu, data = mediation_combined)
summary(mediation_act_struc)

#  mediate
mediation_stru_glu_act <- mediate(mediation_struc_glu, mediation_act_struc, 
                                  treat = "mean_structure", mediator = "mean_glu", 
                                  boot = TRUE, sims = 1000)
plot(mediation_stru_glu_act)
summary(mediation_stru_glu_act)




############# 4 models #########
# M <- X
m1_med <- lm(m_m_precuneus ~ mean_structure, data = mediation_combined)
# Y ~ X + M
m1_out <- lm(activation_parietal_sup_l ~ mean_structure + m_m_precuneus, data = mediation_combined)

# Mediate
med_prec_parietal <- mediate(m1_med, m1_out, 
                             treat = "mean_structure", mediator = "m_m_precuneus", 
                             boot = TRUE, sims = 1000)
summary(med_prec_parietal)
plot(med_prec_parietal)


# M <- X
m1_med <- lm(m_m_precuneus ~ mean_structure, data = mediation_combined)
# Y ~ X + M
m1_out <- lm(activation_parietal_sup_l ~ mean_structure + m_m_precuneus, data = mediation_combined)

# Mediate
med_prec_parietal <- mediate(m1_med, m1_out, 
                             treat = "mean_structure", mediator = "m_m_precuneus", 
                             boot = TRUE, sims = 1000)
summary(med_prec_parietal)
plot(med_prec_parietal)






# M <- X
m3_med <- lm(m_m_acc ~ mean_structure, data = mediation_combined)
# Y ~ X + M
m3_out <- lm(activation_parietal_sup_l ~ mean_structure + m_m_acc, data = mediation_combined)

# Mediate
med_acc_parietal <- mediate(m3_med, m3_out, 
                            treat = "mean_structure", mediator = "m_m_acc", 
                            boot = TRUE, sims = 1000)
summary(med_acc_parietal)
plot(med_acc_parietal)



# M <- X
m4_med <- lm(m_m_acc ~ mean_structure, data = mediation_combined)
# Y ~ X + M
m4_out <- lm(hipp_mean_act ~ mean_structure + m_m_acc, data = mediation_combined)

# Mediate
med_acc_hipp <- mediate(m4_med, m4_out, 
                        treat = "mean_structure", mediator = "m_m_acc", 
                        boot = TRUE, sims = 1000)
summary(med_acc_hipp)
plot(med_acc_hipp)




# Model M: Does Structure predict Precuneus Glu?
m_prec_med <- lm(m_m_precuneus ~ mean_structure, data = mediation_combined)

# Model Y: Do Structure + Precuneus Glu predict Hippocampal Activation?
m_prec_out <- lm(hipp_mean_act ~ mean_structure + m_m_precuneus, data = mediation_combined)

# Run Mediation
med_prec_hipp <- mediate(m_prec_med, m_prec_out, 
                         treat = "mean_structure", mediator = "m_m_precuneus", 
                         boot = TRUE, sims = 1000)


summary(med_prec_hipp)



# other hypotheses

# Model 3 Prec --> act --> stru
# M <- X (Does Glu predict Activation?)
mediation_act_prec <- lm(activation_parietal_sup_l ~ m_m_precuneus, data = MRS_Prec)
summary(mediation_act_prec)

# Y ~ X + M (Do Glu + Activation predict Structure?)
mediation_struc_act_prec <- lm(mean_structure ~ m_m_precuneus + activation_parietal_sup_l, data = MRS_Prec)
summary(mediation_struc_act_prec)

#  3 mediate
mediation_prec_act_struc <- mediate(mediation_act_prec, mediation_struc_act_prec, 
                                    treat = "m_m_precuneus", mediator = "activation_parietal_sup_l", 
                                    boot = TRUE, sims = 5000)
summary(mediation_prec_act_struc)

# Model 4 ACC --> act --> stru
mediation_act_acc <- lm(activation_parietal_sup_l ~ m_m_acc, data = MRS_ACC)
summary(mediation_act_acc)

# Y ~ X + M (Do Glu + Activation predict Structure?)
mediation_struc_act_acc <- lm(mean_structure ~ m_m_acc + activation_parietal_sup_l, data = MRS_ACC)
summary(mediation_struc_act_acc)

#  mediate
mediation_acc_act_struc <- mediate(mediation_act_acc, mediation_struc_act_acc, 
                                   treat = "m_m_acc", mediator = "activation_parietal_sup_l", 
                                   boot = TRUE, sims = 5000)
summary(mediation_acc_act_struc)





# --- Model 5: Structure --> Precuneus Glu --> Memory ---

# Step 1: Model M (Does Structure predict Glu?)
med_glu_prec_mem <- lm(m_m_precuneus ~ mean_structure, data = MRS_Prec)
summary(med_glu_prec_mem)

# Step 2: Model Y (Do Structure + Glu predict Memory?)
out_mem_prec <- lm(memoria_libre_correcte ~ mean_structure + m_m_precuneus, data = MRS_Prec)
summary(out_mem_prec)

# Step 3: Run Mediation
res_struc_prec_mem <- mediate(med_glu_prec_mem, out_mem_prec, 
                              treat = "mean_structure", mediator = "m_m_precuneus", 
                              boot = TRUE, sims = 5000)
summary(res_struc_prec_mem)


# --- Model 6: Structure --> ACC Glu --> Memory ---

# Step 1: Model M (Does Structure predict Glu?)
med_glu_acc_mem <- lm(m_m_acc ~ mean_structure, data = MRS_ACC)
summary(med_glu_acc_mem)

# Step 2: Model Y (Do Structure + Glu predict Memory?)
out_mem_acc <- lm(memoria_libre_correcte ~ mean_structure + m_m_acc, data = MRS_ACC)
summary(out_mem_acc)

# Step 3: Run Mediation
res_struc_acc_mem <- mediate(med_glu_acc_mem, out_mem_acc, 
                             treat = "mean_structure", mediator = "m_m_acc", 
                             boot = TRUE, sims = 5000)
summary(res_struc_acc_mem)




# ==============================================================================
# MODEL 7: Precuneus Glu --> Memory --> Structure
# ==============================================================================

# Step 1: Model M (Does Glu predict Memory?)
med_mem_prec <- lm(memoria_libre_correcte ~ m_m_precuneus, data = MRS_Prec)
summary(med_mem_prec)

# Step 2: Model Y (Do Glu + Memory predict Structure?)
out_struc_prec <- lm(mean_structure ~ m_m_precuneus + memoria_libre_correcte, data = MRS_Prec)
summary(out_struc_prec)

# Step 3: Run Mediation
res_glu_mem_struc_prec <- mediate(med_mem_prec, out_struc_prec, 
                                  treat = "m_m_precuneus", mediator = "memoria_libre_correcte", 
                                  boot = TRUE, sims = 5000)
summary(res_glu_mem_struc_prec)


# ==============================================================================
# MODEL 8: ACC Glu --> Memory --> Structure
# ==============================================================================

# Step 1: Model M (Does Glu predict Memory?)
med_mem_acc <- lm(memoria_libre_correcte ~ m_m_acc, data = MRS_ACC)
summary(med_mem_acc)

# Step 2: Model Y (Do Glu + Memory predict Structure?)
out_struc_acc <- lm(mean_structure ~ m_m_acc + memoria_libre_correcte, data = MRS_ACC)
summary(out_struc_acc)

# Step 3: Run Mediation
res_glu_mem_struc_acc <- mediate(med_mem_acc, out_struc_acc, 
                                 treat = "m_m_acc", mediator = "memoria_libre_correcte", 
                                 boot = TRUE, sims = 5000)
summary(res_glu_mem_struc_acc)

