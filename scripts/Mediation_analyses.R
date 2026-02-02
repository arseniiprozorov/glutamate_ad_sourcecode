## Mediation analysis

install.packages("mediation")
library(mediation)

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
  boot = TRUE, sims = 5000)
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

