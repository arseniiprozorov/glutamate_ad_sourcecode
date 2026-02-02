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


# Models structure -- Glu -- activation
# Model 1 hipp -- prec -- parietal
#  M <- X
mediation_hipp_prec <- lm(m_m_precuneus ~ hipp_mean, data = MRS_Prec)
summary(mediation_hipp_prec)
#M Y˜X + M
mediation_parietal_hipp <- lm(activation_parietal_sup_l ~ hipp_mean + m_m_precuneus, data = MRS_Prec)
summary(mediation_parietal_hipp)

# Model 3 mediate
mediation_stru_glu_act <- mediate(mediation_hipp_prec, mediation_parietal_hipp, 
  treat = "hipp_mean", mediator = "m_m_precuneus", 
  boot = TRUE, sims = 5000)
summary(mediation_stru_glu_act)




