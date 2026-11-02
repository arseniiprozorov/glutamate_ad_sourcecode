
## Glutamate Relates to Structural and Functional markers of Disease Severity in early Alzheimer’s disease  ######
## Arsenii Prozorov 
#Création d’une banque de données
import_data_arsenii_20250609_hippR <- read_excel("C:/Users/okkam/Desktop/labo/article 1/rencontre_Sylvie_20250923/surface/import_data_arsenii_20250609_hippR_surface.xlsx")
MRS_surface <- import_data_arsenii_20250609_hippR

# Clean the column name
MRS_surface <- janitor::clean_names(MRS_surface)
names(MRS_surface)
# Creer Hipp mean
MRS_surface$hipp_mean <- (MRS_surface$hip_l_nor_icv + MRS_surface$righ_hip_vol)/2
MRS_surface$hipp_mean_act <- (MRS_surface$activation_hippocampus_l + MRS_surface$activation_hippocampus_r)/2



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
MRS_surface$memoria_libre_correcte <- winsorize_iqr(MRS_surface$memoria_libre_correcte, 1.5)
MRS_surface$face_name_rappel_differe_spectro <- winsorize_iqr(MRS_surface$face_name_rappel_differe_spectro, 1.5)
MRS_surface$activation_parietal_sup_l <- winsorize_iqr(MRS_surface$activation_parietal_sup_l, 1.5)
MRS_surface$hipp_mean_act <- winsorize_iqr(MRS_surface$hipp_mean_act, 1.5)

MRS_surface$hipp_mean <- winsorize_iqr(MRS_surface$hipp_mean, 1.5)
MRS_surface$cortical_thickness_adsignature_dickson <- winsorize_iqr(MRS_surface$cortical_thickness_adsignature_dickson, 1.5)
MRS_surface$thick_normalized_surface <- winsorize_iqr(MRS_surface$thick_normalized_surface, 1.5)
MRS_surface$m_m_acc <- winsorize_iqr(MRS_surface$m_m_acc, 1.5)
MRS_surface$m_m_precuneus <- winsorize_iqr(MRS_surface$m_m_precuneus, 1.5)


#jmv::descriptives(data = MRS_surface, vars = vars(hipp_mean_act),sd = TRUE,range = TRUE,skew = TRUE,kurt = TRUE)


#### Centrer et mettre au carré
# m_m_precuneus
MRS_surface$m_m_precuneus_c  <- scale(MRS_surface$m_m_precuneus, center = TRUE, scale = FALSE)
MRS_surface$m_m_precuneus_sq <- MRS_surface$m_m_precuneus_c^2
# m_m_acc
MRS_surface$m_m_acc_c  <- scale(MRS_surface$m_m_acc, center = TRUE, scale = FALSE)
MRS_surface$m_m_acc_sq <- MRS_surface$m_m_acc_c^2
# memoria_libre_correcte
MRS_surface$memoria_libre_correcte_c  <- scale(MRS_surface$memoria_libre_correcte, center = TRUE, scale = FALSE)
MRS_surface$memoria_libre_correcte_sq <- MRS_surface$memoria_libre_correcte_c^2
# face_name_rappel_differe_spectro
MRS_surface$face_name_rappel_differe_spectro_c  <- scale(MRS_surface$face_name_rappel_differe_spectro, center = TRUE, scale = FALSE)
MRS_surface$face_name_rappel_differe_spectro_sq <- MRS_surface$face_name_rappel_differe_spectro_c^2
# hipp_mean
MRS_surface$hipp_mean_c  <- scale(MRS_surface$hipp_mean, center = TRUE, scale = FALSE)
MRS_surface$hipp_mean_sq <- MRS_surface$hipp_mean_c^2
# cortical_thickness_adsignature_dickson
MRS_surface$cortical_thickness_adsignature_dickson_c  <- scale(MRS_surface$cortical_thickness_adsignature_dickson, center = TRUE, scale = FALSE)
MRS_surface$cortical_thickness_adsignature_dickson_sq <- MRS_surface$cortical_thickness_adsignature_dickson_c^2
# thick_normalized_surface
MRS_surface$thick_normalized_surface_c  <- scale(MRS_surface$thick_normalized_surface, center = TRUE, scale = FALSE)
MRS_surface$thick_normalized_surface_sq <- MRS_surface$thick_normalized_surface_c^2
#  activation_hipp mean
MRS_surface$hipp_mean_act_c  <- scale(MRS_surface$hipp_mean_act, center = TRUE, scale = FALSE)
MRS_surface$hipp_mean_act_sq <- MRS_surface$hipp_mean_act_c^2
#  activation_parietal_sup_l 
MRS_surface$activation_parietal_sup_l_c  <- scale(MRS_surface$activation_parietal_sup_l, center = TRUE, scale = FALSE)
MRS_surface$activation_parietal_sup_l_sq <- MRS_surface$activation_parietal_sup_l_c^2




##### Creer de banques de donees pour les analyses
#### Structure ####
MRS_S_Prec <- MRS_surface[, c("m_m_precuneus", "m_m_precuneus_c", "m_m_precuneus_sq",
                              "hipp_mean", "hipp_mean_c", "hipp_mean_sq",
                              "cortical_thickness_adsignature_dickson", "cortical_thickness_adsignature_dickson_c", "cortical_thickness_adsignature_dickson_sq",
                              "thick_normalized_surface", "thick_normalized_surface_c", "thick_normalized_surface_sq")] |> na.omit()

MRS_S_ACC <- MRS_surface[, c("m_m_acc", "m_m_acc_c", "m_m_acc_sq",
                             "hipp_mean", "hipp_mean_c", "hipp_mean_sq",
                             "cortical_thickness_adsignature_dickson", "cortical_thickness_adsignature_dickson_c", "cortical_thickness_adsignature_dickson_sq",
                             "thick_normalized_surface", "thick_normalized_surface_c", "thick_normalized_surface_sq")] |> na.omit()



# Memory 
MRS_M_Prec <- MRS_surface[, c("m_m_precuneus", "m_m_precuneus_c", "m_m_precuneus_sq",
                              "face_name_rappel_differe_spectro", "face_name_rappel_differe_spectro_c", "face_name_rappel_differe_spectro_sq",
                              "memoria_libre_correcte", "memoria_libre_correcte_c", "memoria_libre_correcte_sq")] |> na.omit()

MRS_M_ACC <- MRS_surface[, c("m_m_acc", "m_m_acc_c", "m_m_acc_sq",
                             "face_name_rappel_differe_spectro", "face_name_rappel_differe_spectro_c", "face_name_rappel_differe_spectro_sq",
                             "memoria_libre_correcte", "memoria_libre_correcte_c", "memoria_libre_correcte_sq")] |> na.omit()



# Model 1 m_m_ACC ~ hipp mean
model_lin_acc_hipp_mean <- lm(m_m_acc ~ hipp_mean_c, data = MRS_S_ACC)
summary(model_lin_acc_hipp_mean)
model_quad_acc_hipp_mean <- lm(m_m_acc ~ hipp_mean_c + hipp_mean_sq, data = MRS_S_ACC)
summary(model_quad_acc_hipp_mean)

AIC(model_lin_acc_hipp_mean)
AIC(model_quad_acc_hipp_mean)
anova(model_lin_acc_hipp_mean, model_quad_acc_hipp_mean)

# Model 2 m_m_Precuneus ~ hipp mean 
model_lin_prec_hipp_mean <- lm(m_m_precuneus ~ hipp_mean_c, data = MRS_S_Prec)
summary(model_lin_prec_hipp_mean)
model_quad_prec_hipp_mean <- lm(m_m_precuneus ~ hipp_mean_c + hipp_mean_sq, data = MRS_S_Prec)
summary(model_quad_prec_hipp_mean)

AIC(model_lin_prec_hipp_mean)
AIC(model_quad_prec_hipp_mean)
anova(model_lin_prec_hipp_mean, model_quad_prec_hipp_mean)

# Model 3 m_m_ACC ~ thickness
model_lin_acc_thick <- lm(m_m_acc ~ thick_normalized_surface_c, data = MRS_S_ACC)
summary(model_lin_acc_thick)
model_quad_acc_thick <- lm(m_m_acc ~ thick_normalized_surface_c + thick_normalized_surface_sq, data = MRS_S_ACC)
summary(model_quad_acc_thick)

AIC(model_lin_acc_thick)
AIC(model_quad_acc_thick)
anova(model_lin_acc_thick, model_quad_acc_thick)

# Model 4 m_m_Precuneus ~ thickness
model_lin_prec_thick <- lm(m_m_precuneus ~ thick_normalized_surface_c , data = MRS_S_Prec)
summary(model_lin_prec_thick)
model_quad_prec_thick <- lm(m_m_precuneus ~ thick_normalized_surface_c + thick_normalized_surface_sq, data = MRS_S_Prec)
summary(model_quad_prec_thick)

AIC(model_lin_prec_thick)
AIC(model_quad_prec_thick)
BIC(model_lin_prec_thick)
BIC(model_quad_prec_thick)
anova(model_lin_prec_thick, model_quad_prec_thick)

