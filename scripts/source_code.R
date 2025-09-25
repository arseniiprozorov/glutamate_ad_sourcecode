library(afex)
library(readxl)  
library(haven)
library(jmv)
library(sjmisc) 
library(sjstats) 
library(labelled)
library(olsrr)
library(methods)
library(lme4)
library(ggplot2)
library(emmeans)
library(tseries)
library(janitor)
library(dplyr)
library(interactions)
library(tidyr)
library(ggpubr)



## Glutamate Relates to Structural and Functional markers of Disease Severity in early Alzheimer’s disease  ######

#ANALYSES PRÉLIMINAIRES :
#Création d’une banque de données

import_data_arsenii_20250609 <- read_excel("supplementary/import_data_arsenii_20250609.xlsx")
MRS_full <- import_data_arsenii_20250609

# Clean the column names
MRS_full <- janitor::clean_names(MRS_full)

# Convertir  en numérique
MRS_full$education <- as.numeric(MRS_full$education)
MRS_full$moca_corr_spectro <- as.numeric(MRS_full$moca_corr_spectro)
MRS_full$dprime_hit_fa <- as.numeric(MRS_full$dprime_hit_fa)
MRS_full$associative_memory_performance <- as.numeric(MRS_full$associative_memory_performance)
MRS_full$moca_corr_spectro <- as.numeric(MRS_full$moca_corr_spectro)
# Regrouper les groupes 'HC' et 'SCD' dans un seul groupe 'SCD'
MRS_full$diagnostic_nick <- ifelse(MRS_full$diagnostic_nick %in% c("HC", "SCD"), "HC", MRS_full$diagnostic_nick)

# Convertir diagnostic_nick en facteur avec l’ordre HC, SCD+, MCI
MRS_full$diagnostic_nick <- factor(MRS_full$diagnostic_nick, levels = c("HC", "SCD+", "MCI"))

# Convertir sex en facteur
MRS_full$sex <- as.factor(MRS_full$sex)


#Vérifier les varibales
lapply(MRS_full,class)   
levels(MRS_full$diagnostic_nick)
table(MRS_full$diagnostic_nick)


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
MRS_full$memoria_libre_correcte <- winsorize_iqr(MRS_full$memoria_libre_correcte, 1.5)
MRS_full$face_name_rappel_differe_spectro <- winsorize_iqr(MRS_full$face_name_rappel_differe_spectro, 1.5)
MRS_full$associative_memory_performance <- winsorize_iqr(MRS_full$associative_memory_performance, 1.5)
MRS_full$activation_temporal_inf_r <- winsorize_iqr(MRS_full$activation_temporal_inf_r, 1.5)
MRS_full$activation_parietal_sup_l <- winsorize_iqr(MRS_full$activation_parietal_sup_l, 1.5)
MRS_full$activation_hippocampus_l <- winsorize_iqr(MRS_full$activation_hippocampus_l, 1.5)
MRS_full$moca_corr_spectro <- winsorize_iqr(MRS_full$moca_corr_spectro, 1.5)
MRS_full$dprime_hit_fa <- winsorize_iqr(MRS_full$dprime_hit_fa, 1.5)
MRS_full$associative_memory_performance_sylvie <- winsorize_iqr(MRS_full$associative_memory_performance_sylvie, 1.5)


# Apply winsorization for 1.7 IQR variables
MRS_full$hip_l_nor_icv <- winsorize_iqr(MRS_full$hip_l_nor_icv, 1.7)
MRS_full$cortical_thickness_adsignature_dickson <- winsorize_iqr(MRS_full$cortical_thickness_adsignature_dickson, 1.7)
MRS_full$m_m_acc <- winsorize_iqr(MRS_full$m_m_acc, 1.7)
MRS_full$m_m_precuneus <- winsorize_iqr(MRS_full$m_m_precuneus, 1.7)


#jmv::descriptives(data = MRS_full, vars = vars(),sd = TRUE,range = TRUE,skew = TRUE,kurt = TRUE)



#### Centrer et mettre au carré
# education
MRS_full$education_c  <- scale(MRS_full$education, center = TRUE, scale = FALSE)
MRS_full$education_sq <- MRS_full$education_c^2
# age_spectro
MRS_full$age_spectro_c  <- scale(MRS_full$age_spectro, center = TRUE, scale = FALSE)
MRS_full$age_spectro_sq <- MRS_full$age_spectro_c^2
# m_m_precuneus
MRS_full$m_m_precuneus_c  <- scale(MRS_full$m_m_precuneus, center = TRUE, scale = FALSE)
MRS_full$m_m_precuneus_sq <- MRS_full$m_m_precuneus_c^2
# m_m_acc
MRS_full$m_m_acc_c  <- scale(MRS_full$m_m_acc, center = TRUE, scale = FALSE)
MRS_full$m_m_acc_sq <- MRS_full$m_m_acc_c^2
# moca_corr_spectro
MRS_full$moca_corr_spectro_c  <- scale(MRS_full$moca_corr_spectro, center = TRUE, scale = FALSE)
MRS_full$moca_corr_spectro_sq <- MRS_full$moca_corr_spectro_c^2
# memoria_libre_correcte
MRS_full$memoria_libre_correcte_c  <- scale(MRS_full$memoria_libre_correcte, center = TRUE, scale = FALSE)
MRS_full$memoria_libre_correcte_sq <- MRS_full$memoria_libre_correcte_c^2
# face_name_rappel_differe_spectro
MRS_full$face_name_rappel_differe_spectro_c  <- scale(MRS_full$face_name_rappel_differe_spectro, center = TRUE, scale = FALSE)
MRS_full$face_name_rappel_differe_spectro_sq <- MRS_full$face_name_rappel_differe_spectro_c^2
# hip_l_nor_icv
MRS_full$hip_l_nor_icv_c  <- scale(MRS_full$hip_l_nor_icv, center = TRUE, scale = FALSE)
MRS_full$hip_l_nor_icv_sq <- MRS_full$hip_l_nor_icv_c^2
# cortical_thickness_adsignature_dickson
MRS_full$cortical_thickness_adsignature_dickson_c  <- scale(MRS_full$cortical_thickness_adsignature_dickson, center = TRUE, scale = FALSE)
MRS_full$cortical_thickness_adsignature_dickson_sq <- MRS_full$cortical_thickness_adsignature_dickson_c^2
#  dprime_hit_fa
MRS_full$dprime_hit_fa_c  <- scale(MRS_full$dprime_hit_fa, center = TRUE, scale = FALSE)
MRS_full$dprime_hit_fa_sq <- MRS_full$dprime_hit_fa_c^2
#  associative_memory_performance
MRS_full$associative_memory_performance_c  <- scale(MRS_full$associative_memory_performance, center = TRUE, scale = FALSE)
MRS_full$associative_memory_performance_sq <- MRS_full$associative_memory_performance_c^2
#  associative_memory_performance_sylvie
MRS_full$associative_memory_performance_sylvie_c  <- scale(MRS_full$associative_memory_performance_sylvie, center = TRUE, scale = FALSE)
MRS_full$associative_memory_performance_sylvie_sq <- MRS_full$associative_memory_performance_sylvie_c^2
#  activation_hippocampus_l
MRS_full$activation_hippocampus_l_c  <- scale(MRS_full$activation_hippocampus_l, center = TRUE, scale = FALSE)
MRS_full$activation_hippocampus_l_sq <- MRS_full$activation_hippocampus_l_c^2
#  activation_parietal_sup_l 
MRS_full$activation_parietal_sup_l_c  <- scale(MRS_full$activation_parietal_sup_l, center = TRUE, scale = FALSE)
MRS_full$activation_parietal_sup_l_sq <- MRS_full$activation_parietal_sup_l_c^2
#  activation_temporal_inf_r
MRS_full$activation_temporal_inf_r_c  <- scale(MRS_full$activation_temporal_inf_r, center = TRUE, scale = FALSE)
MRS_full$activation_temporal_inf_r_sq <- MRS_full$activation_temporal_inf_r_c^2





################## Table 1 participants charactertistics  ################
# Sex
# chi-square
participant_sex <- table(MRS_full$diagnostic_nick, MRS_full$sex)
participant_sex
chisq.test(participant_sex)
fisher.test(participant_sex[c("HC","SCD+"), ])
fisher.test(participant_sex[c("HC","MCI"), ])
fisher.test(participant_sex[c("SCD+","MCI"), ]) 

# Age
anova_age <- aov(age_spectro ~ diagnostic_nick, data = MRS_full)
summary(anova_age)               
TukeyHSD(anova_age) 

jmv::descriptives(data = MRS_full,vars = vars(age_spectro), splitBy = "diagnostic_nick",sd = TRUE,hist = TRUE)

# Education
anova_edu <- aov(education ~ diagnostic_nick, data = MRS_full)
summary(anova_edu)
TukeyHSD(anova_edu)

jmv::descriptives(data = MRS_full,vars = vars(education), splitBy = "diagnostic_nick",sd = TRUE, hist = TRUE)

# MoCA
anova_moca <- aov(moca_corr_spectro ~ diagnostic_nick, data = MRS_full)
summary(anova_moca)
TukeyHSD(anova_moca)

jmv::descriptives(data = MRS_full,vars = vars(moca_corr_spectro), splitBy = "diagnostic_nick",sd = TRUE, hist = TRUE)

# Memoria Free Word recall
anova_free <- aov(memoria_libre_correcte ~ diagnostic_nick, data = MRS_full)
summary(anova_free)
TukeyHSD(anova_free)

jmv::descriptives(data = MRS_full,vars = vars(memoria_libre_correcte), splitBy = "diagnostic_nick",sd = TRUE, hist = TRUE)

# Face-name recall
anova_face <- aov(face_name_rappel_differe_spectro ~ diagnostic_nick, data = MRS_full)
summary(anova_face)
TukeyHSD(anova_face)

jmv::descriptives(data = MRS_full,vars = vars(face_name_rappel_differe_spectro), splitBy = "diagnostic_nick",sd = TRUE, hist = TRUE)


# d prime
anova_dprime <- aov(dprime_hit_fa ~ diagnostic_nick, data = MRS_full)
summary(anova_dprime)
TukeyHSD(anova_dprime)

jmv::descriptives(data = MRS_full,vars = vars(dprime_hit_fa), splitBy = "diagnostic_nick",sd = TRUE, hist = TRUE)


# Associative memory
anova_assoc <- aov(associative_memory_performance ~ diagnostic_nick, data = MRS_full)
summary(anova_assoc)
TukeyHSD(anova_assoc)

jmv::descriptives(data = MRS_full,vars = vars(associative_memory_performance), splitBy = "diagnostic_nick",sd = TRUE, hist = TRUE)



################################################### SEM ####################################################################
###### moderation ######
names(MRS_full)
DATA_REG0 <- DATA_REG[, c(
  "m_m_precuneus", "m_m_precuneus_c", "m_m_precuneus_sq",
  "memoria_libre_correcte", "memoria_libre_correcte_c", "memoria_libre_correcte_sq",
  "hip_l_nor_icv", "hip_l_nor_icv_c", "hip_l_nor_icv_sq",
  "cortical_thickness_adsignature_dickson", "cortical_thickness_adsignature_dickson_c", "cortical_thickness_adsignature_dickson_sq"
)] |> na.omit()


### Proxies of disease severity ######

model_quad_proxy_memory <- lm(m_m_precuneus ~ memoria_libre_correcte_c + memoria_libre_correcte_sq, data = DATA_REG0)
summary(model_quad_proxy_memory)
model_lin_proxy_memory <- lm(m_m_precuneus ~ memoria_libre_correcte_c, data = DATA_REG0)
summary(model_lin_proxy_memory)


model_quad_proxy_memory <- lm(m_m_precuneus ~ memoria_libre_correcte_c + memoria_libre_correcte_sq, data = DATA_REG0)
summary(model_quad_proxy_memory)
model_lin_proxy_memory <- lm(m_m_precuneus ~ memoria_libre_correcte_c, data = DATA_REG0)
summary(model_lin_proxy_memory)



### mediation leg 1 ###
## S~G + G^2
#hip
model_quad_hip <- lm(hip_l_nor_icv ~ m_m_precuneus_c + m_m_precuneus_sq, data = DATA_REG0)
summary(model_quad_hip)
model_lin_hip <- lm(hip_l_nor_icv ~ m_m_precuneus_c, data = DATA_REG0)
summary(model_lin_hip)

AIC(model_quad_hip)   
AIC(model_lin_hip)   
#thick
model_quad_thick <- lm(cortical_thickness_adsignature_dickson ~ m_m_precuneus_c + m_m_precuneus_sq, data = DATA_REG0)
summary(model_quad_hip)
model_lin_thick <- lm(cortical_thickness_adsignature_dickson ~ m_m_precuneus_c, data = DATA_REG0)
summary(model_lin_thick)

AIC(model_quad_thick)   
AIC(model_lin_thick)   


### mediation leg 2 ###

## M~S + G^2
model_quad_hip2 <- lm(memoria_libre_correcte ~ hip_l_nor_icv_c  + hip_l_nor_icv_sq, data = DATA_REG0)
summary(model_quad_hip2)
model_lin_hip2 <- lm(memoria_libre_correcte ~ hip_l_nor_icv_c, data = DATA_REG0)
summary(model_lin_hip2)

AIC(model_quad_hip2)   
AIC(model_lin_hip2)   


## direct path G -->M
## M~G + G^2
model_quad_MG <- lm(memoria_libre_correcte ~ m_m_precuneus_c + m_m_precuneus_sq, data = DATA_REG0)
summary(model_quad_MG)
model_lin_MG <- lm(memoria_libre_correcte ~ m_m_precuneus_c, data = DATA_REG0)
summary(model_lin_MG)

AIC(model_quad_MG)   
AIC(model_lin_MG)   


### moderation ####
model_mod_hipp <- lm(memoria_libre_correcte ~ m_m_precuneus_c + hip_l_nor_icv_c:m_m_precuneus_c, data = DATA_REG0)
summary(model_mod_hipp)

model_mod_thick <- lm(memoria_libre_correcte ~ m_m_precuneus_c + cortical_thickness_adsignature_dickson_c:m_m_precuneus_c, data = DATA_REG0)
summary(model_mod_thick)
### nested model ####
model_nested_hipp <- lm(memoria_libre_correcte ~ m_m_precuneus_c + m_m_precuneus_sq + hip_l_nor_icv_c:m_m_precuneus_sq, data = DATA_REG0)
summary(model_nested_hipp)

model_nested_thick <- lm(memoria_libre_correcte ~ m_m_precuneus_c + m_m_precuneus_sq + cortical_thickness_adsignature_dickson_c:m_m_precuneus_sq, data = DATA_REG0)
summary(model_nested_thick)

names(DATA_REG0)


## mediation
library(mediation)

# Pick the centered versions consistently
DATA1 <- DATA_REG0[, c("memoria_libre_correcte",
                       "m_m_precuneus_c",
                       "hip_l_nor_icv_c")]
DATA1 <- na.omit(DATA1)

# Make sure they are plain numeric
DATA1$m_m_precuneus_c  <- as.numeric(DATA1$m_m_precuneus_c)
DATA1$hip_l_nor_icv_c  <- as.numeric(DATA1$hip_l_nor_icv_c)

# 1. Mediator model: S ~ G
med.fit <- lm(hip_l_nor_icv_c ~ m_m_precuneus_c, data = DATA1)

# 2. Outcome model: M ~ G + S
out.fit <- lm(memoria_libre_correcte ~ m_m_precuneus_c + hip_l_nor_icv_c, data = DATA1)

# 3. Mediation analysis
set.seed(123)
med.out <- mediate(med.fit, out.fit,
                   treat = "m_m_precuneus_c",
                   mediator = "hip_l_nor_icv_c",
                   boot = TRUE, sims = 5000)

summary(med.out)




### mediation allowing for quadratic terms ####
DATA1$m_m_precuneus_sq <- DATA1$m_m_precuneus_c^2

# Mediator model: S ~ G (+ G^2)
med.fit <- lm(hip_l_nor_icv_c ~ m_m_precuneus_c + m_m_precuneus_sq, data = DATA1)

# Outcome model: M ~ G + G^2 + S
out.fit <- lm(memoria_libre_correcte ~ m_m_precuneus_c + m_m_precuneus_sq + hip_l_nor_icv_c, data = DATA1)

library(mediation)
set.seed(123)
med.out <- mediate(med.fit, out.fit,
                   treat = "m_m_precuneus_c",
                   mediator = "hip_l_nor_icv_c",
                   boot = TRUE, sims = 5000)
summary(med.out)



###################### AIC quadratic > 5.6 than linear ####################
## mediation (linear)
library(mediation)

# Pick the centered versions consistently
DATA1 <- DATA_REG0[, c("memoria_libre_correcte",
                       "m_m_precuneus_c",
                       "hip_l_nor_icv_c")]
DATA1 <- na.omit(DATA1)

# Make sure they are plain numeric
DATA1$m_m_precuneus_c  <- as.numeric(DATA1$m_m_precuneus_c)
DATA1$hip_l_nor_icv_c  <- as.numeric(DATA1$hip_l_nor_icv_c)

# 1. Mediator model: S ~ G
med.fit.lin <- lm(hip_l_nor_icv_c ~ m_m_precuneus_c, data = DATA1)

# 2. Outcome model: M ~ G + S
out.fit.lin <- lm(memoria_libre_correcte ~ m_m_precuneus_c + hip_l_nor_icv_c, data = DATA1)

# AICs
cat("\nAIC mediator (linear): ", AIC(med.fit.lin), "\n")
cat("AIC outcome  (linear): ", AIC(out.fit.lin), "\n")

# 3. Mediation analysis
set.seed(123)
med.out.lin <- mediate(med.fit.lin, out.fit.lin,
                       treat = "m_m_precuneus_c",
                       mediator = "hip_l_nor_icv_c",
                       boot = TRUE, sims = 5000)
summary(med.out.lin)



### mediation allowing for quadratic terms ####
DATA1$m_m_precuneus_sq <- DATA1$m_m_precuneus_c^2

# Mediator model: S ~ G (+ G^2)
med.fit.quad <- lm(hip_l_nor_icv_c ~ m_m_precuneus_c + m_m_precuneus_sq, data = DATA1)

# Outcome model: M ~ G + G^2 + S
out.fit.quad <- lm(memoria_libre_correcte ~ m_m_precuneus_c + m_m_precuneus_sq + hip_l_nor_icv_c, data = DATA1)

# AICs
cat("\nAIC mediator (quad): ", AIC(med.fit.quad), "\n")
cat("AIC outcome  (quad): ", AIC(out.fit.quad), "\n")

# 3. Mediation analysis
set.seed(123)
med.out.quad <- mediate(med.fit.quad, out.fit.quad,
                        treat = "m_m_precuneus_c",
                        mediator = "hip_l_nor_icv_c",
                        boot = TRUE, sims = 5000)
summary(med.out.quad)





### ACC ######
# Create DATA_REG2 for the ACC analyses
DATA_REG2 <- MRS_full[, c(
  "m_m_acc", "memoria_libre_correcte", "memoria_libre_correcte_c", "memoria_libre_correcte_sq",
  "hip_l_nor_icv", "hip_l_nor_icv_c", "hip_l_nor_icv_sq",
  "cortical_thickness_adsignature_dickson", "cortical_thickness_adsignature_dickson_c", "cortical_thickness_adsignature_dickson_sq"
)] |> na.omit()

# ACC centered & squared (ensure plain numeric)
DATA_REG2$m_m_acc_c  <- as.numeric(scale(DATA_REG2$m_m_acc, center = TRUE, scale = FALSE))
DATA_REG2$m_m_acc_sq <- DATA_REG2$m_m_acc_c^2

### mediation leg 1 ###
## S ~ G + G^2
# hip
model_quad_hip <- lm(hip_l_nor_icv ~ m_m_acc_c + m_m_acc_sq, data = DATA_REG2)
summary(model_quad_hip)
model_lin_hip  <- lm(hip_l_nor_icv ~ m_m_acc_c, data = DATA_REG2)
summary(model_lin_hip)

AIC(model_quad_hip)
AIC(model_lin_hip)

# thick
model_quad_thick <- lm(cortical_thickness_adsignature_dickson ~ m_m_acc_c + m_m_acc_sq, data = DATA_REG2)
summary(model_quad_thick)
model_lin_thick  <- lm(cortical_thickness_adsignature_dickson ~ m_m_acc_c, data = DATA_REG2)
summary(model_lin_thick)

AIC(model_quad_thick)
AIC(model_lin_thick)


### mediation leg 2 ###
## M ~ S (+ S^2 if desired; here you used hip^2)
model_quad_hip2 <- lm(memoria_libre_correcte ~ hip_l_nor_icv_c + hip_l_nor_icv_sq, data = DATA_REG2)
summary(model_quad_hip2)
model_lin_hip2  <- lm(memoria_libre_correcte ~ hip_l_nor_icv_c, data = DATA_REG2)
summary(model_lin_hip2)

AIC(model_quad_hip2)
AIC(model_lin_hip2)


## direct path G --> M
## M ~ G + G^2
model_quad_MG <- lm(memoria_libre_correcte ~ m_m_acc_c + m_m_acc_sq, data = DATA_REG2)
summary(model_quad_MG)
model_lin_MG  <- lm(memoria_libre_correcte ~ m_m_acc_c, data = DATA_REG2)
summary(model_lin_MG)

AIC(model_quad_MG)
AIC(model_lin_MG)


### moderation ####
model_mod_hipp  <- lm(memoria_libre_correcte ~ m_m_acc_c + hip_l_nor_icv_c:m_m_acc_c, data = DATA_REG2)
summary(model_mod_hipp)

model_mod_thick <- lm(memoria_libre_correcte ~ m_m_acc_c + cortical_thickness_adsignature_dickson_c:m_m_acc_c, data = DATA_REG2)
summary(model_mod_thick)

### nested model ####
model_nested_hipp  <- lm(memoria_libre_correcte ~ m_m_acc_c + m_m_acc_sq + hip_l_nor_icv_c:m_m_acc_sq, data = DATA_REG2)
summary(model_nested_hipp)

model_nested_thick <- lm(memoria_libre_correcte ~ m_m_acc_c + m_m_acc_sq + cortical_thickness_adsignature_dickson_c:m_m_acc_sq, data = DATA_REG2)
summary(model_nested_thick)

names(DATA_REG2)






### graph ###
# --- Ajuster le modèle (tel que tu l'as fait) ---
model_nested_thick <- lm(
  memoria_libre_correcte ~ m_m_acc_c + m_m_acc_sq +
    cortical_thickness_adsignature_dickson_c:m_m_acc_sq,
  data = DATA_REG2
)
summary(model_nested_thick)

# --- Préparer une grille de valeurs de G centré (ACC glutamate) ---
rng_G  <- range(DATA_REG2$m_m_acc_c, na.rm = TRUE)
grid_G <- seq(rng_G[1], rng_G[2], length.out = 200)
grid_G2 <- grid_G^2

# --- Définir "faible" et "élevée" épaisseur (±1 SD) ---
S_mean <- mean(DATA_REG2$cortical_thickness_adsignature_dickson_c, na.rm = TRUE)
S_sd   <- sd(DATA_REG2$cortical_thickness_adsignature_dickson_c, na.rm = TRUE)
S_low  <- S_mean - S_sd
S_high <- S_mean + S_sd

# --- Construire jeux de données pour la prédiction ---
new_low <- data.frame(
  m_m_acc_c = grid_G,
  m_m_acc_sq = grid_G2,
  cortical_thickness_adsignature_dickson_c = S_low
)

new_high <- data.frame(
  m_m_acc_c = grid_G,
  m_m_acc_sq = grid_G2,
  cortical_thickness_adsignature_dickson_c = S_high
)

# --- Prédictions et IC 95% pour chaque niveau de S ---
pred_low  <- predict(model_nested_thick, newdata = new_low,  se.fit = TRUE)
pred_high <- predict(model_nested_thick, newdata = new_high, se.fit = TRUE)

alpha <- 0.05
crit  <- qnorm(1 - alpha/2)

df_plot <- rbind(
  data.frame(
    S_level = paste0("Faible épaisseur (", round(S_low, 2), ")"),
    G_c = grid_G,
    M_hat = pred_low$fit,
    lo = pred_low$fit - crit * pred_low$se.fit,
    hi = pred_low$fit + crit * pred_low$se.fit
  ),
  data.frame(
    S_level = paste0("Épaisseur élevée (", round(S_high, 2), ")"),
    G_c = grid_G,
    M_hat = pred_high$fit,
    lo = pred_high$fit - crit * pred_high$se.fit,
    hi = pred_high$fit + crit * pred_high$se.fit
  )
)

# --- Graphique (ggplot2) ---
# install.packages("ggplot2") # si nécessaire
library(ggplot2)

ggplot(df_plot, aes(x = G_c, y = M_hat, group = S_level)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = S_level), alpha = 0.15, linewidth = 0) +
  geom_line(aes(linetype = S_level)) +
  labs(
    x = "Glutamate ACC (centré, m_m_acc_c)",
    y = "Mémoire (prédit) : memoria_libre_correcte",
    title = "Interaction S × G² : forme de la courbe selon l’épaisseur corticale",
    subtitle = "Modèle : M ~ G + G² + S:G² (sans effet principal de S)"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.title = element_blank())

# --- Option alternative : split médian (si tu préfères) ---
# seuil <- median(DATA_REG2$cortical_thickness_adsignature_dickson_c, na.rm = TRUE)
# S_low  <- seuil - 1e-8
# S_high <- seuil + 1e-8
# (puis réutiliser le même code de prédiction en remplaçant S_low / S_high)






# 1) Coerce S to a numeric vector and refit
DATA_REG2$S_vec <- as.numeric(DATA_REG2$cortical_thickness_adsignature_dickson_c)

model_nested_thick <- lm(
  memoria_libre_correcte ~ m_m_acc_c + m_m_acc_sq + S_vec:m_m_acc_sq,
  data = DATA_REG2
)
summary(model_nested_thick)

# 2) Grid of G (centered ACC glutamate)
rng_G  <- range(DATA_REG2$m_m_acc_c, na.rm = TRUE)
grid_G <- seq(rng_G[1], rng_G[2], length.out = 200)

# 3) Define low/high S as ±1 SD (or use quantiles if you prefer)
S_mean <- mean(DATA_REG2$S_vec, na.rm = TRUE)
S_sd   <- sd(DATA_REG2$S_vec, na.rm = TRUE)
S_low  <- S_mean - S_sd
S_high <- S_mean + S_sd

# 4) Build newdata for predictions
new_low <- data.frame(
  m_m_acc_c = grid_G,
  m_m_acc_sq = grid_G^2,
  S_vec = S_low
)
new_high <- data.frame(
  m_m_acc_c = grid_G,
  m_m_acc_sq = grid_G^2,
  S_vec = S_high
)

# 5) Predictions + 95% CI
pred_low  <- predict(model_nested_thick, newdata = new_low,  se.fit = TRUE)
pred_high <- predict(model_nested_thick, newdata = new_high, se.fit = TRUE)
crit <- qnorm(0.975)

df_plot <- rbind(
  data.frame(S_level = paste0("Faible épaisseur (", round(S_low,2), ")"),
             G_c = grid_G, M_hat = pred_low$fit,
             lo = pred_low$fit - crit*pred_low$se.fit,
             hi = pred_low$fit + crit*pred_low$se.fit),
  data.frame(S_level = paste0("Épaisseur élevée (", round(S_high,2), ")"),
             G_c = grid_G, M_hat = pred_high$fit,
             lo = pred_high$fit - crit*pred_high$se.fit,
             hi = pred_high$fit + crit*pred_high$se.fit)
)

library(ggplot2)
ggplot(df_plot, aes(x = G_c, y = M_hat, group = S_level)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = S_level), alpha = 0.15, linewidth = 0) +
  geom_line(aes(linetype = S_level)) +
  labs(x = "Glutamate ACC (centré, m_m_acc_c)",
       y = "Mémoire (prédite)",
       title = "Interaction S × G² : courbes pour faible vs élevée épaisseur") +
  theme_minimal(base_size = 12) +
  theme(legend.title = element_blank())





DATA_REG0 <- DATA_REG[, c(
  "m_m_precuneus", "m_m_precuneus_c", "m_m_precuneus_sq",
  "memoria_libre_correcte", "memoria_libre_correcte_c", "memoria_libre_correcte_sq",
  "hip_l_nor_icv", "hip_l_nor_icv_c", "hip_l_nor_icv_sq",
  "cortical_thickness_adsignature_dickson", "cortical_thickness_adsignature_dickson_c", "cortical_thickness_adsignature_dickson_sq"
)] |> na.omit()
##### linear precuneus ######
# --- 1) Coerce to numeric vectors and refit on the SAME data frame you’ll use ---
DATA_REG0$mmprec_c <- as.numeric(DATA_REG0$m_m_precuneus_c)   # predictor (G)
DATA_REG0$hip_vec  <- as.numeric(DATA_REG0$hip_l_nor_icv_c)   # moderator (H)

model_mod_precuneus <- lm(
  memoria_libre_correcte ~ mmprec_c + hip_vec:mmprec_c,
  data = DATA_REG0
)
summary(model_mod_precuneus)

# --- 2) Grid for G (no sign flip here; we’ll reverse the axis in ggplot) ---
rng_G  <- range(DATA_REG0$mmprec_c, na.rm = TRUE)
grid_G <- seq(rng_G[1], rng_G[2], length.out = 200)

# --- 3) Define low/high H as ±1 SD ---
H_mean <- mean(DATA_REG0$hip_vec, na.rm = TRUE)
H_sd   <- sd(DATA_REG0$hip_vec, na.rm = TRUE)
H_low  <- H_mean - H_sd
H_high <- H_mean + H_sd

# --- 4) Build newdata and predict with SEs ---
new_low  <- data.frame(mmprec_c = grid_G, hip_vec = H_low)
new_high <- data.frame(mmprec_c = grid_G, hip_vec = H_high)

pred_low  <- predict(model_mod_precuneus, newdata = new_low,  se.fit = TRUE)
pred_high <- predict(model_mod_precuneus, newdata = new_high, se.fit = TRUE)

crit <- qnorm(0.975)
df_plot <- rbind(
  data.frame(H_level = paste0("Hippocampe faible (", round(H_low,2), ")"),
             G_c = grid_G, M_hat = pred_low$fit,
             lo = pred_low$fit - crit*pred_low$se.fit,
             hi = pred_low$fit + crit*pred_low$se.fit),
  data.frame(H_level = paste0("Hippocampe élevé (", round(H_high,2), ")"),
             G_c = grid_G, M_hat = pred_high$fit,
             lo = pred_high$fit - crit*pred_high$se.fit,
             hi = pred_high$fit + crit*pred_high$se.fit)
)

# --- 5) Plot (visual x-axis inversion only) ---
library(ggplot2)
ggplot(df_plot, aes(x = G_c, y = M_hat, group = H_level)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = H_level), alpha = 0.15, linewidth = 0) +
  geom_line(aes(linetype = H_level)) +
  scale_x_reverse() +
  labs(x = "Glutamate précunéus (centré, axe inversé)",
       y = "Mémoire prédite",
       title = "Interaction linéaire G × Hippocampe (précunéus)") +
  theme_minimal(base_size = 12) +
  theme(legend.title = element_blank())
