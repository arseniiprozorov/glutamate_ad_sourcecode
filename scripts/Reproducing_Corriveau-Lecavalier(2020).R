
## reproduce Nick 2020

##### Creer de banques de donees pour les analyses

# Structure ####
Structure_thick_hipp <- MRS_full[, c("hipp_mean","hipp_mean_c","hipp_mean_sq","cortical_thickness_adsignature_dickson",
                                     "cortical_thickness_adsignature_dickson_c","cortical_thickness_adsignature_dickson_sq", "age_spectro_c")] |> na.omit()

# Activaiton 
Activation_parietal_hipp  <- MRS_full[, c("hipp_mean_act","hipp_mean_act_c","hipp_mean_act_sq","activation_parietal_sup_l", 
                                          "activation_parietal_sup_l_c","activation_parietal_sup_l_sq","age_spectro_c")] |> na.omit()

# Memory 
Memory_assos_memoria <- MRS_full[, c("memoria_libre_correcte", "memoria_libre_correcte_c", "memoria_libre_correcte_sq",
                           "associative_memory_performance","associative_memory_performance_c",
                           "associative_memory_performance_sq", "age_spectro_c")] |> na.omit()

# Combined
Structure_activation <- MRS_full[, c("hipp_mean","hipp_mean_c","hipp_mean_sq","cortical_thickness_adsignature_dickson",
                                     "cortical_thickness_adsignature_dickson_c","cortical_thickness_adsignature_dickson_sq", 
                                     "hipp_mean_act","hipp_mean_act_c","hipp_mean_act_sq","activation_parietal_sup_l", 
                                     "activation_parietal_sup_l_c","activation_parietal_sup_l_sq","age_spectro_c")] |> na.omit()


# Models (activation ~ structure)
model_lin_act_parietal_thick <- lm(activation_parietal_sup_l ~ cortical_thickness_adsignature_dickson_c, data= Structure_activation)
summary(model_lin_act_parietal_thick)
model_quad_act_parietal_thick <- lm(activation_parietal_sup_l ~ cortical_thickness_adsignature_dickson_c + cortical_thickness_adsignature_dickson_sq, data= Structure_activation)
summary(model_quad_act_parietal_thick)  

model_lin_act_parietal_hipp <- lm(activation_parietal_sup_l ~ hipp_mean_c, data= Structure_activation)
summary(model_lin_act_parietal_hipp)
model_quad_act_parietal_hipp <- lm(activation_parietal_sup_l ~ hipp_mean_c + hipp_mean_sq, data= Structure_activation)
summary(model_quad_act_parietal_hipp)  


model_lin_act_hipp_thick <- lm(hipp_mean_act ~ cortical_thickness_adsignature_dickson_c, data= Structure_activation)
summary(model_lin_act_hipp_thick)
model_quad_act_hipp_thick <- lm(hipp_mean_act ~ cortical_thickness_adsignature_dickson_c + cortical_thickness_adsignature_dickson_sq, data= Structure_activation)
summary(model_quad_act_hipp_thick)  


model_lin_act_hipp_hipp <- lm(hipp_mean_act ~ hipp_mean_c, data= Structure_activation)
summary(model_lin_act_hipp_hipp)
model_quad_act_hipp_hipp <- lm(hipp_mean_act ~ hipp_mean_c + hipp_mean_sq, data= Structure_activation)
summary(model_quad_act_hipp_hipp)  



# Combined MCI et SCD+
# 1. Create a logical filter to find MCI and SCD+ rows
# This creates a 'true/false' list for every row
rows_to_keep <- MRS_full$diagnostic_nick == "MCI" | MRS_full$diagnostic_nick == "SCD+"

# 2. Filter the rows and select the columns simultaneously
# Syntax: data[rows, columns]
Structure_activation_MCI_SCD <- MRS_full[rows_to_keep, c("hipp_mean", "hipp_mean_c", "hipp_mean_sq", 
                                                         "cortical_thickness_adsignature_dickson", 
                                                         "cortical_thickness_adsignature_dickson_c", 
                                                         "cortical_thickness_adsignature_dickson_sq", 
                                                         "hipp_mean_act", "hipp_mean_act_c", "hipp_mean_act_sq", 
                                                         "activation_parietal_sup_l", "activation_parietal_sup_l_c", 
                                                         "activation_parietal_sup_l_sq", "age_spectro_c")]

# 3. Remove any rows that have missing data (NA)
Structure_activation_MCI_SCD <- na.omit(Structure_activation_MCI_SCD)


# 4. Same models without HC

model_lin_act_parietal_thick <- lm(activation_parietal_sup_l ~ cortical_thickness_adsignature_dickson_c, data= Structure_activation_MCI_SCD)
summary(model_lin_act_parietal_thick)
model_quad_act_parietal_thick <- lm(activation_parietal_sup_l ~ cortical_thickness_adsignature_dickson_c + cortical_thickness_adsignature_dickson_sq, data= Structure_activation_MCI_SCD)
summary(model_quad_act_parietal_thick)  
model_lin_act_parietal_hipp <- lm(activation_parietal_sup_l ~ hipp_mean_c, data= Structure_activation_MCI_SCD)
summary(model_lin_act_parietal_hipp)
model_quad_act_parietal_hipp <- lm(activation_parietal_sup_l ~ hipp_mean_c + hipp_mean_sq, data= Structure_activation_MCI_SCD)
summary(model_quad_act_parietal_hipp)  


model_lin_act_hipp_thick <- lm(hipp_mean_act ~ cortical_thickness_adsignature_dickson_c, data= Structure_activation_MCI_SCD)
summary(model_lin_act_hipp_thick)
model_quad_act_hipp_thick <- lm(hipp_mean_act ~ cortical_thickness_adsignature_dickson_c + cortical_thickness_adsignature_dickson_sq, data= Structure_activation_MCI_SCD)
summary(model_quad_act_hipp_thick)  

model_lin_act_hipp_hipp <- lm(hipp_mean_act ~ hipp_mean_c, data= Structure_activation_MCI_SCD)
summary(model_lin_act_hipp_hipp)
model_quad_act_hipp_hipp <- lm(hipp_mean_act ~ hipp_mean_c + hipp_mean_sq, data= Structure_activation_MCI_SCD)
summary(model_quad_act_hipp_hipp)  




################ Reproducing memory results

# Activaiton 
Activation_parietal_hipp  <- MRS_full[, c("hipp_mean_act","hipp_mean_act_c","hipp_mean_act_sq","activation_parietal_sup_l", 
                                          "activation_parietal_sup_l_c","activation_parietal_sup_l_sq","age_spectro_c")] |> na.omit()

# Memory 
MRS_M_Prec <- MRS_full[, c("memoria_libre_correcte", "memoria_libre_correcte_c", "memoria_libre_correcte_sq",
                           "associative_memory_performance","associative_memory_performance_c",
                           "associative_memory_performance_sq", "age_spectro_c")] |> na.omit()

# Combined
Memory_activation <- MRS_full[, c("memoria_libre_correcte", "memoria_libre_correcte_c", "memoria_libre_correcte_sq",
                                  "associative_memory_performance","associative_memory_performance_c",
                                  "associative_memory_performance_sq","hipp_mean_act","hipp_mean_act_c","hipp_mean_act_sq","activation_parietal_sup_l", 
                                  "activation_parietal_sup_l_c","activation_parietal_sup_l_sq","age_spectro_c")] |> na.omit()




# Models 
model_lin_act_parietal_memoria <- lm(activation_parietal_sup_l ~ memoria_libre_correcte_c, data= Memory_activation)
summary(model_lin_act_parietal_memoria)
model_quad_act_parietal_memoria <- lm(activation_parietal_sup_l ~ memoria_libre_correcte_c + memoria_libre_correcte_sq, data= Memory_activation)
summary(model_quad_act_parietal_memoria) 

model_lin_act_parietal_assos <- lm(activation_parietal_sup_l ~ associative_memory_performance_c, data= Memory_activation)
summary(model_lin_act_parietal_assos)
model_quad_act_parietal_assos <- lm(activation_parietal_sup_l ~ associative_memory_performance_c + associative_memory_performance_sq, data= Memory_activation)
summary(model_quad_act_parietal_assos)  

model_lin_act_hipp_assos <- lm(hipp_mean_act ~ associative_memory_performance_c, data= Memory_activation)
summary(model_lin_act_hipp_assos)
model_quad_act_hipp_assos<- lm(hipp_mean_act ~ associative_memory_performance_c + associative_memory_performance_sq, data= Memory_activation)
summary(model_quad_act_hipp_assos)  









