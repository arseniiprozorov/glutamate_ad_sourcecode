library(here)
source(file.path("scripts", "source_code.R"))
library(car)

### Normality ###





##### Participants characteristics postulats ######

## Age
shapiro.test(residuals(anova_age))
ks.test(scale(residuals(anova_age)), "pnorm")
leveneTest(model.frame(anova_age)[[1]] ~ model.frame(anova_age)[[2]])

## Education
shapiro.test(residuals(anova_edu))
ks.test(scale(residuals(anova_edu)), "pnorm")
leveneTest(model.frame(anova_edu)[[1]] ~ model.frame(anova_edu)[[2]])

## MoCA
shapiro.test(residuals(anova_moca))
ks.test(scale(residuals(anova_moca)), "pnorm")
leveneTest(model.frame(anova_moca)[[1]] ~ model.frame(anova_moca)[[2]])

## Memoria Free Word recall
shapiro.test(residuals(anova_free))
ks.test(scale(residuals(anova_free)), "pnorm")
leveneTest(model.frame(anova_free)[[1]] ~ model.frame(anova_free)[[2]])





#### ANOVA postulats  #####
# ACC
shapiro.test(residuals(anova_acc))
ks.test(scale(residuals(anova_acc)), "pnorm")
leveneTest(model.frame(anova_acc)[[1]] ~ model.frame(anova_acc)[[2]])

# Precuneus
shapiro.test(residuals(anova_precuneus))
ks.test(scale(residuals(anova_precuneus)), "pnorm")
leveneTest(model.frame(anova_precuneus)[[1]] ~ model.frame(anova_precuneus)[[2]])




### Regression postulats ###############
models <- c(
  "model_lin_acc_hipp","model_lin_acc_hipp_act","model_lin_acc_memor","model_lin_acc_sup_act","model_lin_acc_thick",
  "model_lin_prec_hipp","model_lin_prec_hipp_act","model_lin_prec_memor","model_lin_prec_sup_act","model_lin_prec_thick",
  "model_mod_hipp_acc","model_mod_hipp_prec","model_mod_thick_acc","model_mod_thick_prec",
  "model_quad_acc_hipp","model_quad_acc_hipp_act","model_quad_acc_memor","model_quad_acc_sup_act","model_quad_acc_thick",
  "model_quad_prec_hipp","model_quad_prec_hipp_act","model_quad_prec_memor","model_quad_prec_sup_act","model_quad_prec_thick"
)


for (m in models) {{
    mod <- get(m, inherits = FALSE)                         
    mod <- tryCatch(update(mod, na.action = na.exclude),    
                    error = function(e) mod)
    res <- residuals(mod)                                   
    print(shapiro.test(res))                                
    suppressWarnings(print(ks.test(scale(res), "pnorm")))   
    
    par(mfrow = c(1,2))
    plot(fitted(mod), res, main = paste("Residuals vs Fitted:", m),
         xlab = "Fitted values", ylab = "Residuals"); abline(h = 0, lty = 2)
    qqnorm(res, main = paste("Q–Q Plot:", m)); qqline(res)
    par(mfrow = c(1,1))
    
  }, error = function(e) {
    cat("Could not process", m, "→", e$message, "\n")
  }
  flush.console()
}

