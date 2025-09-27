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





