library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(purrr)
library(ragg)





#### Figure 2 Neurodegeneration and GLU
#Old figure 4 quadratic

# 1) Z-scores 
MRS_full <- MRS_full %>%
  mutate(
    hipp_mean_z = as.numeric(scale(hipp_mean)),
    ct_ad_z     = as.numeric(scale(cortical_thickness_adsignature_dickson)),
    group       = factor(diagnostic_nick, levels = c("HC", "SCD+", "MCI"))
  )

# 2) Common Y limits
y_rng <- range(MRS_full$m_m_acc, MRS_full$m_m_precuneus, na.rm = TRUE)
y_lim <- c(floor(y_rng[1] - 0.2), ceiling(y_rng[2] + 0.2))

# 3) Helper (quadratic by default)
base_scatter <- function(df, x, y, xlab, ylab, reverse_x = FALSE) {
  p <- ggplot(df, aes(x = .data[[x]], y = .data[[y]])) +
    geom_point(aes(shape = group), color = "grey30", size = 1.4, alpha = 0.8) +
    stat_smooth(method = "lm", formula = y ~ x + I(x^2),
                se = FALSE, color = "black", linewidth = 0.8) +
    scale_shape_manual(values = c(16,17,15), name = "Group",
                       labels = c("HC","SCD+","MCI")) +
    labs(x = xlab, y = ylab) +
    coord_cartesian(ylim = y_lim, clip = "off") +
    theme_classic(base_size = 12, base_family = "Arial") +
    theme(
      plot.margin = margin(8, 8, 8, 8),
      legend.position = "right",
      legend.title = element_text(size = 12),
      legend.text  = element_text(size = 11)
    )
  if (reverse_x) p <- p + scale_x_reverse()
  p
}

# 4) Data 
df_acc  <- MRS_full %>%
  select(m_m_acc, m_m_precuneus, hipp_mean_z, ct_ad_z, group) %>%
  tidyr::drop_na()
df_prec <- df_acc

# 5) Panels (ALL quadratic)
p1 <- base_scatter(df_acc,  "hipp_mean_z", "m_m_acc",
                   "Mean hippocampal volume", "ACC Glu (mM)",
                   reverse_x = TRUE)

p2 <- base_scatter(df_prec, "hipp_mean_z", "m_m_precuneus",
                   "Mean hippocampal volume", "Precuneus Glu (mM)",
                   reverse_x = TRUE)

p3 <- base_scatter(df_acc,  "ct_ad_z", "m_m_acc",
                   "Cortical thickness AD", "ACC Glu (mM)",
                   reverse_x = TRUE)

p4 <- base_scatter(df_prec, "ct_ad_z", "m_m_precuneus",
                   "Cortical thickness AD", "Precuneus Glu (mM)",
                   reverse_x = TRUE)

# 6) Arrange 2x2
combo <- (p1 | p2) / (p3 | p4) + plot_layout(guides = "collect") &
  theme(
    legend.position = "right",
    plot.tag = element_text(face = "bold", size = 12),
    plot.tag.position = c(0.02, 0.98)
  )

# optional thin right spacer, then apply tags
final_plot <- (combo | plot_spacer()) + plot_layout(widths = c(1, 0.06))
final_plot <- final_plot + plot_annotation(tag_levels = "A") &
  theme(legend.box.margin = margin(10, 8, 0, 0))

# 7) Export to your project figures folder
out_dir <- "C:/Users/okkam/Desktop/labo/article 1/Brain/methodes_resultats"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ragg::agg_tiff(file.path(out_dir, "fig_glu_four_panels_mean_hipp_1200dpi.tiff"),
               width = 7.0, height = 6.0, units = "in",
               res = 1200, compression = "lzw")
print(final_plot)
dev.off()







### New figure 3 linear 1 quadratic 

# 1) Z-scores & group factor
MRS_fig2 <- MRS_full %>%
  mutate(
    hipp_mean_z = as.numeric(scale(hipp_mean)),
    ct_ad_z     = as.numeric(scale(cortical_thickness_adsignature_dickson)),
    group       = factor(diagnostic_nick, levels = c("HC", "SCD+", "MCI"))
  )

# 2) Common Y limits
y_rng <- range(MRS_fig2$m_m_acc, MRS_fig2$m_m_precuneus, na.rm = TRUE)
y_lim <- c(y_rng[1] - 0.3, y_rng[2] + 0.3)

# 3) Equal-proportion scatter helper with 95% CI band
base_scatter <- function(df, x, y, xlab, ylab,
                         fit_formula = y ~ x) {
  df_clean <- df %>% tidyr::drop_na(all_of(c(x, y, "group")))
  
  ggplot(df_clean, aes(x = .data[[x]], y = .data[[y]])) +
    geom_point(aes(shape = group), color = "grey35", size = 1.5, alpha = 0.8) +
    stat_smooth(method = "lm", formula = fit_formula,
                se = TRUE, fill = "grey75", alpha = 0.30,
                color = "black", linewidth = 0.9) +
    scale_shape_manual(values = c(16, 17, 15), name = "Group",
                       labels = c("CU", "SCD+", "MCI")) +
    scale_x_reverse() +
    scale_y_continuous(breaks = c(12.5, 15.0, 17.5)) +
    labs(x = xlab, y = ylab) +
    coord_cartesian(ylim = y_lim, clip = "off") +
    theme_classic(base_size = 11, base_family = "Arial") +
    theme(
      plot.margin  = margin(t = 6, r = 8, b = 6, l = 8),
      legend.title = element_text(size = 11),
      legend.text  = element_text(size = 10)
    )
}

# 4) Panels: A, B, C linear | D quadratic
pA <- base_scatter(MRS_fig2, "hipp_mean_z", "m_m_acc",
                   "Mean hippocampal volume (z-score)", "ACC Glu (mM)",
                   fit_formula = y ~ x)

pB <- base_scatter(MRS_fig2, "hipp_mean_z", "m_m_precuneus",
                   "Mean hippocampal volume (z-score)", "Precuneus Glu (mM)",
                   fit_formula = y ~ x)

pC <- base_scatter(MRS_fig2, "ct_ad_z", "m_m_acc",
                   "AD-signature cortical thickness (z-score)", "ACC Glu (mM)",
                   fit_formula = y ~ x)

pD <- base_scatter(MRS_fig2, "ct_ad_z", "m_m_precuneus",
                   "AD-signature cortical thickness (z-score)", "Precuneus Glu (mM)",
                   fit_formula = y ~ x + I(x^2))

# 5) Balanced 2x2 grid
final_plot <- (pA | pB) / (pC | pD) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    legend.position = "right",
    plot.tag        = element_text(face = "bold", size = 12)
  )

# 6) Export (1200-DPI TIFF)
out_dir <- "C:/Users/okkam/Desktop/labo/article 1/NeuroImage Clinical/figures"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ragg::agg_tiff(file.path(out_dir, "fig_glu_four_panels_mean_hipp_1200dpi.tiff"),
               width = 7.2, height = 6.0, units = "in",
               res = 1200, compression = "lzw")
print(final_plot)
dev.off()














##### Activaitons #############
library(dplyr)
library(ggplot2)
library(patchwork)

# 1) Create z-score for parietal and keep group
MRS_full <- MRS_full %>%
  mutate(
    act_parietal_z   = as.numeric(scale(activation_parietal_sup_l)),
    hipp_mean_act_z  = as.numeric(scale(hipp_mean_act_c)),
    group            = factor(diagnostic_nick, levels = c("HC","SCD+","MCI"))
  )

# 2) Common Y-range
y_rng <- range(MRS_full$m_m_acc, MRS_full$m_m_precuneus, na.rm = TRUE)
y_lim <- c(floor(y_rng[1] - 0.2), ceiling(y_rng[2] + 0.2))

# 3) Base scatter (now LINEAR by default)
base_scatter <- function(df, x, y, xlab, ylab,
                         fit_formula = y ~ x) {  # ← Linear Default
  ggplot(df, aes(x = .data[[x]], y = .data[[y]])) +
    geom_point(aes(shape = group), color = "grey30", size = 1.4, alpha = 0.8) +
    stat_smooth(method = "lm", formula = fit_formula,
                se = FALSE, color = "black", linewidth = 0.8) +
    scale_shape_manual(values = c(16, 17, 15), name = "Group",
                       labels = c("CU","SCD+","MCI")) +
    labs(x = xlab, y = ylab) +
    coord_cartesian(ylim = y_lim, clip = "off") +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "right",
      plot.margin = margin(8, 8, 8, 8)
    )
}

# 4) Select needed variables
df_plot <- MRS_full %>%
  select(m_m_acc, m_m_precuneus,
         act_parietal_z, hipp_mean_act_z, group) %>%
  tidyr::drop_na()

# 5) Create plots (all LINEAR)
p_acc_parietal <- base_scatter(
  df_plot, "act_parietal_z", "m_m_acc",
  "Superior parietal activation", "ACC Glu (mM)"
)

p_pcn_parietal <- base_scatter(
  df_plot, "act_parietal_z", "m_m_precuneus",
  "Superior parietal activation", "Precuneus Glu (mM)"
)

p_pcn_meanhipp <- base_scatter(
  df_plot, "hipp_mean_act_z", "m_m_precuneus",
  "Mean hippocampal activation", "Precuneus Glu (mM)"
)

# 6) Arrange in one figure with tags A/B/C
combo <- (p_acc_parietal | p_pcn_parietal | p_pcn_meanhipp) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A")

combo <- combo & theme(
  plot.tag = element_text(face = "bold", size = 12),
  plot.tag.position = c(0.02, 0.98)
)

# 7) Export high resolution
out_dir <- "C:/Users/okkam/Desktop/labo/article 1/Brain/methodes_resultats"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ragg::agg_tiff(file.path(out_dir, "fig_activation_Glu_linear_3panels_1200dpi.tiff"),
               width = 8, height = 2.56, units = "in",
               res = 1200, compression = "lzw")
print(combo)
dev.off()





### New figure 
# 1) Create z-score for parietal and keep group
MRS_full <- MRS_full %>%
  mutate(
    act_parietal_z   = as.numeric(scale(activation_parietal_sup_l)),
    hipp_mean_act_z  = as.numeric(scale(hipp_mean_act_c)),
    group            = factor(diagnostic_nick, levels = c("HC","SCD+","MCI"))
  )

# 2) Common Y-range
y_rng <- range(MRS_full$m_m_acc, MRS_full$m_m_precuneus, na.rm = TRUE)
y_lim <- c(floor(y_rng[1] - 0.2), ceiling(y_rng[2] + 0.2))

# 3) Base scatter (LINEAR by default + 95% CI ribbon)
base_scatter <- function(df, x, y, xlab, ylab,
                         fit_formula = y ~ x) {
  ggplot(df, aes(x = .data[[x]], y = .data[[y]])) +
    geom_point(aes(shape = group), color = "grey30", size = 1.4, alpha = 0.8) +
    stat_smooth(method = "lm", formula = fit_formula,
                se = TRUE, fill = "grey75", alpha = 0.30,
                color = "black", linewidth = 0.8) +
    scale_shape_manual(values = c(16, 17, 15), name = "Group",
                       labels = c("CU","SCD+","MCI")) +
    labs(x = xlab, y = ylab) +
    coord_cartesian(ylim = y_lim, clip = "off") +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "right",
      plot.margin = margin(8, 8, 8, 8)
    )
}

# 4) Select needed variables
df_plot <- MRS_full %>%
  select(m_m_acc, m_m_precuneus,
         act_parietal_z, hipp_mean_act_z, group) %>%
  tidyr::drop_na()

# 5) Create plots (all LINEAR)
p_acc_parietal <- base_scatter(
  df_plot, "act_parietal_z", "m_m_acc",
  "Superior parietal activation", "ACC Glu (mM)"
)

p_pcn_parietal <- base_scatter(
  df_plot, "act_parietal_z", "m_m_precuneus",
  "Superior parietal activation", "Precuneus Glu (mM)"
)

p_pcn_meanhipp <- base_scatter(
  df_plot, "hipp_mean_act_z", "m_m_precuneus",
  "Mean hippocampal activation", "Precuneus Glu (mM)"
)

# 6) Arrange in one figure with tags A/B/C
combo <- (p_acc_parietal | p_pcn_parietal | p_pcn_meanhipp) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A")

combo <- combo & theme(
  plot.tag = element_text(face = "bold", size = 12),
  plot.tag.position = c(0.02, 0.98)
)

# 7) Export high resolution
out_dir <- "C:/Users/okkam/Desktop/labo/article 1/NeuroImage Clinical/figures"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ragg::agg_tiff(file.path(out_dir, "fig_activation_Glu_linear_3panels_1200dpi.tiff"),
               width = 8, height = 2.56, units = "in",
               res = 1200, compression = "lzw")
print(combo)
dev.off()



### Memoria libre as DV

# -----------------------------
# Data preparation
# -----------------------------
DF <- MRS_full %>%
  transmute(
    memoria_libre_correcte = as.numeric(memoria_libre_correcte),
    mmprec_z               = as.numeric(scale(m_m_precuneus)),
    group                  = factor(diagnostic_nick, levels = c("HC","SCD+","MCI"))
  ) %>%
  drop_na()

# -----------------------------
# Quadratic fit
# -----------------------------
mod_quad <- lm(memoria_libre_correcte ~ mmprec_z + I(mmprec_z^2), data = DF)
summary(mod_quad)

# -----------------------------
# Axis setup
# -----------------------------
y_min    <- floor(min(DF$memoria_libre_correcte, na.rm = TRUE))
y_max    <- ceiling(max(DF$memoria_libre_correcte, na.rm = TRUE))
y_breaks <- seq(y_min, y_max, by = 1)

# -----------------------------
# Plot (non-inverted x-axis)
# -----------------------------
p_quad <- ggplot(DF, aes(x = mmprec_z, y = memoria_libre_correcte)) +
  geom_point(aes(shape = group), color = "grey30", size = 1.8, alpha = 0.85) +
  stat_smooth(
    method = "lm", formula = y ~ x + I(x^2),
    se = FALSE, color = "black", linewidth = 0.9
  ) +
  scale_shape_manual(values = c(16, 17, 15), name = "Group",
                     labels = c("HC", "SCD+", "MCI")) +
  scale_y_continuous(
    breaks = y_breaks,
    limits = c(y_min, y_max),
    expand = c(0.02, 0.02)
  ) +
  labs(
    x = "Precuneus Glu (z-score)",
    y = "Memoria free word recall (correct)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 11),
    plot.margin  = margin(8, 8, 8, 8)
  )

# -----------------------------
# Export (1200-DPI TIFF)
# -----------------------------
out_dir <- "C:/Users/okkam/Desktop/labo/article 1/Brain/methodes_resultats"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ragg::agg_tiff(
  filename = file.path(out_dir, "fig_memoria_vs_precuneusGluZ_quad_1200dpi.tiff"),
  width = 6.5, height = 4.3, units = "in", res = 1200, compression = "lzw"
)
print(p_quad)
dev.off()








## New figure

# -----------------------------
# Data preparation
# -----------------------------
DF <- MRS_full %>%
  transmute(
    memoria_libre_correcte = as.numeric(memoria_libre_correcte),
    mmprec_z               = as.numeric(scale(m_m_precuneus)),
    group                  = factor(diagnostic_nick, levels = c("HC","SCD+","MCI"))
  ) %>%
  drop_na()

# -----------------------------
# Axis setup
# -----------------------------
y_min    <- floor(min(DF$memoria_libre_correcte, na.rm = TRUE))
y_max    <- ceiling(max(DF$memoria_libre_correcte, na.rm = TRUE))
y_breaks <- seq(y_min, y_max, by = 1)

# -----------------------------
# Plot (non-inverted x-axis + 95% CI ribbon)
# -----------------------------
p_quad <- ggplot(DF, aes(x = mmprec_z, y = memoria_libre_correcte)) +
  geom_point(aes(shape = group), color = "grey30", size = 1.8, alpha = 0.85) +
  stat_smooth(
    method = "lm", formula = y ~ x + I(x^2),
    se = TRUE, fill = "grey75", alpha = 0.30,
    color = "black", linewidth = 0.9
  ) +
  scale_shape_manual(values = c(16, 17, 15), name = "Group",
                     labels = c("CU", "SCD+", "MCI")) +
  scale_y_continuous(
    breaks = y_breaks,
    expand = c(0.02, 0.02)
  ) +
  coord_cartesian(ylim = c(y_min, y_max)) +
  labs(
    x = "Precuneus Glu (z-score)",
    y = "Memoria free word recall (correct)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 11),
    plot.margin  = margin(8, 8, 8, 8)
  )

# -----------------------------
# Export (1200-DPI TIFF)
# -----------------------------
out_dir <- "C:/Users/okkam/Desktop/labo/article 1/NeuroImage Clinical/figures"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ragg::agg_tiff(
  filename = file.path(out_dir, "fig_memoria_vs_precuneusGluZ_quad_1200dpi.tiff"),
  width = 6.5, height = 4.3, units = "in", res = 1200, compression = "lzw"
)
print(p_quad)
dev.off()


