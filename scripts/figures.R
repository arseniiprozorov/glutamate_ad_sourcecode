library(here)
source(file.path("scripts", "source_code.R"))
ls()
names(MRS_full)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(ragg)

#### Bar plot #####

# --- long format for ACC & Precuneus ---
long_data_mM <- MRS_full %>%
  dplyr::select(diagnostic_nick, m_m_acc, m_m_precuneus) %>%
  dplyr::rename(ACC = m_m_acc, Precuneus = m_m_precuneus) %>%
  tidyr::pivot_longer(c(ACC, Precuneus),
                      names_to = "Region",
                      values_to = "Value")

# ----- summarise both regions -----
sum_dat <- long_data_mM %>%
  filter(Region %in% c("ACC","Precuneus")) %>%
  group_by(Region, diagnostic_nick) %>%
  summarise(
    mean = mean(Value, na.rm = TRUE),
    se   = sd(Value, na.rm = TRUE)/sqrt(sum(!is.na(Value))),
    ci   = qt(0.975, df = max(sum(!is.na(Value))-1,1)) * se,
    .groups = "drop"
  ) %>%
  mutate(xi = as.numeric(factor(diagnostic_nick, levels=c("HC","SCD+","MCI"))))

# --- helper to build each plot ---
make_plot <- function(region_name, y_label, sig_xmin, sig_xmax, y_pos, star_label="*",
                      y_min = 13, y_max = 16) {
  
  dat <- sum_dat %>% filter(Region == region_name)
  
  sig_df <- data.frame(
    xmin = sig_xmin,
    xmax = sig_xmax,
    y = y_pos,
    label = rep(star_label, length(y_pos))
  )
  
  ggplot(dat, aes(x = xi, fill = diagnostic_nick)) +
    # bars from baseline to mean
    geom_rect(aes(xmin = xi - 0.35,
                  xmax = xi + 0.35,
                  ymin = y_min,
                  ymax = mean,
                  fill = diagnostic_nick),
              colour = "black") +
    # error bars (95% CI)
    geom_errorbar(aes(x = xi, ymin = mean, ymax = mean + ci),
                  width = 0.2, color = "black") +
    # significance bars
    geom_segment(data = sig_df,
                 aes(x = xmin, xend = xmax, y = y, yend = y),
                 inherit.aes = FALSE, linewidth = 0.8) +
    geom_text(data = sig_df,
              aes(x = (xmin + xmax)/2, y = y + 0.05, label = label),
              inherit.aes = FALSE, size = 5) +
    scale_fill_manual(values = c("HC"="#C8C8C8","SCD+"="#9A9A9A","MCI"="#6B6B6B"),
                      name = "Group") +
    scale_x_continuous(breaks = 1:3,
                       labels = c("HC","SCD+","MCI"),
                       limits = c(0.6,3.4),
                       expand = c(0,0)) +
    scale_y_continuous(breaks = seq(y_min, y_max, 1),
                       limits = c(y_min, y_max),
                       expand = c(0,0)) +
    labs(x = "", y = y_label) +
    theme_classic(base_size = 14) +
    theme(
      axis.title.y = element_text(size=16),
      axis.text.x  = element_text(size=14),
      axis.text.y  = element_text(size=14),
      plot.margin  = margin(t=26, r=18, b=10, l=22),
      legend.position = "bottom"
    )
}

# --- build panels  ---
p_acc  <- make_plot("ACC","ACC Glu (mM)",
                    sig_xmin = c(1,1), sig_xmax = c(2,3),
                    y_pos = c(15.5, 15.7))

p_prec <- make_plot("Precuneus","Precuneus Glu (mM)",
                    sig_xmin = 1, sig_xmax = 3,
                    y_pos = 15.8)

panel <- p_acc + p_prec + plot_layout(guides = "collect") &
  theme(legend.position = "right")

# add A/B tags
panel <- panel + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 12),
        plot.tag.position = c(0.02, 0.98))

# --- export to your project root figures folder ---
out_dir <- "C:/Users/okkam/Desktop/labo/article 1/code/Glut_project/figures"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ragg::agg_tiff(file.path(out_dir, "fig_bar_ACC_Precuneus_panel_1200dpi.tiff"),
               width = 7, height = 4, units = "in", res = 1200, compression = "lzw")
print(panel)
dev.off()






####### Structure and Memory #####################
##########  quadratic structure  ##########

# 1) Z-scores + group factor
MRS_full <- MRS_full %>%
  mutate(
    hip_l_nor_icv_z = as.numeric(scale(hip_l_nor_icv)),
    ct_ad_z         = as.numeric(scale(cortical_thickness_adsignature_dickson)),
    group           = factor(diagnostic_nick, levels = c("HC","SCD+","MCI"))
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

# 4) Data (ensure no NAs in any plotted var)
df_acc  <- MRS_full %>%
  select(m_m_acc, m_m_precuneus, hip_l_nor_icv_z, ct_ad_z, group) %>%
  tidyr::drop_na()
df_prec <- df_acc

# 5) Panels (ALL quadratic; keep reversed x as you had)
p1 <- base_scatter(df_acc,  "hip_l_nor_icv_z", "m_m_acc",
                   "Hippocampal volume (Z score)", "ACC Glu (mM)",
                   reverse_x = TRUE)

p2 <- base_scatter(df_prec, "hip_l_nor_icv_z", "m_m_precuneus",
                   "Hippocampal volume (Z score)", "Precuneus Glu (mM)",
                   reverse_x = TRUE)

p3 <- base_scatter(df_acc,  "ct_ad_z", "m_m_acc",
                   "Cortical thickness (Z score)", "ACC Glu (mM)",
                   reverse_x = TRUE)

p4 <- base_scatter(df_prec, "ct_ad_z", "m_m_precuneus",
                   "Cortical thickness (Z score)", "Precuneus Glu (mM)",
                   reverse_x = TRUE)

# 6) Arrange 2x2, ONE legend
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
out_dir <- "C:/Users/okkam/Desktop/labo/article 1/code/Glut_project/figures"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ragg::agg_tiff(file.path(out_dir, "fig_glu_four_panels_quad_1200dpi.tiff"),
               width = 7.0, height = 6.0, units = "in",
               res = 1200, compression = "lzw")
print(final_plot)
dev.off()







##### Activaitons #############
# 1) Z-scores + group factor (parietal activation) — unchanged
MRS_full <- MRS_full %>%
  mutate(
    act_parietal_z = as.numeric(scale(activation_parietal_sup_l)),
    group          = factor(diagnostic_nick, levels = c("HC","SCD+","MCI"))
  )

# 2) Common Y limits across metabolites — unchanged
y_rng <- range(MRS_full$m_m_acc, MRS_full$m_m_precuneus, na.rm = TRUE)
y_lim <- c(floor(y_rng[1] - 0.2), ceiling(y_rng[2] + 0.2))

# 3) Helper: accepts a formula (default quadratic). Set se=FALSE to match your style.
base_scatter <- function(df, x, y, xlab, ylab,
                         fit_formula = y ~ x + I(x^2)) {
  ggplot(df, aes(x = .data[[x]], y = .data[[y]])) +
    geom_point(aes(shape = group), color = "grey30", size = 1.4, alpha = 0.8) +
    stat_smooth(method = "lm", formula = fit_formula,
                se = FALSE, color = "black", linewidth = 0.8) +
    scale_shape_manual(values = c(16, 17, 15), name = "Group",
                       labels = c("HC","SCD+","MCI")) +
    labs(x = xlab, y = ylab) +
    coord_cartesian(ylim = y_lim, clip = "off") +
    theme_classic(base_size = 12) +
    theme(
      plot.margin   = margin(8, 8, 8, 8),
      legend.position = "right",
      legend.title  = element_text(size = 12),
      legend.text   = element_text(size = 11)
    )
}

# 4) Data (drop NAs on needed cols) — unchanged
df_act <- MRS_full %>%
  select(m_m_acc, m_m_precuneus, act_parietal_z, group) %>%
  tidyr::drop_na()

# 5) Two panels (ACC = linear; Precuneus = quadratic)
p_acc  <- base_scatter(
  df_act, "act_parietal_z", "m_m_acc",
  "Superior parietal activation (Z score)", "ACC Glu (mM)",
  fit_formula = y ~ x                    # linear
)

p_prec <- base_scatter(
  df_act, "act_parietal_z", "m_m_precuneus",
  "Superior parietal activation (Z score)", "Precuneus Glu (mM)",
  fit_formula = y ~ x + I(x^2)           # quadratic
)

# 6) Arrange side-by-side, one legend + A/B tags — unchanged
combo <- (p_acc | p_prec) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A")

combo <- combo & theme(
  legend.position = "right",
  plot.tag = element_text(face = "bold", size = 12),
  plot.tag.position = c(0.02, 0.98)
)

# 7) Export
out_dir <- "C:/Users/okkam/Desktop/labo/article 1/code/Glut_project/figures"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ragg::agg_tiff(file.path(out_dir, "fig_glu_vs_actpar_ACClin_PRECquad_1200dpi.tiff"),
               width = 7.0, height = 2.6, units = "in",
               res = 1200, compression = "lzw")
print(combo)
dev.off()






#### Moderation #####

DATA_REG0 <- MRS_full[, c(
  "m_m_precuneus",          
  "memoria_libre_correcte", 
  "hip_l_nor_icv_c"         
)] |> tidyr::drop_na()

# ---- 2) Transform: G in Z; H centered (already) ----
DATA_REG0$mmprec_z <- as.numeric(scale(DATA_REG0$m_m_precuneus))  
DATA_REG0$hip_vec  <- as.numeric(DATA_REG0$hip_l_nor_icv_c)       

# ---- 3) Linear moderation model: M ~ Gz + H:Gz ----
model_mod_precuneus_z <- lm(
  memoria_libre_correcte ~ mmprec_z  + hip_vec:mmprec_z,
  data = DATA_REG0
)
print(summary(model_mod_precuneus_z))

# ---- 4) Prediction grid over G (Z) with H at mean ±1 SD and mean ----
rng_z  <- range(DATA_REG0$mmprec_z, na.rm = TRUE)
grid_z <- seq(rng_z[1], rng_z[2], length.out = 200)

H_mean <- mean(DATA_REG0$hip_vec, na.rm = TRUE)  # ~0 by centering
H_sd   <- sd(DATA_REG0$hip_vec,   na.rm = TRUE)

mk_curve <- function(label, Hval) {
  pr <- predict(model_mod_precuneus_z,
                newdata = data.frame(mmprec_z = grid_z, hip_vec = Hval),
                se.fit = TRUE)
  crit <- qnorm(0.975)
  data.frame(
    H_level = label,
    G_z     = grid_z,
    M_hat   = pr$fit,
    lo      = pr$fit - crit*pr$se.fit,
    hi      = pr$fit + crit*pr$se.fit
  )
}

df_plot <- dplyr::bind_rows(
  mk_curve("High (+1 SD)", H_mean + H_sd),
  mk_curve("Mean (0 SD)",  H_mean),
  mk_curve("Low (−1 SD)",  H_mean - H_sd)
)
df_plot$H_level <- factor(
  df_plot$H_level,
  levels = c("High (+1 SD)", "Mean (0 SD)", "Low (−1 SD)")
)

# ---- 5) Y axis ticks at whole numbers ----
y_min    <- floor(min(DATA_REG0$memoria_libre_correcte, na.rm = TRUE))
y_max    <- ceiling(max(DATA_REG0$memoria_libre_correcte, na.rm = TRUE))
y_breaks <- seq(y_min, y_max, by = 1)

# ---- 6) Plot ----
p_mod <- ggplot() +
  # raw points
  geom_point(
    data = DATA_REG0,
    aes(x = mmprec_z, y = memoria_libre_correcte),
    color = "grey30", size = 1.6, alpha = 0.8, shape = 16
  ) +
  # ribbons (95% CI) and lines
  geom_ribbon(
    data = df_plot,
    aes(x = G_z, ymin = lo, ymax = hi, fill = H_level),
    alpha = 0.12, linewidth = 0
  ) +
  geom_line(
    data = df_plot,
    aes(x = G_z, y = M_hat, linetype = H_level, linewidth = H_level),
    colour = "black"
  ) +
  # reverse x to match your convention
  scale_x_reverse() +
  scale_y_continuous(breaks = y_breaks, limits = c(y_min, y_max), expand = c(0.02, 0.02)) +
  # line styles: High dashed, Mean solid, Low dotted
  scale_linetype_manual(values = c(
    "High (+1 SD)" = "dashed",  
    "Mean (0 SD)"  = "solid",
    "Low (−1 SD)"  = "dotted"
  )) +
  scale_linewidth_manual(values = c(
    "High (+1 SD)" = 0.9,
    "Mean (0 SD)"  = 1.0,
    "Low (−1 SD)"  = 0.9
  )) +
  # neutral greys for ribbons
  scale_fill_manual(values = c(
    "High (+1 SD)" = "grey60",
    "Mean (0 SD)"  = "grey75",
    "Low (−1 SD)"  = "grey88"
  )) +
  labs(x = "Precuneus Glu (Z score)", y = "Free recall (correct)") +
  theme_classic(base_size = 12) +
  theme(
    plot.margin     = margin(8, 8, 8, 8),
    legend.position = "right",
    legend.title    = element_blank(),
    legend.text     = element_text(size = 11)
  )

# ---- 7) Export ----
out_dir <- "C:/Users/okkam/Desktop/labo/article 1/code/Glut_project/figures"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ragg::agg_tiff(file.path(out_dir, "fig_mod_precuneus_Z_reversed_high_dashed_1200dpi.tiff"),
               width = 6.5, height = 4.3, units = "in", res = 1200, compression = "lzw")
print(p_mod)
dev.off()
