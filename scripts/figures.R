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
                   "Hippocampal volume (score z)", "ACC Glu (mM)",
                   reverse_x = TRUE)

p2 <- base_scatter(df_prec, "hip_l_nor_icv_z", "m_m_precuneus",
                   "Hippocampal volume (score z)", "Precuneus Glu (mM)",
                   reverse_x = TRUE)

p3 <- base_scatter(df_acc,  "ct_ad_z", "m_m_acc",
                   "Cortical thickness AD (score z)", "ACC Glu (mM)",
                   reverse_x = TRUE)

p4 <- base_scatter(df_prec, "ct_ad_z", "m_m_precuneus",
                   "Cortical thickness AD (score z)", "Precuneus Glu (mM)",
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
  "Superior parietal activation (score z)", "ACC Glu (mM)",
  fit_formula = y ~ x                    # linear
)

p_prec <- base_scatter(
  df_act, "act_parietal_z", "m_m_precuneus",
  "Superior parietal activation (score z)", "Precuneus Glu (mM)",
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


# -----------------------------
# 1) Data prep (from your MRS_SM_Prec)
# -----------------------------
DATA_REG0 <- MRS_SM_Prec[, c(
  "m_m_precuneus_c",         # centered Glu (we'll z-score this)
  "memoria_libre_correcte",  # outcome
  "hip_l_nor_icv_c"          # centered hippocampal volume (moderator)
)] |>
  tidyr::drop_na() |>
  dplyr::mutate(
    m_m_precuneus_c        = as.numeric(m_m_precuneus_c),
    hip_l_nor_icv_c        = as.numeric(hip_l_nor_icv_c),
    memoria_libre_correcte = as.numeric(memoria_libre_correcte),
    # z-score the (centered) Glu — same z as from raw
    mmprec_z               = as.numeric(scale(m_m_precuneus_c))
  )

# -----------------------------
# 2) Fit moderation model in z-units
# -----------------------------
model_mod_hipp_prec_z <- lm(
  memoria_libre_correcte ~ mmprec_z + hip_l_nor_icv_c:mmprec_z,
  data = DATA_REG0
)
print(summary(model_mod_hipp_prec_z))

# -----------------------------
# 3) Prediction grid over z(Glu) with H at mean, ±1 SD
# -----------------------------
z_rng  <- range(DATA_REG0$mmprec_z, na.rm = TRUE)
z_grid <- seq(z_rng[1], z_rng[2], length.out = 200)

H_mean <- mean(DATA_REG0$hip_l_nor_icv_c, na.rm = TRUE)  # ~0 by centering
H_sd   <- sd(DATA_REG0$hip_l_nor_icv_c,   na.rm = TRUE)

mk_curve <- function(label, Hval) {
  nd <- data.frame(
    mmprec_z        = as.numeric(z_grid),
    hip_l_nor_icv_c = rep(as.numeric(Hval), length(z_grid)),
    stringsAsFactors = FALSE
  )
  pr <- predict(model_mod_hipp_prec_z, newdata = nd, se.fit = TRUE)
  crit <- qnorm(0.975)
  data.frame(
    H_level = label,
    Xz      = z_grid,
    Yhat    = pr$fit,
    lo      = pr$fit - crit * pr$se.fit,
    hi      = pr$fit + crit * pr$se.fit
  )
}

df_plot <- dplyr::bind_rows(
  mk_curve("High (+1 SD)", H_mean + H_sd),
  mk_curve("Mean (0 SD)",  H_mean),
  mk_curve("Low (−1 SD)",  H_mean - H_sd)
)

# -----------------------------
# 4) Bin raw points into Low/Mean/High by hippocampal volume (±1 SD)
# -----------------------------
DATA_REG0 <- DATA_REG0 |>
  mutate(
    Group = dplyr::case_when(
      hip_l_nor_icv_c <= (H_mean - H_sd) ~ "Low (−1 SD)",
      hip_l_nor_icv_c >= (H_mean + H_sd) ~ "High (+1 SD)",
      TRUE                                ~ "Mean (0 SD)"
    )
  )

# Harmonize factor levels for single legend
lvl <- c("High (+1 SD)", "Mean (0 SD)", "Low (−1 SD)")
DATA_REG0$Group <- factor(DATA_REG0$Group, levels = lvl)
df_plot$Group   <- factor(df_plot$H_level, levels = lvl)

# -----------------------------
# 5) Axis ticks (Y)
# -----------------------------
y_min    <- floor(min(DATA_REG0$memoria_libre_correcte, na.rm = TRUE))
y_max    <- ceiling(max(DATA_REG0$memoria_libre_correcte, na.rm = TRUE))
y_breaks <- seq(y_min, y_max, by = 1)

# -----------------------------
# Plot
# -----------------------------
p_mod <- ggplot() +
  geom_point(
    data = DATA_REG0,
    aes(x = mmprec_z, y = memoria_libre_correcte, shape = Group),
    color = "grey30", size = 1.9, alpha = 0.85
  ) +
  geom_ribbon(
    data = df_plot,
    aes(x = Xz, ymin = lo, ymax = hi, fill = Group),
    alpha = 0.12, linewidth = 0
  ) +
  geom_line(
    data = df_plot,
    aes(x = Xz, y = Yhat, linetype = Group, linewidth = Group),
    colour = "black"
  ) +
  scale_x_reverse() +
  scale_y_continuous(breaks = y_breaks, limits = c(y_min, y_max),
                     expand = c(0.02, 0.02)) +
  
  # two legends: (1) lines+ribbons, (2) points/shapes
  scale_linetype_manual(
    name = "  ",
    values = c("High (+1 SD)"="dashed","Mean (0 SD)"="dotted","Low (−1 SD)"="solid"),
    breaks = lvl
  ) +
  scale_linewidth_manual(
    name = "  ",
    values = c("High (+1 SD)"=0.9,"Mean (0 SD)"=1.0,"Low (−1 SD)"=0.9),
    breaks = lvl
  ) +
  scale_fill_manual(
    name = "  ",
    values = c("High (+1 SD)"="grey60","Mean (0 SD)"="grey75","Low (−1 SD)"="grey88"),
    breaks = lvl
  ) +
  scale_shape_manual(
    name = "",
    values = c("High (+1 SD)"=16,"Mean (0 SD)"=17,"Low (−1 SD)"=15),
    breaks = lvl
  ) +
  
  labs(x = "Precuneus Glu (score z)", y = "Memoria free word recall (correct)") +
  theme_classic(base_size = 12) +
  theme(
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 11),
    legend.box = "vertical",
    plot.margin = margin(8, 8, 8, 8)
  ) +
  guides(
    shape = guide_legend(override.aes = list(size = 3, alpha = 1, colour = "grey30")),
    linewidth = "none"
  )

# -----------------------------
# 7) Export (1200-DPI TIFF)
# -----------------------------
out_dir <- "C:/Users/okkam/Desktop/labo/article 1/code/Glut_project/figures"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ragg::agg_tiff(
  filename = file.path(out_dir, "fig_mod_precuneus_Z_reversed_shapes_lines_1200dpi.tiff"),
  width = 6.5, height = 4.3, units = "in", res = 1200, compression = "lzw"
)
print(p_mod)
dev.off()



















