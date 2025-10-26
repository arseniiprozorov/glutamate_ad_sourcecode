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

# long format for ACC and Precuneus
long_data_mM <- MRS_full %>%
  dplyr::select(diagnostic_nick, m_m_acc, m_m_precuneus) %>%
  dplyr::rename(ACC = m_m_acc, Precuneus = m_m_precuneus) %>%
  tidyr::pivot_longer(c(ACC, Precuneus),
                      names_to = "Region",
                      values_to = "Value")

#  summarise both regions 
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

#  helper 
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

#  build panels  
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

#  export  
out_dir <- "C:\\Users\\okkam\\Desktop\\labo\\article 1\\rencontre_Sylvie_20251017"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ragg::agg_tiff(file.path(out_dir, "fig_bar_ACC_Precuneus_panel_1200dpi.tiff"),
               width = 7, height = 4, units = "in", res = 1200, compression = "lzw")
print(panel)
dev.off()






#### Violin plot version #######################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(purrr)
library(ragg)

# --- Long format for ACC and Precuneus (unchanged) ---
long_data_mM <- MRS_full %>%
  dplyr::select(diagnostic_nick, m_m_acc, m_m_precuneus) %>%
  dplyr::rename(ACC = m_m_acc, Precuneus = m_m_precuneus) %>%
  tidyr::pivot_longer(c(ACC, Precuneus),
                      names_to = "Region",
                      values_to = "Value")


# --- Helper: pointy violins + EXACT Tukey overlays (no caps) ---
make_violin_plot_A <- function(region_name, y_label,
                               y_min = 13, y_max = 16,
                               show_points = TRUE,
                               violin_bw = 0.35,
                               violin_adjust = 1,
                               headroom = 0.45,
                               box_halfwidth = 0.12,
                               median_width   = 0.26) {
  
  axis_pad  <- 0.02; box_eps <- 1e-6
  
  dat_long <- long_data_mM %>%
    dplyr::filter(Region == region_name) %>%
    dplyr::mutate(group_f = factor(diagnostic_nick, levels = c("HC","SCD+","MCI")),
                  xi_num  = as.numeric(group_f))
  
  tallest <- suppressWarnings(max(dat_long$Value, na.rm = TRUE))
  y_upper  <- max(y_max, tallest) + headroom
  y_lower  <- y_min - axis_pad
  
  # Tukey box/whisker stats (no artificial minimum whiskers)
  box_df <- dat_long %>%
    dplyr::group_by(group_f) %>%
    dplyr::summarise(
      xi   = as.numeric(first(factor(group_f, levels = c("HC","SCD+","MCI")))),
      vals = list(sort(stats::na.omit(Value))),
      .groups = "drop"
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      q1  = as.numeric(stats::quantile(vals, 0.25, names = FALSE, type = 7)),
      med = as.numeric(stats::quantile(vals, 0.50, names = FALSE, type = 7)),
      q3  = as.numeric(stats::quantile(vals, 0.75, names = FALSE, type = 7)),
      iqr = q3 - q1,
      loF = q1 - 1.5 * iqr,
      hiF = q3 + 1.5 * iqr,
      wlow  = if (length(vals) == 0) q1 else {
        vv <- vals[vals >= loF]
        if (length(vv)) min(vv) else q1
      },
      whigh = if (length(vals) == 0) q3 else {
        vv <- vals[vals <= hiF]
        if (length(vv)) max(vv) else q3
      },
      eps = 0.015,
      draw_q1 = if (iqr == 0) med - eps else q1,
      draw_q3 = if (iqr == 0) med + eps else q3,
      draw_q3d = if ((draw_q3 - draw_q1) < 1e-8) draw_q3 + box_eps else draw_q3,
      xmin = xi - box_halfwidth, xmax = xi + box_halfwidth,
      mmin = xi - median_width/2, mmax = xi + median_width/2
    ) %>% dplyr::ungroup()
  
  ggplot(dat_long, aes(x = group_f, y = Value, fill = group_f)) +
    geom_violin(trim = FALSE, scale = "width", width = 0.8,
                bw = violin_bw, adjust = violin_adjust,
                color = "black", linewidth = 0.4, na.rm = TRUE) +
    { if (show_points)
      geom_point(position = position_jitter(width = 0.07, height = 0),
                 size = 1.0, alpha = 0.45, shape = 16, stroke = 0, na.rm = TRUE)
      else NULL } +
    geom_rect(data = box_df,
              aes(xmin = xmin, xmax = xmax, ymin = draw_q1, ymax = draw_q3d),
              inherit.aes = FALSE, fill = NA, color = "black", linewidth = 0.65) +
    geom_segment(data = box_df,
                 aes(x = mmin, xend = mmax, y = med, yend = med),
                 inherit.aes = FALSE, linewidth = 0.75, color = "black") +
    geom_segment(data = box_df,
                 aes(x = xi, xend = xi, y = wlow, yend = draw_q1),
                 inherit.aes = FALSE, linewidth = 0.7, lineend = "butt", color = "black") +
    geom_segment(data = box_df,
                 aes(x = xi, xend = xi, y = draw_q3, yend = whigh),
                 inherit.aes = FALSE, linewidth = 0.7, lineend = "butt", color = "black") +
    scale_fill_manual(values = c("HC"="#D9D9D9","SCD+"="#8F8F8F","MCI"="#4A4A4A"),
                      name = "Group") +
    scale_x_discrete(labels = c("HC","SCD+","MCI"), drop = FALSE) +
    scale_y_continuous(breaks = seq(y_min, y_max, 1),
                       limits = c(y_lower, y_upper),
                       expand = expansion(mult = c(0, 0.01))) +
    coord_cartesian(clip = "off") +
    labs(x = "", y = y_label) +
    theme_classic(base_size = 14) +
    theme(
      axis.title.y = element_text(size = 16),
      axis.text.x  = element_text(size = 14),
      axis.text.y  = element_text(size = 14),
      plot.margin  = margin(t = 26, r = 18, b = 10, l = 22),
      legend.position = "bottom"
    )
}


# --- Build panels (no significance elements) ---
p_acc  <- make_violin_plot("ACC", "ACC Glu (mM)", y_min = 12, y_max = 16, show_points = TRUE)
p_prec <- make_violin_plot("Precuneus", "Precuneus Glu (mM)", y_min = 12, y_max = 16, show_points = TRUE)

panel <- p_acc + p_prec + patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "right")

# A/B tags
panel <- panel + patchwork::plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 12),
        plot.tag.position = c(0.02, 0.98))

# --- Export ---
out_dir <- "C:/Users/okkam/Desktop/labo/article 1/Brain/methodes_resultats"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ragg::agg_tiff(file.path(out_dir, "fig_violin_ACC_Precuneus_panel_1200dpi.tiff"),
               width = 7, height = 4, units = "in", res = 1200, compression = "lzw")
print(panel)
dev.off()




####### Structure and Memory #####################
##########  quadratic structure  ##########

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






## --- Peak locations (z-scale and raw-scale) for your four panels ---
# Peaks (vertex x, z-scale only)
peaks_z <- c(
  ACC_vs_Hippocampus     = vertex_x(mod_acc_hipp_z,  "hipp_mean_z", "I(hipp_mean_z^2)"),
  Precuneus_vs_Hippocampus = vertex_x(mod_prec_hipp_z, "hipp_mean_z", "I(hipp_mean_z^2)"),
  ACC_vs_Thickness       = vertex_x(mod_acc_ct_z,    "ct_ad_z",     "I(ct_ad_z^2)"),
  Precuneus_vs_Thickness = vertex_x(mod_prec_ct_z,   "ct_ad_z",     "I(ct_ad_z^2)")
)

round(peaks_z, 3)








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
                       labels = c("HC","SCD+","MCI")) +
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
               width = 10, height = 3.2, units = "in",
               res = 1200, compression = "lzw")
print(combo)
dev.off()






### Memoria libre as DV
# -----------------------------
# Libraries
# -----------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(ragg)

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
# Plot (matches your figure)
# -----------------------------
p_quad_inverted <- ggplot(DF, aes(x = mmprec_z, y = memoria_libre_correcte)) +
  geom_point(aes(shape = group), color = "grey30", size = 1.8, alpha = 0.85) +
  stat_smooth(
    method = "lm", formula = y ~ x + I(x^2),
    se = FALSE, color = "black", linewidth = 0.9
  ) +
  scale_shape_manual(values = c(16, 17, 15), name = "Group",
                     labels = c("HC", "SCD+", "MCI")) +
  scale_x_reverse() +
  scale_y_continuous(
    breaks = y_breaks,
    limits = c(y_min, y_max),
    expand = c(0.02, 0.02)
  ) +
  labs(
    x = "Precuneus Glu",
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
  filename = file.path(out_dir, "fig_memoria_vs_precuneusGluZ_quad_inverted_1200dpi.tiff"),
  width = 6.5, height = 4.3, units = "in", res = 1200, compression = "lzw"
)
print(p_quad_inverted)
dev.off()








# --- Helpers: vertex of y = b0 + b1*x + b2*x^2 ---
vertex_x <- function(mod, x_name, x2_name) {
  b <- coef(mod)
  stopifnot(all(c(x_name, x2_name) %in% names(b)))
  -b[[x_name]] / (2 * b[[x2_name]])
}

vertex_xy <- function(mod, x_name, x2_name) {
  b <- coef(mod)
  x0 <- -b[[x_name]] / (2 * b[[x2_name]])
  y0 <- b[[1]] + b[[x_name]] * x0 + b[[x2_name]] * x0^2
  c(x = x0, y = y0)
}


# --- Compute peak for your model (Precuneus Glu z -> Memoria libre) ---
vertex_x (mod_quad, "mmprec_z", "I(mmprec_z^2)")
