## figures.R
##
## Project : Hidden overweight and adverse birth outcomes among Asian mothers
##           in the United States, 2016-2023
## Purpose : Figure 2, eFigure 2 and eFigure 3
## Input   : efig2_spline_grid.csv from 02_main_models.sas (Section 7).
##           Figure 2 and eFigure 3 use values from Table 3 and from
##           02_main_models.sas Section 8a, entered below.
## Output  : Figure2, eFigure2 and eFigure3 (.tiff at 600 dpi, .png)
## Software: R 4.6.0, ggplot2 (>= 3.4.0)
## Author  : Jabin Tasnim Jahan, University at Albany, SUNY
## Contact : jtasnimjahan@albany.edu
## Version : 1.0 (November 2026), code for the submitted manuscript

library(ggplot2)

## Set this to the folder that holds efig2_spline_grid.csv
out_dir <- "path/to/your/output"

save_fig <- function(p, name, w, h) {
  ggsave(file.path(out_dir, paste0(name, ".tiff")), p, width = w, height = h,
         units = "in", dpi = 600, compression = "lzw", bg = "white")
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = w, height = h,
         units = "in", dpi = 300, bg = "white")
}

outcome_levels <- c("Hypertensive disorders", "Caesarean delivery",
                    "Preterm birth", "Low birthweight")


## ---------------------------------------------------------------------------
## Figure 2. Adjusted RRs across Asian-specific BMI categories (Table 3)
## ---------------------------------------------------------------------------

cats <- c("<18.5", "18.5–22.9\n(Reference)", "23.0–24.9",
          "25.0–27.4", "≥27.5")

f2 <- data.frame(
  outcome = rep(outcome_levels, each = 5),
  bmi = rep(cats, times = 4),
  rr  = c(0.74, 1, 1.44, 1.83, 2.79,
          0.89, 1, 1.14, 1.23, 1.41,
          1.07, 1, 1.08, 1.16, 1.37,
          1.36, 1, 0.95, 0.97, 1.02),
  lo  = c(0.70, NA, 1.40, 1.79, 2.73,
          0.88, NA, 1.13, 1.22, 1.41,
          1.04, NA, 1.06, 1.13, 1.35,
          1.33, NA, 0.93, 0.95, 1.00),
  hi  = c(0.78, NA, 1.47, 1.87, 2.85,
          0.90, NA, 1.14, 1.23, 1.42,
          1.10, NA, 1.10, 1.18, 1.39,
          1.39, NA, 0.97, 0.99, 1.04)
)
f2$bmi     <- factor(f2$bmi, levels = cats)
f2$outcome <- factor(f2$outcome, levels = outcome_levels)

pd <- position_dodge(width = 0.45)

p2 <- ggplot(f2, aes(x = bmi, y = rr, colour = outcome, shape = outcome,
                     linetype = outcome, group = outcome)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50",
             linewidth = 0.4) +
  geom_line(position = pd, linewidth = 0.7) +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = pd, width = 0.15,
                linewidth = 0.5, linetype = "solid", na.rm = TRUE,
                show.legend = FALSE) +
  geom_point(position = pd, size = 2.4, fill = "white") +
  scale_y_log10(breaks = c(0.7, 0.8, 1, 1.5, 2, 3), limits = c(0.65, 3.1)) +
  scale_colour_manual(values = c("#000000", "#D55E00", "#0072B2", "#009E73")) +
  scale_shape_manual(values = c(16, 17, 15, 18)) +
  scale_linetype_manual(values = c("solid", "longdash", "dotted", "dotdash")) +
  labs(x = expression("Prepregnancy body mass index (kg/m"^2*")"),
       y = "Adjusted risk ratio (log scale)",
       colour = NULL, shape = NULL, linetype = NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position  = "right",
        legend.key.width = unit(2, "cm"),
        axis.text        = element_text(colour = "black"),
        plot.background  = element_rect(colour = NA, fill = "white"))

save_fig(p2, "Figure2", 7.5, 4.5)


## ---------------------------------------------------------------------------
## eFigure 2. Asian and White RR curves across BMI (restricted cubic spline)
## Reference: NH White mothers at BMI 25 kg/m2
## ---------------------------------------------------------------------------

e2 <- read.csv(file.path(out_dir, "efig2_spline_grid.csv"))
names(e2) <- tolower(names(e2))

ref <- data.frame(outcome = unique(e2$outcome), group = "White",
                  bmi = 25, rr = 1, lo = 1, hi = 1)
e2 <- rbind(e2, ref)
e2 <- e2[order(e2$outcome, e2$group, e2$bmi), ]

e2$outcome <- factor(e2$outcome,
  levels = c("htn_disorder", "cesarean", "ptb", "lbw"),
  labels = outcome_levels)
e2$group <- factor(e2$group, levels = c("Asian", "White"),
  labels = c("Non-Hispanic Asian", "Non-Hispanic White"))

pe2 <- ggplot(e2, aes(x = bmi, y = rr, colour = group, fill = group,
                      linetype = group)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50",
             linewidth = 0.4) +
  geom_vline(xintercept = c(23, 27.5), linetype = "dotted",
             colour = "grey40", linewidth = 0.5) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, colour = NA,
              show.legend = FALSE) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ outcome, ncol = 2, scales = "free_y") +
  scale_y_log10(breaks = c(0.3, 0.4, 0.5, 0.6, 0.8, 1, 1.2, 1.5, 2, 2.5),
                labels = function(x) sprintf("%.1f", x)) +
  scale_x_continuous(breaks = seq(20, 40, by = 5)) +
  scale_colour_manual(values = c("#D55E00", "#000000")) +
  scale_fill_manual(values = c("#D55E00", "#000000")) +
  scale_linetype_manual(values = c("solid", "longdash")) +
  labs(x = expression("Prepregnancy body mass index (kg/m"^2*")"),
       y = "Adjusted risk ratio (log scale)",
       colour = NULL, fill = NULL, linetype = NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position  = "bottom",
        legend.key.width = unit(1.5, "cm"),
        strip.background = element_blank(),
        strip.text       = element_text(face = "bold", size = 11),
        axis.text        = element_text(colour = "black"),
        panel.spacing    = unit(1, "lines"),
        plot.background  = element_rect(colour = NA, fill = "white"))

save_fig(pe2, "eFigure2", 7.5, 6.5)


## ---------------------------------------------------------------------------
## eFigure 3. % of standard-normal-weight mothers with BMI 23.0-24.9, by age
## Values from 02_main_models.sas, Section 8a
## ---------------------------------------------------------------------------

ages <- c("20–24", "25–29", "30–34", "35–39", "≥40")

e3 <- data.frame(
  age   = rep(ages, times = 2),
  group = rep(c("All Asian", "Asian Indian"), each = 5),
  pct   = c(29.1, 29.7, 30.4, 31.2, 31.8,
            32.9, 37.1, 39.1, 39.9, 41.1)
)
e3$age   <- factor(e3$age, levels = ages)
e3$group <- factor(e3$group, levels = c("All Asian", "Asian Indian"))

pe3 <- ggplot(e3, aes(x = age, y = pct, colour = group, shape = group,
                      linetype = group, group = group)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_y_continuous(limits = c(25, 45), breaks = seq(25, 45, by = 5),
                     labels = function(x) paste0(x, "%")) +
  scale_colour_manual(values = c("#0072B2", "#D55E00")) +
  scale_shape_manual(values = c(16, 17)) +
  scale_linetype_manual(values = c("solid", "longdash")) +
  labs(x = "Maternal age (years)",
       y = "Proportion classified as overweight\nunder Asian-specific thresholds",
       colour = NULL, shape = NULL, linetype = NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position  = "top",
        legend.key.width = unit(1.5, "cm"),
        axis.text        = element_text(colour = "black"),
        plot.background  = element_rect(colour = NA, fill = "white"))

save_fig(pe3, "eFigure3", 7, 4.5)
