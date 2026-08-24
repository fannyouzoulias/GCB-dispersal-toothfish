################################################################################
# 06_retention_front_glm.R
#
# Relationship between front intensity and larval retention.
#
# The response is the annual number of retained eggs, R, summed over the four
# recruitment sectors. A Gaussian GLM on log R is a lognormal model on R, so the
# fitted relationship is exponential on the original scale: R = exp(a + b z).
#
#   - log R against standardized PF intensity  (main figure, and Table A3)
#   - log R against standardized SAF intensity (Table A4)
#   - regression diagnostics of the PF model
#
# Requires: 00_setup.R, and the CSVs written by 04_retention_recruitment.R
#           and 05_front_indices.R
#
# Author: Fanny Ouzoulias
# Date:   2026-08-19
################################################################################

## Annual table ----------------------------------------------------------------
dat <- read_csv(file.path(out_dir, "recruited_annual.csv"), show_col_types = FALSE) %>%
  left_join(read_csv(file.path(out_dir, "front_indices_annual.csv"),
                     show_col_types = FALSE), by = "Year") %>%
  drop_na(Recruited_total, PF_mean_cm_s, SAF_mean_cm_s) %>%
  mutate(
    logR  = log(Recruited_total),
    z_PF  = as.numeric(scale(PF_mean_cm_s)),
    z_SAF = as.numeric(scale(SAF_mean_cm_s))
  )

## Polar Front model -----------------------------------------------------------
m_pf <- glm(logR ~ z_PF, data = dat, family = gaussian())
summary(m_pf)

# Share of the variance of log R explained.
r2     <- 1 - deviance(m_pf) / m_pf$null.deviance
r2_adj <- 1 - (1 - r2) * ((nrow(dat) - 1) / (nrow(dat) - length(coef(m_pf))))
message("PF model: R2 = ", round(r2, 3), ", adjusted R2 = ", round(r2_adj, 3))

## Subantarctic Front model ----------------------------------------------------
m_saf <- glm(logR ~ z_SAF, data = dat, family = gaussian())
summary(m_saf)

## Main figure: effect of PF intensity on retention ----------------------------
grid <- tibble(z_PF = seq(min(dat$z_PF), max(dat$z_PF), length.out = 200))
pred <- predict(m_pf, newdata = grid, se.fit = TRUE)

grid <- grid %>%
  mutate(fit_R = exp(pred$fit),
         lo_R  = exp(pred$fit - 1.96 * pred$se.fit),
         hi_R  = exp(pred$fit + 1.96 * pred$se.fit))

p_pf_effect <- ggplot(dat, aes(x = z_PF, y = Recruited_total)) +
  geom_ribbon(data = grid, aes(x = z_PF, ymin = lo_R, ymax = hi_R),
              alpha = 0.20, inherit.aes = FALSE) +
  geom_line(data = grid, aes(x = z_PF, y = fit_R),
            linewidth = 1.1, inherit.aes = FALSE) +
  geom_point(size = 2.6, alpha = 0.9) +
  labs(x = "Polar Front intensity", y = "Number of \nlarvae retained") +
  theme_paper() +
  theme(panel.grid = element_blank())

save_fig("pf_effect_on_retention.png", p_pf_effect, width = 8, height = 6)

## Diagnostics of the PF model -------------------------------------------------
diag_df <- dat %>%
  mutate(fit       = fitted(m_pf),
         resid_std = rstandard(m_pf),
         sqrt_abs  = sqrt(abs(resid_std)),
         fit_R     = exp(fitted(m_pf))) %>%
  filter(is.finite(fit), is.finite(resid_std))

p1 <- ggplot(diag_df, aes(fit, resid_std)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE, colour = "black") +
  labs(title = "Residuals vs fitted", x = "Fitted values (logR)",
       y = "Standardized residuals") + theme_bw()

p2 <- ggplot(diag_df, aes(sample = resid_std)) +
  stat_qq(size = 1.8, alpha = 0.8) + stat_qq_line(colour = "black") +
  labs(title = "Normal Q-Q plot", x = "Theoretical quantiles",
       y = "Standardized residuals") + theme_bw()

p3 <- ggplot(diag_df, aes(fit, sqrt_abs)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE, colour = "black") +
  labs(title = "Scale-location", x = "Fitted values (logR)",
       y = expression(sqrt("|Standardized residuals|"))) + theme_bw()

p4 <- ggplot(diag_df, aes(fit_R, Recruited_total)) +
  geom_point(size = 2, alpha = 0.85) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_x_log10() + scale_y_log10() +
  labs(title = "Observed vs fitted (original scale)",
       x = "Fitted recruited", y = "Observed recruited") + theme_bw()

save_fig("diagnostics_pf_model.png", (p1 | p2) / (p3 | p4), width = 8, height = 6)

## Coefficient tables ----------------------------------------------------------
# These are the appendix tables of the paper.
coef_table <- function(model, label) {
  as.data.frame(summary(model)$coefficients) %>%
    tibble::rownames_to_column("Parameter") %>%
    mutate(Model = label)
}

readr::write_csv(bind_rows(coef_table(m_pf, "logR ~ z(I_PF)"),
                           coef_table(m_saf, "logR ~ z(I_SAF)")),
                 file.path(out_dir, "glm_coefficients.csv"))
