################################################################################
# 03_fecundity_length.R
#
# Fecundity-length relationship for female Patagonian toothfish: the number of
# eggs carried by a stage-4 female grows exponentially with total length,
#
#     N_eggs = a * exp(b * TL)
#
# fitted by non-linear least squares on two sources: an onboard sampling
# programme run in 2024 (Albius) and 2025 (Ile de La Reunion 2), and a 1990
# dataset.
#
# Requires: 00_setup.R
#
# Author: Fanny Ouzoulias
# Date:   2026-08-19
################################################################################

# protocole is a label, not a number: read back from CSV, its 1990 / 2024 / 2025
# values would otherwise come in as doubles.
data_PCE <- readr::read_csv(dpath("fecundity_samples.csv"),
                            show_col_types = FALSE) %>%
  mutate(protocole = as.character(protocole))

# Mature females with both a length and an egg count.
data_fit <- data_PCE %>%
  filter(SEXE == "F", STADE == "4") %>%
  na.omit() %>%
  filter(is.finite(NB_OEUF), is.finite(LT))

## Fit -------------------------------------------------------------------------
fit_eggs <- nls(NB_OEUF ~ a * exp(b * LT),
                data  = data_fit,
                start = list(a = min(data_fit$NB_OEUF), b = 0.01))

summary(fit_eggs)

## Prediction with a confidence band -------------------------------------------
# The band comes from the delta method: the gradient of a*exp(b*TL) with respect
# to (a, b), combined with the parameter covariance matrix.
new_data <- tibble(LT = seq(min(data_fit$LT), max(data_fit$LT), length.out = 100))
new_data$fit <- predict(fit_eggs, newdata = new_data)

a <- coef(fit_eggs)["a"]
b <- coef(fit_eggs)["b"]
grad <- cbind(exp(b * new_data$LT), a * new_data$LT * exp(b * new_data$LT))
new_data$se    <- sqrt(rowSums((grad %*% vcov(fit_eggs)) * grad))
new_data$lower <- new_data$fit - 1.96 * new_data$se
new_data$upper <- new_data$fit + 1.96 * new_data$se

## Figure ----------------------------------------------------------------------
# The two 2024 and 2025 campaigns follow the same protocol and are shown together.
data_plot <- data_PCE %>%
  mutate(protocole = if_else(protocole %in% c("2024", "2025"), "2024-2025", protocole))

p_fecundity <- ggplot() +
  geom_point(data = data_plot, aes(x = LT, y = NB_OEUF / 1e6, colour = protocole),
             alpha = 0.6) +
  geom_ribbon(data = new_data, aes(x = LT, ymin = lower / 1e6, ymax = upper / 1e6),
              alpha = 0.2, fill = "cadetblue4") +
  geom_line(data = new_data, aes(x = LT, y = fit / 1e6),
            colour = "cadetblue4", linewidth = 1) +
  scale_x_continuous(breaks = seq(60, 200, 20)) +
  scale_y_continuous(breaks = seq(0, 1.2, 0.1)) +
  labs(x = "TL (cm)", y = "Number of eggs (millions)", colour = "Sampling") +
  theme_paper()

save_fig("fecundity_length_relationship.png", p_fecundity, width = 15, height = 10)

saveRDS(fit_eggs, file.path(out_dir, "fecundity_model.rds"))
