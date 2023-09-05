### Frequentist regression ###

# MSE for freq

mse_freq <- function (model, test.data, target.var) {
  y <- predict(model, test.data)
  if (hasName(y, "fit")) y <- y$fit
  t <- test.data[[target.var]]
  return (mean((t - y)^2))
}

# Fitting model for hospH8

freq.intcar <- lm(intcarH8 ~ arancione + intcar + newpos_av7D, data=training.set.cat)


# Normality test
shapiro.test(residuals(freq.intcar))
summary(freq.intcar)

# Diagnostics


# Fitting best bayesian model for hospH8
par(mfrow = c(2,1))
plot(freq.intcar)

# Now hospH8


freq.hospH8 <- lm(hospH8 ~ rossa + hosp + newpos_av7D, data=training.set.cat)

# Normality test
shapiro.test(residuals(freq.hospH8))
summary(freq.hospH8)

# Diagnostics


# Fitting best bayesian model for hospH8
par(mfrow = c(2,1))
plot(freq.hospH8)

mse_freq(freq.hospH8, test.set.cat, 'hospH8')
