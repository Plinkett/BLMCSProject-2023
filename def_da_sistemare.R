# Bayesian Learning and Montecarlo Simulation Project

# Group 14: Guerrini Michele                  10607256
#           Mozzi Davide                      10674709
#           Santillan Moreno Carlos Alberto   10659783

# Dataset 2.8: Covid Data (covidLom2020_21.csv)

library(BAS)

#### Functions ####

# List of useful functions used in the project

standardize <- function (data, names, centers, scales) {
  new.data <- data
  for (name in names) {
    if (is.numeric(data[1, name])) {
      new.data[[name]] <- scale(data[[name]],
                                ifelse(missing(centers), TRUE, centers[[name]]),
                                ifelse(missing(scales), TRUE, scales[[name]]))
    }
  }
  return (new.data)
}

shuffle <- function (data) {
  return (data[sample(seq(1, nrow(data)), nrow(data)),])
}

mse <- function (model, test.data, target.var, estimator = "BMA") {
  y <- predict(model, test.data, estimator = estimator)$fit
  t <- test.data[[target.var]]
  return (mean((t - y)^2))
}

color.to.indicators <- function (data) {
  new.data <- data
  new.data$bianca <- as.integer(data$color == "Bianca")
  new.data$gialla <- as.integer(data$color == "Gialla")
  new.data$arancione <- as.integer(data$color == "Arancione")
  new.data$rossa <- as.integer(data$color == "Rossa")
  return (new.data)
}

color.to.numerical <- function (data) {
  f <- function (c) {
    return (switch(c,
                   Bianca = 1,
                   Gialla = 2,
                   Arancione = 3,
                   Rossa = 4))
  }
  new.data <- data
  new.data$color.num <- sapply(data$color, f)
  return (new.data)
}

day.to.season <- function (data) {
  f <- function (d) {
    # This returns the meteorological season corresponding to date 'd':
    # winter=1, spring=2, summer=3, autumn=4
    return (((as.POSIXlt(d)$mon + 1) %/% 3) %% 4 + 1)
  }
  s <- c("winter", "spring", "summer", "autumn")
  new.data <- data
  for (i in seq(s)) {
    new.data[s[i]] <- as.integer(f(data$day) == i)
  }
  return (new.data)
}


#### Data acquisition ####

# Read data from CSV file
covid.data <- read.csv("covidLom2020_21.csv")

# Remove columns that are never used
covid.data$X <- NULL

# Set useful constants and parameters
n <- nrow(covid.data)
k <- 10  # Number of "folds" in k-fold cross-validation

# Shuffle the data
shuffled.data <- shuffle(covid.data)

# Standardize the data
# The targets are not standardized. 
to.standardize <- c("newpos", "intcar", "hosp", "newpos_av7D")
standardized.data <- standardize(shuffled.data, to.standardize)

# Transform the "color" categorical variable in several indicator variables
data.cat <- color.to.indicators(standardized.data)

# Transform the "color" categorical variable in a numerical variable, assuming that
# "Bianca" < "Gialla" < "Arancione" < "Rossa" and that the value increases linearly
data.num <- color.to.numerical(standardized.data)

# Extract the "season" categorical value (and conver it into indicator variables) from "day"
data.cat.szn <- day.to.season(data.cat)


#### Testing different priors ####
# Note: here we are training only the full model each time (the one with all the covariates).

## Default prior of the BAS library ("ZS-null"), which is an approximation of the Zellner-Siow prior

# Target: hospH8

lm.ZSnull.hosp <- bas.lm(hospH8 ~ gialla + arancione + rossa + newpos + intcar + hosp + newpos_av7D,
                    data.cat,
                    include.always = ~ .,
                    n.models = 1)
summary(lm.ZSnull.hosp)

coef(lm.ZSnull.hosp)

par(mfrow = c(2, 4))
plot(coef(lm.ZSnull.hosp), ask = FALSE)
par(mfrow = c(1, 1))

plot(confint(coef(lm.ZSnull.hosp)))

# Target: intcarH8

lm.ZSnull.intcar <- bas.lm(intcarH8 ~ gialla + arancione + rossa + newpos + intcar + hosp + newpos_av7D,
                           data.cat,
                           include.always = ~ .,
                           n.models = 1)
summary(lm.ZSnull.intcar)

coef(lm.ZSnull.intcar)

par(mfrow = c(2, 4))
plot(coef(lm.ZSnull.intcar), ask = FALSE)
par(mfrow = c(1, 1))

plot(confint(coef(lm.ZSnull.intcar)))


## Zellner's informative g-prior

# Target: hospH8

lm.gprior.hosp <- bas.lm(hospH8 ~ gialla + arancione + rossa + newpos + intcar + hosp + newpos_av7D,
                         data.cat,
                         prior = "g-prior",
                         alpha = 100,
                         include.always = ~ .,
                         n.models = 1)
summary(lm.gprior.hosp)

coef(lm.gprior.hosp)

par(mfrow = c(2, 4))
plot(coef(lm.gprior.hosp), ask = FALSE) # TODO: subset = 2:8 ?
par(mfrow = c(1, 1))

plot(confint(coef(lm.gprior.hosp)))

# Target: intcarH8

lm.gprior.intcar <- bas.lm(intcarH8 ~ gialla + arancione + rossa + newpos + intcar + hosp + newpos_av7D,
                           data.cat,
                           prior = "g-prior",
                           alpha = 100,
                           include.always = ~ .,
                           n.models = 1)
summary(lm.gprior.intcar)

coef(lm.gprior.intcar)

par(mfrow = c(2, 4))
plot(coef(lm.gprior.intcar), ask = FALSE)
par(mfrow = c(1, 1))

plot(confint(coef(lm.gprior.intcar)))


#### Model selection using BIC ####

# Target: hospH8

lm.bic.hosp <- bas.lm(hospH8 ~ gialla + arancione + rossa + newpos + intcar + hosp + newpos_av7D,
                      data.cat,
                      prior = "BIC")
summary(lm.bic.hosp)

plot(seq(coef(lm.bic.hosp)$probne0), coef(lm.bic.hosp)$probne0,
     type = "h",
     lwd = 4,
     xaxt = "n",
     main = "Feature inclusion probabilities",
     xlab = "",
     ylab = "post p(B != 0)")
axis(1, seq(coef(lm.bic.hosp)$probne0), labels = coef(lm.bic.hosp)$namesx)

image(lm.bic.hosp, rotate = FALSE)

# Target: intcarH8

lm.bic.intcar <- bas.lm(intcarH8 ~ gialla + arancione + rossa + newpos + intcar + hosp + newpos_av7D,
                      data.cat,
                      prior = "BIC")
summary(lm.bic.intcar)

plot(seq(coef(lm.bic.intcar)$probne0), coef(lm.bic.intcar)$probne0,
     type = "h",
     lwd = 4,
     xaxt = "n",
     main = "Feature inclusion probabilities",
     xlab = "",
     ylab = "post p(B != 0)")
axis(1, seq(coef(lm.bic.intcar)$probne0), labels = coef(lm.bic.intcar)$namesx)

image(lm.bic.intcar, rotate = FALSE)


#### Posterior analysis ####
# Here we analyze the models trained in the model selection paragraph

# Target: hospH8 (using BMA)

par(mfrow = c(2, 4))
plot(coef(lm.bic.hosp), ask = FALSE)
par(mfrow = c(1, 1))

plot(confint(coef(lm.bic.hosp)))

# Target: intcarH8 (using HPM)

par(mfrow = c(2, 4))
plot(coef(lm.bic.intcar, estimator = "HPM"), ask = FALSE)
par(mfrow = c(1, 1))

plot(confint(coef(lm.bic.intcar, estimator = "HPM")))


#### Appendix: different representations for some of the covariates ####

## Color as a numerical variable

# Target: hospH8

lm.numcolor.hosp <- bas.lm(hospH8 ~ color.num + newpos + intcar + hosp + newpos_av7D,
                      data.num,
                      prior = "BIC")
summary(lm.numcolor.hosp)

plot(seq(coef(lm.numcolor.hosp)$probne0), coef(lm.numcolor.hosp)$probne0,
     type = "h",
     lwd = 4,
     xaxt = "n",
     main = "Feature inclusion probabilities",
     xlab = "",
     ylab = "post p(B != 0)")
axis(1, seq(coef(lm.numcolor.hosp)$probne0), labels = coef(lm.numcolor.hosp)$namesx)

image(lm.numcolor.hosp, rotate = FALSE)

# Target: intcarH8

lm.numcolor.intcar <- bas.lm(intcarH8 ~ color.num + newpos + intcar + hosp + newpos_av7D,
                           data.num,
                           prior = "BIC")
summary(lm.numcolor.intcar)

plot(seq(coef(lm.numcolor.intcar)$probne0), coef(lm.numcolor.intcar)$probne0,
     type = "h",
     lwd = 4,
     xaxt = "n",
     main = "Feature inclusion probabilities",
     xlab = "",
     ylab = "post p(B != 0)")
axis(1, seq(coef(lm.numcolor.intcar)$probne0), labels = coef(lm.numcolor.intcar)$namesx)

image(lm.numcolor.intcar, rotate = FALSE)

## Season categorical variable (derived from day)

# Note: we can for sure ignore the "autumn" indicator variable since in our data it's always 0.
# Then we can remove one of the remaining three in order to get independent conditions.

# Target: hospH8

lm.szn.hosp <- bas.lm(hospH8 ~ gialla + arancione + rossa + newpos + intcar + hosp + newpos_av7D + winter + summer,
                      data.cat.szn,
                      prior = "BIC")
summary(lm.szn.hosp)

plot(seq(coef(lm.szn.hosp)$probne0), coef(lm.szn.hosp)$probne0,
     type = "h",
     lwd = 4,
     xaxt = "n",
     main = "Feature inclusion probabilities",
     xlab = "",
     ylab = "post p(B != 0)")
axis(1, seq(coef(lm.szn.hosp)$probne0), labels = coef(lm.szn.hosp)$namesx)

image(lm.szn.hosp, rotate = FALSE)

# Target: intcarH8

lm.szn.intcar <- bas.lm(intcarH8 ~ gialla + arancione + rossa + newpos + intcar + hosp + newpos_av7D + winter + summer,
                      data.cat.szn,
                      prior = "BIC")
summary(lm.szn.intcar)

plot(seq(coef(lm.szn.intcar)$probne0), coef(lm.szn.intcar)$probne0,
     type = "h",
     lwd = 4,
     xaxt = "n",
     main = "Feature inclusion probabilities",
     xlab = "",
     ylab = "post p(B != 0)")
axis(1, seq(coef(lm.szn.intcar)$probne0), labels = coef(lm.szn.intcar)$namesx)

image(lm.szn.intcar, rotate = FALSE)


#### Appendix: comparison wiht frequentist linear regression ####

# # TODO: check this part
# 
# # MSE for freq
# 
# mse_freq <- function (model, test.data, target.var) {
#   y <- predict(model, test.data)
#   t <- test.data[[target.var]]
#   return (mean((t - y)^2))
# }
# 
# # Fitting model for intcarH8
# 
# freq.intcar <- lm(intcarH8 ~ arancione + intcar + newpos_av7D, data=training.set.cat)
# 
# 
# # Normality test
# shapiro.test(residuals(freq.intcar))
# summary(freq.intcar)
# 
# # Diagnostics
# 
# 
# # Fitting best bayesian model for intcarH8
# par(mfrow = c(2,1))
# plot(freq.intcar)
# 
# # Now hospH8
# 
# 
# freq.hospH8 <- lm(hospH8 ~ rossa + hosp + newpos_av7D, data=training.set.cat)
# 
# # Normality test
# shapiro.test(residuals(freq.hospH8))
# summary(freq.hospH8)
# 
# # Diagnostics
# 
# 
# # Fitting best bayesian model for hospH8
# par(mfrow = c(2,1))
# plot(freq.hospH8)
# 
# mse_freq(freq.hospH8, test.set.cat, 'hospH8')


#### Cross-validation ####

mse.ZSnull.hosp <- numeric(k)
mse.ZSnull.intcar <- numeric(k)
mse.g.hosp <- numeric(k)
mse.g.intcar <- numeric(k)
mse.bic.bma.hosp <- numeric(k)
mse.bic.hpm.hosp <- numeric(k)
mse.bic.bma.intcar <- numeric(k)
mse.bic.hpm.intcar <- numeric(k)
mse.numcolor.hosp <- numeric(k)
mse.numcolor.intcar <- numeric(k)
mse.szn.hosp <- numeric(k)
mse.szn.intcar <- numeric(k)

idx <- 1
n.val <- n %/% k
for (i in 1:k) {
  test.seq <- seq(idx, idx + n.val - 1)
  idx <- idx + n.val
  if (i == k - n %% k) n.val <- n.val + 1
  
  training.set <- shuffled.data[-test.seq,]
  validation.set <- shuffled.data[test.seq,]
  
  # Standardize the training and validation sets
  # The validation set is standardized using the same center and scale values used for the training set.
  
  training.set <- standardize(training.set, to.standardize)
  centers <- list()
  scales <- list()
  for (name in to.standardize) {
    centers[[name]] <- attr(training.set[[name]], "scaled:center")
    scales[[name]] <- attr(training.set[[name]], "scaled:scale")
  }
  validation.set <- standardize(validation.set, to.standardize, centers, scales)
  
  # Transform the "color" categorical variable in several indicator variables
  training.set.cat <- color.to.indicators(training.set)
  validation.set.cat <- color.to.indicators(validation.set)
  
  # Transform the "color" categorical variable in a numerical variable, assuming that
  # "Bianca" < "Gialla" < "Arancione" < "Rossa" and that the value increases linearly
  training.set.num <- color.to.numerical(training.set)
  validation.set.num <- color.to.numerical(validation.set)
  
  # Extract the "season" categorical value (and conver it into indicator variables) from "day"
  training.set.cat.szn <- day.to.season(training.set.cat)
  validation.set.cat.szn <- day.to.season(validation.set.cat)
  
  # Train the models and compute the Mean Square Errors
  
  lm.ZSnull.hosp <- bas.lm(hospH8 ~ gialla + arancione + rossa + newpos + intcar + hosp + newpos_av7D,
                           training.set.cat,
                           include.always = ~ .,
                           n.models = 1)
  mse.ZSnull.hosp[i] <- mse(lm.ZSnull.hosp, validation.set.cat, "hospH8")
  
  lm.ZSnull.intcar <- bas.lm(intcarH8 ~ gialla + arancione + rossa + newpos + intcar + hosp + newpos_av7D,
                             training.set.cat,
                             include.always = ~ .,
                             n.models = 1)
  mse.ZSnull.intcar[i] <- mse(lm.ZSnull.intcar, validation.set.cat, "intcarH8")
  
  lm.gprior.hosp <- bas.lm(hospH8 ~ gialla + arancione + rossa + newpos + intcar + hosp + newpos_av7D,
                           training.set.cat,
                           prior = "g-prior",
                           alpha = 100,
                           include.always = ~ .,
                           n.models = 1)
  mse.g.hosp[i] <- mse(lm.gprior.hosp, validation.set.cat, "hospH8")
  
  lm.gprior.intcar <- bas.lm(intcarH8 ~ gialla + arancione + rossa + newpos + intcar + hosp + newpos_av7D,
                             training.set.cat,
                             prior = "g-prior",
                             alpha = 100,
                             include.always = ~ .,
                             n.models = 1)
  mse.g.intcar[i] <- mse(lm.gprior.intcar, validation.set.cat, "intcarH8")
  
  lm.bic.hosp <- bas.lm(hospH8 ~ gialla + arancione + rossa + newpos + intcar + hosp + newpos_av7D,
                        training.set.cat,
                        prior = "BIC")
  mse.bic.bma.hosp[i] <- mse(lm.bic.hosp, validation.set.cat, "hospH8", "BMA")
  mse.bic.hpm.hosp[i] <- mse(lm.bic.hosp, validation.set.cat, "hospH8", "HPM")
  
  lm.bic.intcar <- bas.lm(intcarH8 ~ gialla + arancione + rossa + newpos + intcar + hosp + newpos_av7D,
                          training.set.cat,
                          prior = "BIC")
  mse.bic.bma.intcar[i] <- mse(lm.bic.intcar, validation.set.cat, "intcarH8", "BMA")
  mse.bic.hpm.intcar[i] <- mse(lm.bic.intcar, validation.set.cat, "intcarH8", "HPM")
  
  lm.numcolor.hosp <- bas.lm(hospH8 ~ color.num + newpos + intcar + hosp + newpos_av7D,
                             training.set.num,
                             prior = "BIC")
  mse.numcolor.hosp[i] <- mse(lm.numcolor.hosp, validation.set.num, "hospH8")
  
  lm.numcolor.intcar <- bas.lm(intcarH8 ~ color.num + newpos + intcar + hosp + newpos_av7D,
                               training.set.num,
                               prior = "BIC")
  mse.numcolor.intcar[i] <- mse(lm.numcolor.intcar, validation.set.num, "intcarH8")
  
  lm.szn.hosp <- bas.lm(hospH8 ~ gialla + arancione + rossa + newpos + intcar + hosp + newpos_av7D + winter + summer,
                        training.set.cat.szn,
                        prior = "BIC")
  mse.szn.hosp[i] <- mse(lm.szn.hosp, validation.set.cat.szn, "hospH8")
  
  lm.szn.intcar <- bas.lm(intcarH8 ~ gialla + arancione + rossa + newpos + intcar + hosp + newpos_av7D + winter + summer,
                          training.set.cat.szn,
                          prior = "BIC")
  mse.szn.intcar[i] <- mse(lm.szn.intcar, validation.set.cat.szn, "intcarH8")
}

mean(mse.ZSnull.hosp)
mean(mse.g.hosp)
mean(mse.bic.bma.hosp)
mean(mse.bic.hpm.hosp)
mean(mse.numcolor.hosp)
mean(mse.szn.hosp)

mean(mse.ZSnull.intcar)
mean(mse.g.intcar)
mean(mse.bic.bma.intcar)
mean(mse.bic.hpm.intcar)
mean(mse.numcolor.intcar)
mean(mse.szn.intcar)
