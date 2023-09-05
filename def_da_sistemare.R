# Bayesian Learning and Montecarlo Simulation Project

# Group 14: Guerrini Michele                  10607256
#           Mozzi Davide                      10674709
#           Santillan Moreno Carlos Alberto   10659783

# Dataset 2.8: Covid Data (covidLom2020_21.csv)


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
  y <- predict(model, test.data, estimator = estimator)
  if (hasName(y, "fit")) y <- y$fit
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
n         <- nrow(covid.data)
n.train   <- 180

# Shuffle the data
shuffled.data <- shuffle(covid.data)

# Create training and test sets
training.set <- shuffled.data[seq(1, n.train),]
test.set <- shuffled.data[seq(n.train + 1, n),]

# Standardize the data
# The targets are not standardized. The test set is standardized using the same center and scale values used for the
# training set.
to.standardize <- c("newpos", "intcar", "hosp", "newpos_av7D")
training.set <- standardize(training.set, to.standardize)

centers <- list()
scales <- list()
for (name in to.standardize) {
  centers[[name]] <- attr(training.set[[name]], "scaled:center")
  scales[[name]] <- attr(training.set[[name]], "scaled:scale")
}

test.set <- standardize(test.set, to.standardize, centers, scales)

# Transform the "color" categorical variable in several indicator variables
training.set.cat <- color.to.indicators(training.set)
test.set.cat <- color.to.indicators(test.set)

# Transform the "color" categorical variable in a numerical variable, assuming that
# "Bianca" < "Gialla" < "Arancione" < "Rossa" and that the value increases linearly.
training.set.num <- color.to.numerical(training.set)
test.set.num <- color.to.numerical(test.set)


#### Testing different priors ####
# Note: here we are training only the full model each time (the one with all the covariates).

## Default prior of the BAS library ("ZS-null"), which is an approximation of the Zellner-Siow prior

# Target: hospH8

lm.ZSnull.hosp <- bas.lm(hospH8 ~ gialla + arancione + rossa + newpos + intcar + hosp + newpos_av7D,
                    training.set.cat,
                    include.always = ~ .,
                    n.models = 1)
summary(lm.ZSnull.hosp)

mse(lm.ZSnull.hosp, test.set.cat, "hospH8")

coef(lm.ZSnull.hosp)

par(mfrow = c(2, 4))
plot(coef(lm.ZSnull.hosp), ask = FALSE) # TODO: subset = 2:8 ?
par(mfrow = c(1, 1))

plot(confint(coef(lm.ZSnull.hosp)))

# Target: intcarH8

lm.ZSnull.intcar <- bas.lm(intcarH8 ~ gialla + arancione + rossa + newpos + intcar + hosp + newpos_av7D,
                           training.set.cat,
                           include.always = ~ .,
                           n.models = 1)
summary(lm.ZSnull.intcar)

mse(lm.ZSnull.intcar, test.set.cat, "intcarH8")

coef(lm.ZSnull.intcar)

par(mfrow = c(2, 4))
plot(coef(lm.ZSnull.intcar), ask = FALSE)
par(mfrow = c(1, 1))

plot(confint(coef(lm.ZSnull.intcar)))


## Zellner's informative g-prior

# Target: hospH8

lm.gprior.hosp <- bas.lm(hospH8 ~ gialla + arancione + rossa + newpos + intcar + hosp + newpos_av7D,
                         training.set.cat,
                         prior = "g-prior",
                         alpha = 100,
                         include.always = ~ .,
                         n.models = 1)
summary(lm.gprior.hosp)

mse(lm.gprior.hosp, test.set.cat, "hospH8")

coef(lm.gprior.hosp)

par(mfrow = c(2, 4))
plot(coef(lm.gprior.hosp), ask = FALSE) # TODO: subset = 2:8 ?
par(mfrow = c(1, 1))

plot(confint(coef(lm.gprior.hosp)))

# Target: intcarH8

lm.gprior.intcar <- bas.lm(intcarH8 ~ gialla + arancione + rossa + newpos + intcar + hosp + newpos_av7D,
                           training.set.cat,
                           prior = "g-prior",
                           alpha = 100,
                           include.always = ~ .,
                           n.models = 1)
summary(lm.gprior.intcar)

mse(lm.gprior.intcar, test.set.cat, "intcarH8")

coef(lm.gprior.intcar)

par(mfrow = c(2, 4))
plot(coef(lm.gprior.intcar), ask = FALSE)
par(mfrow = c(1, 1))

plot(confint(coef(lm.gprior.intcar)))


#### Model selection using BIC

# Target: hospH8

lm.bic.hosp <- bas.lm(hospH8 ~ gialla + arancione + rossa + newpos + intcar + hosp + newpos_av7D,
                      training.set.cat,
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

par(mfrow = c(2, 4))
plot(coef(lm.bic.hosp), ask = FALSE)
par(mfrow = c(1, 1))

mse(lm.bic.hosp, test.set.cat, "hospH8", "BMA")
mse(lm.bic.hosp, test.set.cat, "hospH8", "HPM")

lm.best.hosp <- bas.lm(hospH8 ~ gialla + arancione + rossa + newpos + intcar + hosp + newpos_av7D,
                      training.set.cat,
                      prior = "BIC")

# Target: intcarH8

lm.bic.intcar <- bas.lm(intcarH8 ~ gialla + arancione + rossa + newpos + intcar + hosp + newpos_av7D,
                      training.set.cat,
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

par(mfrow = c(2, 4))
plot(coef(lm.bic.intcar), ask = FALSE)
par(mfrow = c(1, 1))

mse(lm.bic.intcar, test.set.cat, "intcarH8", "BMA")
mse(lm.bic.intcar, test.set.cat, "intcarH8", "HPM")


#### Appendix: different representations for some of the covariates ####

## Color as a numerical variable

# TODO

## Season categorical variable (derived from day)

# TODO
