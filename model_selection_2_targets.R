library(BAS)

# TODO:
# - Add models to predict intcarH8
# - Use the parameters of the training set to "standardize" the test set
# - Do not standardize targets
# - Try "splitting" some covariates and combining them with others
#   (e.g. intcar.bianca, intcar.gialla, ...)
# - Try different priors
# - Sensitivity analysis on prior parameters

#### Functions ####

standardize <- function (data) {
  new.data <- data
  for (name in names(data)) {
    if (is.numeric(data[1, name])) {
      new.data[[name]] <- scale(data[[name]])
    }
  }
  return (new.data)
}

shuffle <- function (data) {
  return (data[sample(seq(1, nrow(data)), nrow(data)),])
}

test <- function (model, test.data, target.var) {
  y <- predict(model, test.data)
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

# Set useful constants and parameters
n         <- nrow(covid.data)
n.train   <- 180

# Standardize and shuffle the data
standardized.data <- standardize(covid.data)
prepared.data <- shuffle(standardized.data)

# Create training and test sets
training.set <- prepared.data[seq(1, n.train),]
test.set <- prepared.data[seq(n.train + 1, n),]


#### Model 1 ####
# "color" is not used
# "day" is not used
# "newpos_av7D" is not used
# "hospH8" is the only target variable
# This is mainly to experiment with the training/test sets and to compare the
# frequentist and bayesian approaches

# Note: model 1a uses a portion of the unshuffled data set as the training set,
# which is not good, as we'll see.

# OLS linear regression
covid.lm.1a <- lm(hospH8 ~ newpos + intcar + hosp, standardized.data, 1:n.train)
covid.lm.1b <- lm(hospH8 ~ newpos + intcar + hosp, training.set)
summary(covid.lm.1a)
summary(covid.lm.1b)

# Testing linear model 1a (not bayesian)
test(covid.lm.1a, standardized.data[seq(n.train + 1, n),], "hospH8")

# Testing linear model 1b (not bayesian)
test(covid.lm.1b, test.set, "hospH8")

# Comparing the models 1a and 1b, changing the size of the training set and
# the portion of the data that it is made of, we can see that the average error
# of model 1a is extremely variable compared with the one of 1b (shuffled data).
#
# Examples (note that the results of 1b depend on how the data is shuffled):
#
# With training set 1:180, test set 181:205
# -> avg error of 1a = 0.007546406
# -> avg error of 1b = 0.02546467
#
# With training set 25:205, test set 1:24
# -> avg error of 1a = 0.06450227
# -> avg error of 1b = 0.03603182
#
# With training set 1:100, test set 101:205
# -> avg error of 1a = 0.5090527
# -> avg error of 1b = 0.03330578
#
# This is due to the fact that we are ignoring some covariates, particularly the
# 'color', which has a relevant relationship with the target variables.
# Anyway, in order to better evaluate the models, from this point forward, we'll
# always use a shuffled data set.

# Bayesian linear regression
covid.lm.1c <- bas.lm(hospH8 ~ newpos + intcar + hosp,
                      training.set,
                      # prior = "g-prior",
                      # alpha = 10,
                      n.models = 1,
                      include.always = ~ .)
summary(covid.lm.1c)

# Testing linear model 1c
test(covid.lm.1c, training.set, "hospH8")

# Note: I don't know why this...
predict(covid.lm.1c, test.set[1,])$fit
# ... is different from this...
coef(covid.lm.1c)[[1]] %*% as.matrix(cbind(1, test.set[1, c("newpos", "intcar", "hosp")]))[1,]
# ... considering that we are only training one model, i.e. BMA should not be
# involved. Is it that BAS doesn't use the posterior mean of the coefficients?
# How to decide what estimator of the coefficient to use?


#### Model 2 ####
# "color" is treated as a categorical variable first, and as a numerical value
# (implying an ordering between the colors) later
# "day" is not used
# "newpos_av7D" is not used
# "hospH8" is the only target variable
# Here we are still forcing 'bas.lm' to only train one model. The goal is to
# compare the two ways of representing the "color" variable

# Transform the "color" categorical variable in several indicator variables
training.set.cat <- color.to.indicators(training.set)
test.set.cat <- color.to.indicators(test.set)

# Train and test model 2a
# Note: we can leave out one of the indicator variables (e.g. 'bianca') to get
# full rank conditions.
covid.lm.2aHosp <- bas.lm(hospH8 ~ gialla + arancione + rossa + newpos +
                        intcar + hosp,
                      training.set.cat,
                      n.models = 1,
                      include.always = ~ .)
summary(covid.lm.2aHosp)

test(covid.lm.2a, test.set.cat, "hospH8")


covid.lm.2aIntcar <- bas.lm(intcarH8 ~ gialla + arancione + rossa + newpos +
                        intcar + hosp,
                      training.set.cat,
                      n.models = 1,
                      include.always = ~ .)
summary(covid.lm.2aIntcar)

test(covid.lm.2aIntcar, test.set.cat, "intcarH8")

# Transform the "color" categorical variable in a numerical variable, assuming
# that "Bianca" < "Gialla" < "Arancione" < "Rossa" and that the value increases
# linearly (the second assumption is not really justified).
training.set.num <- color.to.numerical(training.set)
test.set.num <- color.to.numerical(test.set)

# Train and test model 2b

covid.lm.2bHosp <- bas.lm(hospH8 ~ color.num + newpos + intcar + hosp,
                      training.set.num,
                      n.models = 1,
                      include.always = ~ .)
summary(covid.lm.2bHosp)

test(covid.lm.2bHosp, test.set.num, "hospH8")


covid.lm.2bIntcar <- bas.lm(intcarH8 ~ color.num + newpos + intcar + hosp,
                      training.set.num,
                      n.models = 1,
                      include.always = ~ .)
summary(covid.lm.2b)

test(covid.lm.2bIntcar, test.set.num, "intcarH8")

# It seems that model 2a (i.e. with "color" as a categorical variable, without
# ordering) consistently performs slightly better than model 2b.


#### Model 3 ####
# "color" is treated as a categorical variable
# "day" is not used
# "hospH8" is the only target variable
# Here we let 'bas.lm' test out different models, potentially leaving out some
# of the covariates

covid.lm.3Hosp <- bas.lm(hospH8 ~ gialla + arancione + rossa + newpos +
                       intcar + hosp + newpos_av7D,
                     training.set.cat)
summary(covid.lm.3Hosp)

test(covid.lm.3Hosp, test.set.cat, "hospH8")

# Show which features are included in the top models selected by 'bas.lm'
plot(seq(coef(covid.lm.3Hosp)$probne0), coef(covid.lm.3Hosp)$probne0,
     type = "h",
     lwd = 4,
     xaxt = "n",
     main = "Feature inclusion probabilities",
     xlab = "",
     ylab = "post p(B != 0)")
axis(1, seq(coef(covid.lm.3Hosp)$probne0), labels = coef(covid.lm.3Hosp)$namesx)

plot(coef(covid.lm.3Hosp), subset = 2:8)

print(covid.lm.3Hosp)

image(covid.lm.3Hosp, rotate = F)

#intcarH8

covid.lm.3Intcar <- bas.lm(intcarH8 ~ gialla + arancione + rossa + newpos +
                       intcar + hosp + newpos_av7D,
                     training.set.cat)
summary(covid.lm.3Intcar)

test(covid.lm.3Intcar, test.set.cat, "intcarH8")

# plot(covid.lm.3)


plot(seq(coef(covid.lm.3Intcar)$probne0), coef(covid.lm.3Intcar)$probne0,
     type = "h",
     lwd = 4,
     xaxt = "n",
     main = "Feature inclusion probabilities",
     xlab = "",
     ylab = "post p(B != 0)")
axis(1, seq(coef(covid.lm.3Intcar)$probne0), labels = coef(covid.lm.3Intcar)$namesx)

# The model for intcarH8 cares more about Arancione probably because the number of people in ICU was fundamental in the
# decision of the color tìof the zone. The whole point of the Arancione and Rossa zones was to save hospitals from
# being overcrowded and only secondly to diminish the number of infeted people.

#### Model 4 ####
# "color" is treated as a categorical variable
# "day" is converted to a "season" categorical variable (and therefore into four
# indicator variables). We are considering the meteorological seasons for
# simplicity (i.e. spring={mar, apr, may}, summer={jun, jul, aug}, ...)
# "hospH8" is the only target variable

# Prepare the data
training.set.cat.szn <- day.to.season(training.set.cat)
test.set.cat.szn <- day.to.season(test.set.cat)

# Note: we can for sure ignore the "autumn" indicator variable since in our data
# it's always 0. Then we can remove one of the remaining three in order to get
# independent conditions.
covid.lm.4Hosp <- bas.lm(hospH8 ~ gialla + arancione + rossa + newpos +
                       intcar + hosp + newpos_av7D + winter + summer,
                     training.set.cat.szn)
summary(covid.lm.4Hosp)

test(covid.lm.4Hosp, test.set.cat.szn, "hospH8")

plot(seq(coef(covid.lm.4Hosp)$probne0), coef(covid.lm.4Hosp)$probne0,
     type = "h",
     lwd = 4,
     xaxt = "n",
     main = "Feature inclusion probabilities",
     xlab = "",
     ylab = "post p(B != 0)")
axis(1, seq(coef(covid.lm.4Hosp)$probne0), labels = coef(covid.lm.4Hosp)$namesx)


# incarH8
covid.lm.4Intcar <- bas.lm(intcarH8 ~ gialla + arancione + rossa + newpos +
                       intcar + hosp + newpos_av7D + winter + summer,
                     training.set.cat.szn)
summary(covid.lm.4Intcar)

test(covid.lm.4Intcar, test.set.cat.szn, "hospH8")

plot(seq(coef(covid.lm.4Intcar)$probne0), coef(covid.lm.4Intcar)$probne0,
     type = "h",
     lwd = 4,
     xaxt = "n",
     main = "Feature inclusion probabilities",
     xlab = "",
     ylab = "post p(B != 0)")
axis(1, seq(coef(covid.lm.4Intcar)$probne0), labels = coef(covid.lm.4Intcar)$namesx)

# With the data we have (which spans over less than one year) the season seems
# to be relevant, in particular whether it's winter or not.
