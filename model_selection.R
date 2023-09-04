library(BAS)

# TODO:
# [x] Add models to predict intcarH8
# [x] Use the parameters of the training set to "standardize" the test set
# [x] Do not standardize targets
# [ ] Try "splitting" some covariates and combining them with others
#     (e.g. intcar.bianca, intcar.gialla, ...)
# [ ] Try different priors
# [ ] Sensitivity analysis on prior parameters

#### Functions ####

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

mse <- function (model, test.data, target.var) {
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
# The targets are not standardized. The test set is standardized using the same
# center and scale values used for the training set.
to.standardize <- c("newpos", "intcar", "hosp", "newpos_av7D")
training.set <- standardize(training.set, to.standardize)

centers <- list()
scales <- list()
for (name in to.standardize) {
  centers[[name]] <- attr(training.set[[name]], "scaled:center")
  scales[[name]] <- attr(training.set[[name]], "scaled:scale")
}

test.set <- standardize(test.set, to.standardize, centers, scales)


#### Model 1 ####
# "color" is not used
# "day" is not used
# "newpos_av7D" is not used
# "hospH8" is the only target variable
# This is mainly to experiment with the training/test sets to see why they have
# to be taken from the shuffled data (this is done using non-bayesian linear
# models). Then there is also a first test with a bayesian model to compare it
# with the frequentist one.

# Note: here we are using the non-standardized data. Model 1a uses a portion of
# the unshuffled data set as the training set, which is not good, as we'll see.

training.interval <- seq(1, n.train)
# training.interval <- 1:180
# training.interval <- 25:205
# training.interval <- 1:100
test.interval <- seq(n.train + 1, n)
# test.interval <- 181:205
# test.interval <- 1:24
# test.interval <- 101:205

# OLS linear regression
covid.lm.1a <- lm(hospH8 ~ newpos + intcar + hosp, covid.data, training.interval)
covid.lm.1b <- lm(hospH8 ~ newpos + intcar + hosp, shuffled.data, training.interval)
summary(covid.lm.1a)
summary(covid.lm.1b)

# Testing linear model 1a (not bayesian), computing the Mean Square Error on the
# test set
mse(covid.lm.1a, covid.data[test.interval,], "hospH8")

# Testing linear model 1b (not bayesian)
mse(covid.lm.1b, shuffled.data[test.interval,], "hospH8")

# Comparing the models 1a and 1b, changing the size of the training set and
# the portion of the data that it is made of, we can see that the average error
# of model 1a is extremely variable compared with the one of 1b (shuffled data).
#
# Examples (note that the results of 1b depend on how the data is shuffled):
#
# With training set 1:180, test set 181:205
# -> MSE of 1a = 439.1033
# -> MSE of 1b = 2151.441
#
# With training set 25:205, test set 1:24
# -> MSE of 1a = 3753.199
# -> MSE of 1b = 1221.284
#
# With training set 1:100, test set 101:205
# -> MSE of 1a = 29620.29
# -> MSE of 1b = 1756.088
#
# This is due to the fact that we are ignoring some covariates (and possibly
# some other correlations that our data doesn't fully capture), so choosing as
# training set a specific period of time makes it highly biased.
# So, in order to better evaluate the models, from this point forward, we'll
# always use a shuffled data set, as prepared in the "Data acquisition" part.

# Bayesian linear regression
covid.lm.1c <- bas.lm(hospH8 ~ newpos + intcar + hosp,
                      shuffled.data,
                      training.interval,
                      # prior = "g-prior",
                      # alpha = 10,
                      n.models = 1,
                      include.always = ~ .)
summary(covid.lm.1c)

# Testing linear model 1c
mse(covid.lm.1c, shuffled.data[test.interval,], "hospH8")


#### Model 2 ####
# "color" is treated as a categorical variable first, and as a numerical value
# (implying an ordering between the colors) later
# "day" is not used
# "newpos_av7D" is not used
# Here we are still forcing 'bas.lm' to only train one model. The goal is to
# compare the two ways of representing the "color" variable

# Transform the "color" categorical variable in several indicator variables
training.set.cat <- color.to.indicators(training.set)
test.set.cat <- color.to.indicators(test.set)

# Train and test model 2a (one model for "hospH8" and one for "intcarH8")
# Note: we can leave out one of the indicator variables (e.g. 'bianca') to get
# full rank conditions.
covid.lm.2a.hosp <- bas.lm(hospH8 ~ gialla + arancione + rossa + newpos +
                        intcar + hosp,
                      training.set.cat,
                      n.models = 1,
                      include.always = ~ .)
summary(covid.lm.2a.hosp)

mse(covid.lm.2a.hosp, test.set.cat, "hospH8")

covid.lm.2a.intcar <- bas.lm(intcarH8 ~ gialla + arancione + rossa + newpos +
                             intcar + hosp,
                           training.set.cat,
                           n.models = 1,
                           include.always = ~ .)
summary(covid.lm.2a.intcar)

mse(covid.lm.2a.intcar, test.set.cat, "intcarH8")

# Transform the "color" categorical variable in a numerical variable, assuming
# that "Bianca" < "Gialla" < "Arancione" < "Rossa" and that the value increases
# linearly (the second assumption is not really justified).
training.set.num <- color.to.numerical(training.set)
test.set.num <- color.to.numerical(test.set)

# Train and test model 2b

covid.lm.2b.hosp <- bas.lm(hospH8 ~ color.num + newpos + intcar + hosp,
                      training.set.num,
                      n.models = 1,
                      include.always = ~ .)
summary(covid.lm.2b.hosp)

mse(covid.lm.2b.hosp, test.set.num, "hospH8")

covid.lm.2b.intcar <- bas.lm(intcarH8 ~ color.num + newpos + intcar + hosp,
                      training.set.num,
                      n.models = 1,
                      include.always = ~ .)
summary(covid.lm.2b.intcar)

mse(covid.lm.2b.intcar, test.set.num, "intcarH8")

# It seems that model 2a (i.e. with "color" as a categorical variable, without
# ordering) consistently performs slightly better than model 2b.


#### Model 3 ####
# "color" is treated as a categorical variable
# "day" is not used
# Here we let 'bas.lm' test out different models, potentially leaving out some
# of the covariates

covid.lm.3.hosp <- bas.lm(hospH8 ~ gialla + arancione + rossa + newpos +
                       intcar + hosp + newpos_av7D,
                     training.set.cat)
summary(covid.lm.3.hosp)

mse(covid.lm.3.hosp, test.set.cat, "hospH8")

# Show which features are included in the top models selected by 'bas.lm'
plot(seq(coef(covid.lm.3.hosp)$probne0), coef(covid.lm.3.hosp)$probne0,
     type = "h",
     lwd = 4,
     xaxt = "n",
     main = "Feature inclusion probabilities",
     xlab = "",
     ylab = "post p(B != 0)")
axis(1, seq(coef(covid.lm.3.hosp)$probne0), labels = coef(covid.lm.3.hosp)$namesx)

plot(coef(covid.lm.3.hosp), subset = 2:8)

print(covid.lm.3.hosp)

image(covid.lm.3.hosp, rotate = F)

# intcarH8

covid.lm.3.intcar <- bas.lm(intcarH8 ~ gialla + arancione + rossa + newpos +
                            intcar + hosp + newpos_av7D,
                          training.set.cat)
summary(covid.lm.3.intcar)

mse(covid.lm.3.intcar, test.set.cat, "intcarH8")

# Show which features are included in the top models selected by 'bas.lm'
plot(seq(coef(covid.lm.3.intcar)$probne0), coef(covid.lm.3.intcar)$probne0,
     type = "h",
     lwd = 4,
     xaxt = "n",
     main = "Feature inclusion probabilities",
     xlab = "",
     ylab = "post p(B != 0)")
axis(1, seq(coef(covid.lm.3.intcar)$probne0), labels = coef(covid.lm.3.intcar)$namesx)

plot(coef(covid.lm.3.intcar), subset = 2:8)

print(covid.lm.3.intcar)

image(covid.lm.3.intcar, rotate = F)

# The model for intcarH8 cares more about Arancione probably because the number of people in ICU was fundamental in the
# decision of the color of the zone. The whole point of the Arancione and Rossa zones was to save hospitals from
# being overcrowded and only secondly to diminish the number of infected people.

#### Model 4 ####
# "color" is treated as a categorical variable
# "day" is converted to a "season" categorical variable (and therefore into four
# indicator variables). We are considering the meteorological seasons for
# simplicity (i.e. spring={mar, apr, may}, summer={jun, jul, aug}, ...)

# Prepare the data
training.set.cat.szn <- day.to.season(training.set.cat)
test.set.cat.szn <- day.to.season(test.set.cat)

# Note: we can for sure ignore the "autumn" indicator variable since in our data
# it's always 0. Then we can remove one of the remaining three in order to get
# independent conditions.
covid.lm.4.hosp <- bas.lm(hospH8 ~ gialla + arancione + rossa + newpos +
                       intcar + hosp + newpos_av7D + winter + summer,
                     training.set.cat.szn)
summary(covid.lm.4.hosp)

mse(covid.lm.4.hosp, test.set.cat.szn, "hospH8")

plot(seq(coef(covid.lm.4.hosp)$probne0), coef(covid.lm.4.hosp)$probne0,
     type = "h",
     lwd = 4,
     xaxt = "n",
     main = "Feature inclusion probabilities",
     xlab = "",
     ylab = "post p(B != 0)")
axis(1, seq(coef(covid.lm.4.hosp)$probne0), labels = coef(covid.lm.4.hosp)$namesx)

# With the data we have (which spans over less than one year) the season seems
# to be relevant, in particular whether it's winter or not.

# Note: I don't know why this...
predict(covid.lm.4.hosp, test.set.cat.szn[1,])$fit
# ... is different from this...
coef(covid.lm.4.hosp)[[1]] %*% as.matrix(cbind(1, test.set.cat.szn[1, c("gialla", "arancione", "rossa", "newpos", "intcar", "hosp", "newpos_av7D", "winter", "summer")]))[1,]
# For reference, the true value of "hospH8" for this data point is:
test.set.cat.szn[1,]$hospH8

# intcarH8
covid.lm.4.intcar <- bas.lm(intcarH8 ~ gialla + arancione + rossa + newpos +
                            intcar + hosp + newpos_av7D + winter + summer,
                          training.set.cat.szn)
summary(covid.lm.4.intcar)

mse(covid.lm.4.intcar, test.set.cat.szn, "intcarH8")

plot(seq(coef(covid.lm.4.intcar)$probne0), coef(covid.lm.4.intcar)$probne0,
     type = "h",
     lwd = 4,
     xaxt = "n",
     main = "Feature inclusion probabilities",
     xlab = "",
     ylab = "post p(B != 0)")
axis(1, seq(coef(covid.lm.4.intcar)$probne0), labels = coef(covid.lm.4.intcar)$namesx)
