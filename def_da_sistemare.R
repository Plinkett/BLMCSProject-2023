# Bayesian Learning and Montecarlo Simulation Project

# Group 14: Guerrini Michele                  10607256
#           Mozzi David"e                      10674709
#           Santillan Moreno Carlos Alberto   10659783

# Dataset 2.8: Covid Data (covidLom2020_21.csv)

#### Functions ####

# List of useful functions used in the project

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

# Transform the "color" categorical variable in several indicator variables
training.set.cat <- color.to.indicators(training.set)
test.set.cat <- color.to.indicators(test.set)


#### Model 1 ####

# First model for the target hospH8
covid.1a <- bas.lm(hospH8 ~ gialla + arancione + rossa + newpos +
                           intcar + hosp + newpos_av7D,
                         data = training.set.cat,
                         modelprior = uniform())
summary(covid.1a)

test(covid.1a, test.set.cat, "hospH8")

# Show the covariates inclusion probabilities, their values and the marginal posterior distribution of every covariate
print(covid.1a)
par(mfrow=c(2,4))
plot(seq(coef(covid.1a)$probne0), coef(covid.1a)$probne0,
     type = "h",
     lwd = 4,
     xaxt = "n",
     main = "Covariates inclusion probabilities",
     xlab = "",
     ylab = "post p(B != 0)")
axis(1, seq(coef(covid.1a)$probne0), labels = coef(covid.1a)$namesx)

plot(coef(covid.1a), subset = 2:8)

# Show which model includes which covariate 
par(mfrow=c(1,1))
image(covid.1a, rotate = F)


#Same model for intcarH8

covid.1b <- bas.lm(intcarH8 ~ gialla + arancione + rossa + newpos +
                             intcar + hosp + newpos_av7D,
                            data = training.set.cat,
                            modelprior = uniform())
summary(covid.1b)

test(covid.1b, test.set.cat, "intcarH8")

# Show the covariates inclusion probabilities, their values and the marginal posterior distribution of every covariate
print(covid.1a)
par(mfrow=c(2,4))
plot(seq(coef(covid.1b)$probne0), coef(covid.1b)$probne0,
     type = "h",
     lwd = 4,
     xaxt = "n",
     main = "Covariate inclusion probabilities",
     xlab = "",
     ylab = "post p(B != 0)")
axis(1, seq(coef(covid.1b)$probne0), labels = coef(covid.1b)$namesx)

plot(coef(covid.1a), subset = 2:8)

# Show which model includes which covariate 
par(mfrow=c(1,1))
image(covid.1b, rotate = F)

# The model for intcarH8 cares more about Arancione probably because the number of people in ICU was fundamental in the
# decision of the color of the zone. The whole point of the Arancione and Rossa zones was to save hospitals from
# being overcrowded and only secondly to diminish the number of infeted people.
