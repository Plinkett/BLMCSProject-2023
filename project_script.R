library(ggplot2)
library(ggpubr)
library(ggcorrplot)
library(GGally)
library(latex2exp)
library(BAS)

setwd('/home/ridley/Documents/GitHub/BLMCSProject-2023')
# setwd('/home/ridley/Documents/R scripts/BLMCS')

data <- read.csv('covidLom2020_21.csv')
View(data)

# ggplot <- function(...) ggplot2::ggplot(...) + scale_color_manual(values=c('#FF9933', 'gray', '#d8ce27', "#C41E3A")) + scale_fill_manual(values=c('#FF9933', 'gray', '#d8ce27', "#C41E3A"))
# unlockBinding("ggplot",parent.env(asNamespace("GGally")))
# assign("ggplot",ggplot,parent.env(asNamespace("GGally")))

# Scatter plot for covariates WITHOUT day of the week
X11()
p1 <- ggpairs(data,
        columns = c('newpos', 'intcar', 'hosp'),
        columnLabels = c('NewPositiveSubjects', 'PatientsICU', 'PatientsAtHospital'),
        aes(color = color,
               alpha= 0.5)) + scale_color_manual(values=c('#FF9933', 'gray', '#d8ce27', "#C41E3A")) + scale_fill_manual(values=c('#FF9933', 'gray', '#d8ce27', "#C41E3A"))
p1

# Same scatterplot as before but with one of the target variables to see correlation (relationship seems linear "enough")
X11()
p2 <- ggpairs(data,
        columns = c('newpos', 'intcar', 'hosp', 'newpos_av7D','hospH8'),
        columnLabels = c('NewPositiveSubjects', 'PatientsICU', 'PatientsAtHospital', 'NewPositiveAverage','Patients7DayAhead'),
        aes(color = color,
            alpha= 0.5)) + scale_color_manual(values=c('#FF9933', 'gray', '#d8ce27', "#C41E3A")) + scale_fill_manual(values=c('#FF9933', 'gray', '#d8ce27', "#C41E3A"))
p2

# Same as before, but with intensive care instead of patients (after 7 days)
X11()
p3 <- ggpairs(data,
              columns = c('newpos', 'intcar', 'hosp', 'newpos_av7D','intcarH8'),
              columnLabels = c('NewPositiveSubjects', 'PatientsICU', 'PatientsAtHospital', 'NewPositiveAverage','PatientsICU7DayAhead'),
              aes(color = color,
                  alpha= 0.5)) + scale_color_manual(values=c('#FF9933', 'gray', '#d8ce27', "#C41E3A")) + scale_fill_manual(values=c('#FF9933', 'gray', '#d8ce27', "#C41E3A"))
p3

# Boxplot for each numeric column (mean centered)
X11()
df_numeric <- data[, c('newpos', 'intcar', 'hosp', 'newpos_av7D')]
df_numeric <- as.data.frame(scale(df_numeric, scale=FALSE))

boxp <- ggplot(stack(df_numeric), 
            aes(x = ind, y = values)) + geom_boxplot(fill=terrain.colors(4))
boxp
# Since the variance between the features is too large we may want to normalize
# the data before   


# Boxplot with normalized columns
X11()
df_scaled <- as.data.frame(scale(df_numeric, scale=TRUE))
boxp_scaled <- ggplot(stack(df_scaled), 
               aes(x = ind, y = values)) + geom_boxplot(fill=terrain.colors(4))
#boxp_scaled <- boxp_scaled + scale_fill_brewer(palette="BuPu")
boxp_scaled

# Now let's look at the correlation plots (heatmaps) without considering the group
# distinction

corrmat <- cor(df_scaled)
ggcorrplot(corrmat, lab = TRUE, outline.col = 'white')

# Now let's separate them by class

colors <- unique(data$color)

plot_list <- list()
i <- 1
for(color in colors) {
  color_indices <- which(data$color == color)
  corrmat <- cor(df_scaled[color_indices,])
  corrmat[is.na(corrmat)] <- 0
  plottemp <- ggcorrplot(corrmat, lab = TRUE, outline.col = 'white')
  plot_list[[i]] <- plottemp
  i <- i + 1
}
X11()
ggarrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], 
          plot_list[[4]], labels = c('Arancione', 'Gialla', 'Rossa', 'Bianca'),
          ncol = 2, nrow = 2)

# I guess we will have two models, one for the hospital ordinary patients after 7 days
# and another one for the intensive care patients after 7 days

# Let's try some random stuff with intensive care!

output <- data.frame(x=numeric(0),y_hat=numeric(0),alpha=numeric(0))
X <- df_scaled$intcar
new_x <- seq(min(X), max(X), length.out = 100)
for(alpha in c(0.001, 0.01, 0.1, 1, 1.0)) {
  blm <- bas.lm(df_scaled$hospH8 ~ X, prior = 'g-prior', alpha = alpha,
                     modelprior = Bernoulli(1))
  intercept <- coef(blm)$postmean[1]
  beta <- coef(blm)$postmean[2]
  
  y_hat <- intercept + beta * new_x
  output <- rbind(output, data.frame(x = new_x, y_hat = y_hat, alpha = alpha))
}

blm <- bas.lm(df_scaled$hospH8 ~ X, prior = 'BIC', 
              modelprior = Bernoulli(1))
intercept <- coef(blm)$postmean[1]
beta <- coef(blm)$postmean[2]

y_hat <- intercept + beta * new_x
output = rbind(output, data.frame(x = new_x, y_hat = y_hat, 
                                  alpha="Noninformative prior"))
output$alpha <- as.factor(output$alpha)

plot_regression <- ggplot(data = df_scaled, aes(x = intcar, y = hospH8))
plot_regression <- plot_regression + geom_point(color = 'steelblue') + xlab('Patients in IC') +
                    ylab('Hospital patients after 8 days') +
                    geom_line(data = output, aes(x = x, y = y_hat, group = alpha, color = alpha), lty = 1) +
                    scale_colour_manual(values = rainbow(length(levels(output$alpha))), name=TeX('$\\alpha$'))
plot_regression


#################
#################
## REGRESSION  ##
#################
#################


# We skip the first prior and use Zellner's prior









