library(ggplot2)
library(GGally)
setwd('/home/ridley/Documents/R scripts/BLMCS')

data <- read.csv('covidLom2020_21.csv')
View(data)

# ggplot <- function(...) ggplot2::ggplot(...) + scale_color_manual(values=c('#FF9933', 'gray', '#d8ce27', "#C41E3A")) + scale_fill_manual(values=c('#FF9933', 'gray', '#d8ce27', "#C41E3A"))
# unlockBinding("ggplot",parent.env(asNamespace("GGally")))
# assign("ggplot",ggplot,parent.env(asNamespace("GGally")))

# Scatter plot for covariates WITHOUT day of the week
X11()
p1 <- ggpairs(data,
        columns = c('newpos', 'intcar', 'hosp'),
        columnLabels = c('#NewPositives', '#PatientsIntensiveCare', '#PatientsHospital'),
        aes(color = color,
               alpha= 0.5)) + scale_color_manual(values=c('#FF9933', 'gray', '#d8ce27', "#C41E3A")) + scale_fill_manual(values=c('#FF9933', 'gray', '#d8ce27', "#C41E3A"))
p1

# Same scatterplot as before but with one of the target variables to see correlation (relationship seems linear "enough")
X11()
p2 <- ggpairs(data,
        columns = c('newpos', 'intcar', 'hosp', 'hospH8'),
        columnLabels = c('#NewPositives', '#PatientsIntensiveCare', '#PatientsHospital', '#PatientH8'),
        aes(color = color,
            alpha= 0.5)) + scale_color_manual(values=c('#FF9933', 'gray', '#d8ce27', "#C41E3A")) + scale_fill_manual(values=c('#FF9933', 'gray', '#d8ce27', "#C41E3A"))
p2

# Same as before, but with intensive care instead of patients (after 7 days)
X11()
p3 <- ggpairs(data,
        columns = c('newpos', 'intcar', 'hosp', 'intcarH8'),
        columnLabels = c('#NewPositives', '#PatientsIntensiveCare', '#PatientsHospital', '#IntensiveH8'),
        aes(color = color,
            alpha= 0.5)) + scale_color_manual(values=c('#FF9933', 'gray', '#d8ce27', "#C41E3A")) + scale_fill_manual(values=c('#FF9933', 'gray', '#d8ce27', "#C41E3A"))
p3

# Boxplot for each numeric column (mean centered)
X11()
df_numeric <- data[, c('newpos', 'intcar', 'hosp', 'newpos_av7D', 'hospH8', 'intcarH8')]
df_numeric <- as.data.frame(scale(df_numeric, scale=FALSE))

boxp <- ggplot(stack(df_numeric), 
            aes(x = ind, y = values)) + geom_boxplot(fill=terrain.colors(6))
boxp

# Boxplot with normalized columns
X11()
df_scaled <- as.data.frame(scale(df_numeric, scale=TRUE))
boxp_scaled <- ggplot(stack(df_scaled), 
               aes(x = ind, y = values)) + geom_boxplot(fill=terrain.colors(6))
#boxp_scaled <- boxp_scaled + scale_fill_brewer(palette="BuPu")
boxp_scaled



