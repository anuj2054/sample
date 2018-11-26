
library(dplyr)
library(ggplot2)
library(ztable) # for the colors on the contingency table
library(lme4) # for the mixed effect
library(broom) # to tidy the data

##################
## DATA CLEANING
####################

# read in the data
efficacy <- read.csv("efficacy.csv")
randomization <- read.csv("randomization.csv")
subject <- read.csv("subject.csv")

# joining the datasets
complete <- dplyr::full_join(efficacy,randomization)
complete <- dplyr::full_join(complete,subject)

# removing the mucus viscorit that is NA
complete <- complete %>% filter(!(is.na(mucus.viscosity)))

## some transformations to the rate
complete$difference <- (complete$previous.year - complete$nosebleeds)
complete$rate <- 365* ( complete$difference / complete$duration ) 
complete$transformed_rate <- sign(((complete$previous.year - complete$nosebleeds)/complete$duration)*365+3) * abs(((complete$previous.year - complete$nosebleeds)/complete$duration)*365+3)^(1/3)
complete$transformed_rate <- log(((complete$previous.year - complete$nosebleeds)/complete$duration)*365+4.5)

## removing the outlier based on the visuaizations
complete <- complete %>% filter(!(rate == 146))

# some subsets of the data
active_only <- complete %>% filter(arm == "ACTIVE") 
placebo_only <- complete %>% filter(arm == "PLACEBO") 

#######################################################################3
# EXPLORING DISTRIBUTIONS OF RESPONSE DATA and DATA TRANSFORMATIONS
####################################################################33


# creating a dervied variable called rate
complete$rate <- ((complete$previous.year - complete$nosebleeds)/complete$duration)*365
## finding the data distribution of the rate
shapiro.test(complete$rate)
max(complete$rate)
min(complete$rate)
mean(complete$rate)
sd(complete$rate)
qqnorm(complete$rate)
hist(complete$rate)
plot(density(complete$rate, na.rm = TRUE))
ggplot(complete, aes(x = rate)) + geom_density()
plot(density(complete$rate))
ggplot(complete, aes(sample = rate, colour = arm)) +
  stat_qq() +
  stat_qq_line()

# finding the data distribution of the nosebleeds
qqnorm(complete$nosebleeds)
shapiro.test(complete$nosebleeds)
hist(complete$nosebleeds)
max(complete$nosebleeds)
min(complete$nosebleeds)
mean(complete$nosebleeds)
sd(complete$nosebleeds)
plot(density(complete$nosebleeds, na.rm = TRUE))


# transforming rate using cube root
complete$transformed_rate <- sign(((complete$previous.year - complete$nosebleeds)/complete$duration)*365+3) * abs(((complete$previous.year - complete$nosebleeds)/complete$duration)*365+3)^(1/3)
# finding the data distributoin of the transformed variable
plot(density(complete$transformed_rate))
ggplot(complete, aes(sample = transformed_rate, colour = arm)) +
  stat_qq() +
  stat_qq_line()
shapiro.test(complete$transformed_rate)
shapiro.test(complete$transformed_rate[complete$arm=="ACTIVE"])

# transforming rate using log transform
complete$transformed_rate <- log(((complete$previous.year - complete$nosebleeds)/complete$duration)*365+4.5)
plot(density(complete$transformed_rate))
ggplot(complete, aes(sample = transformed_rate, colour = arm)) +
  stat_qq() +
  stat_qq_line()
shapiro.test(complete$transformed_rate)
shapiro.test(complete$transformed_rate[complete$arm=="ACTIVE"])


# in case the transformed variable still does seem right, finding the parameters of the gamma distribution to use non-parametric statistics
# shape of 10.019 and rate of 1.45 was found for the rate
library(fitdistrplus)
fit.gamma <- fitdist(complete$rate, distr = "gamma", method = "mle")
plot(fit.gamma)
summary(fit.gamma)
# finding how well my data fits to the gamma distribution found using above procedure
qqplot(complete$rate,rnorm(300,1.879437,0.2880162))
ggplot(complete, aes(sample = rate, colour = arm)) +
  stat_qq() +
  stat_qq_line()


