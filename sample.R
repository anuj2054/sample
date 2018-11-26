
#########################################
## Author Anuj Guruacharya
## Objective Completed as a work sample for a pharma company
## Date November 25 2018
#########################################

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



#########################################################################
##EXPLORATORY ANALYSIS ( descriptive stats, visualizations, cluster)
########################################################################


### figure 1 this is not needed
ggplot(complete, aes(x= previous.year, fill = nosebleeds)) +
  geom_histogram(binwidth=.5, position="dodge",stat = "count")

### figure 2 this is needed
complete.df <- as.data.frame(table(complete$previous.year,complete$nosebleeds))
ggplot(complete.df, aes(Var2, Var1)) +
  geom_tile(aes(fill = Freq), colour = "black") +
  scale_fill_gradient(low = "white", high = "red") + 
  geom_text(aes(label = round(Freq, 1))) +
  ggtitle("prevoius year vs nosebleed for complete data")

### figure 3 this is needed
complete.df <- as.data.frame(table(complete$previous.year,complete$tissue.use))
ggplot(complete.df, aes(Var2, Var1)) +
  geom_tile(aes(fill = Freq), colour = "black") +
  scale_fill_gradient(low = "white", high = "red") + 
  geom_text(aes(label = round(Freq, 1))) + 
  ggtitle("tissue use vs prevoius year nosebleed for complete data ")


## figure 4 this is not needed
ggplot(complete, aes(x=nosebleeds, y=rate, color = arm)) + 
  geom_boxplot() 

# figure 4 this is not needde
ggplot(complete, aes(x=nosebleeds, y=transformed_rate, color = arm)) + 
  geom_boxplot() 

# figure 4 this is needed
ggplot(complete, aes(x=nosebleeds, y=mucus.viscosity, color = arm)) + 
  geom_boxplot() 

# figure 4 this is not needed
ggplot(complete, aes(x=previous.year, y=rate, color = arm)) + 
  geom_boxplot() 

# figure 4 this is needed
ggplot(complete, aes(x=previous.year, y=transformed_rate, color = arm)) + 
  geom_boxplot() 

# figure 4 this is needed
ggplot(complete, aes(x=previous.year, y=mucus.viscosity)) + 
  geom_boxplot() 

# figure 4 this is needed
ggplot(complete, aes(x=arm, y=mucus.viscosity, color = arm)) + 
  geom_boxplot() + 
  geom_jitter(shape=16, position=position_jitter(0.2))

# figure 4 this is needed
ggplot(complete, aes(x=arm, y=transformed_rate, color = arm)) + 
  geom_boxplot() + 
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
  stat_summary(fun.y=mean, geom="point", shape=23, size=4)

# figure 4 this is needed
ggplot(complete, aes(x=country, y=transformed_rate, color = country)) + 
  geom_boxplot() 

# figure 4 this is needed
ggplot(complete, aes(x=country, y=transformed_rate, color = arm)) + 
  geom_boxplot() 

# figure 4 this is needed
ggplot(complete, aes(x=country, y=mucus.viscosity, color = country)) + 
  geom_boxplot()  

# figure 4 this is needed
ggplot(complete, aes(x=country, y=mucus.viscosity, color = arm)) + 
  geom_boxplot()  

# figure 4 this is needed
ggplot(complete, aes(x=eye.colour, y=mucus.viscosity, color = eye.colour)) + 
  geom_boxplot()

# figure 4 this is needed
ggplot(complete, aes(x=eye.colour, y=transformed_rate, color = arm)) + 
  geom_boxplot() 

# figure 4 this is not needed
ggplot(complete, aes(x=tissue.use, y=transformed_rate)) + 
  geom_boxplot() 

# figure 4 this is needed
ggplot(complete, aes(x=tissue.use, y=transformed_rate, color = arm)) + 
  geom_boxplot() 

# figure 4 this is needed
ggplot(complete, aes(x=tissue.use, y=mucus.viscosity)) + 
  geom_boxplot() 

# figure 4 this is needed
ggplot(complete, aes(x=tissue.use, y=mucus.viscosity, color = arm)) + 
  geom_boxplot() 


# figure 4a this is needed
complete.df <- as.data.frame(table(complete$country,complete$tissue.use))
ggplot(complete.df, aes(Var2, Var1)) +
  geom_tile(aes(fill = Freq), colour = "black") +
  scale_fill_gradient(low = "white", high = "red") + 
  geom_text(aes(label = round(Freq, 1))) 
# figure 4b this is not needed
ggplot(complete, aes(x= country, fill = tissue.use)) +
  geom_histogram(binwidth=.5, position="dodge",stat = "count")
# figure 4c this is not needed
ggplot(complete, aes(x= tissue.use, fill = country)) +
  geom_histogram(binwidth=.5, position="dodge",stat = "count")
# figure 4d this is not needed
ggplot(complete, aes(x= country, fill = tissue.use)) +
  geom_histogram(stat = "count")
# figure 4e this is not needed
ggplot(complete, aes(x= tissue.use, fill = country)) +
  geom_histogram(stat = "count")
# figure 4f this is not needed
ztable(table(complete$country,complete$tissue.use))  %>% makeHeatmap() %>% print(caption="Table 4. Heatmap Table")

# figure 4a this is not needed
ztable(table(complete$country,complete$arm))  %>% makeHeatmap() %>% print(caption="Table 4. Heatmap Table")
# figure 4b this is not needed
ggplot(complete, aes(x= country, fill = arm)) +
  geom_histogram(binwidth=.5, position="dodge",stat = "count")
# figure 4c this is not needed
complete.df <- as.data.frame(table(complete$country,complete$arm))
ggplot(complete.df, aes(Var2, Var1)) +
  geom_tile(aes(fill = Freq), colour = "black") +
  scale_fill_gradient(low = "white", high = "red") + 
  geom_text(aes(label = round(Freq, 1))) 

# figure 4a this is not needed
ztable(table(complete$country,complete$eye.colour))  %>% makeHeatmap() %>% print(caption="Table 4. Heatmap Table")
# figure 4a this is not needed
ggplot(complete, aes(x= country, fill = eye.colour)) +
  geom_histogram(binwidth = 0.5,bins = 3, position="dodge",stat = "count")
# figure 4c this is not needed
ggplot(complete, aes(x= country, fill = eye.colour)) +
  geom_histogram(stat = "count")
# figure 4c this is needed
complete.df <- as.data.frame(table(complete$country,complete$eye.colour))
ggplot(complete.df, aes(Var2, Var1)) +
  geom_tile(aes(fill = Freq), colour = "black") +
  scale_fill_gradient(low = "white", high = "red") + 
  geom_text(aes(label = round(Freq, 1))) 

# figure 4a this is not needd
ztable(table(complete$eye.colour,complete$tissue.use))  %>% makeHeatmap() %>% print(caption="Table 4. Heatmap Table")
# figure 4b this is not needed
complete.df <- as.data.frame(table(complete$tissue.use,complete$eye.colour))
ggplot(complete.df, aes(Var2, Var1)) +
  geom_tile(aes(fill = Freq), colour = "black") +
  scale_fill_gradient(low = "white", high = "red") + 
  geom_text(aes(label = round(Freq, 1))) 
# figure 4c this is needed
ggplot(complete %>% filter(!is.na(eye.colour)), aes(x= tissue.use, fill = eye.colour)) +
  geom_histogram(binwidth = 0.5,bins = 3, position="dodge",stat = "count")

# this is not needed
ztable(table(complete$tissue.use,complete$arm))  %>% makeHeatmap() %>% print(caption="Table 4. Heatmap Table")
# this is needed
complete.df <- as.data.frame(table(complete$arm,complete$tissue.use))
ggplot(complete.df, aes(Var2, Var1)) +
  geom_tile(aes(fill = Freq), colour = "black") +
  scale_fill_gradient(low = "white", high = "red") + 
  geom_text(aes(label = round(Freq, 1))) 



## this is needed
ggplot(data=complete, aes(x=mucus.viscosity, y=transformed_rate, group=arm)) +
  #  geom_line(aes(color=arm))+
  geom_point(aes(color=arm)) +
  ggtitle("viscosity vs transformed rate of complete data ")


## this is  needed
ggplot(data=complete, aes(x=mucus.viscosity, y=transformed_rate, group=country)) +
  #  geom_line(aes(color=arm))+
  geom_point(aes(color=country)) + 
  ggtitle("viscosity vs transformed rate of complete data ")


## this is  needed
ggplot(data=complete, aes(x=mucus.viscosity, y=transformed_rate, group=tissue.use)) +
  #  geom_line(aes(color=arm))+
  geom_point(aes(color=tissue.use)) + 
  ggtitle("viscosity vs transformed rate of complete data ")

## this is  needed
ggplot(data=active_only, aes(x=mucus.viscosity, y=transformed_rate, group=tissue.use)) +
  #  geom_line(aes(color=arm))+
  geom_point(aes(color=tissue.use)) + 
  ggtitle("viscosity vs transformed rate of active arm  data ONLY ")


## is significant only for mucus viscosity values higher than 2, in the active arm
## is NOT significant only for mucus viscosity values higher than 2, in the placebo arm


# this is needed
complete_high_tissue <- complete %>% filter(tissue.use == "HIGH")
ggplot(data=complete_high_tissue, aes(x=mucus.viscosity, y=transformed_rate, group=arm)) +
  #  geom_line(aes(color=arm))+
  geom_point(aes(color=arm)) + 
  ggtitle("tissue paper high use data only")



###########################
## SIMPLE INFERENTIAL TESTS
###########################

chisq.test(table(complete$country,complete$tissue.use)) # this is the one that is significant
chisq.test(table(complete$arm,complete$country))
chisq.test(table(complete$country,complete$eye.colour))
chisq.test(table(complete$arm,complete$eye.colour))
chisq.test(table(complete$tissue.use,complete$eye.colour))
chisq.test(table(complete$arm,complete$tissue.use))


t.test(complete$mucus.viscosity[complete$arm == "ACTIVE"],complete$mucus.viscosity[complete$arm == "PLACEBO"])
t.test(complete$mucus.viscosity[complete$tissue.use== "MEDIUM"],complete$mucus.viscosity[complete$tissue.use == "HIGH"]) # this is significant
t.test(active_only$mucus.viscosity[active_only$tissue.use== "MEDIUM"],active_only$mucus.viscosity[active_only$tissue.use == "HIGH"]) # this is significant
t.test(placebo_only$mucus.viscosity[placebo_only$tissue.use== "MEDIUM"],placebo_only$mucus.viscosity[placebo_only$tissue.use == "HIGH"]) # this is significant


t.test(complete$transformed_rate[complete$arm == "ACTIVE"],complete$transformed_rate[complete$arm == "PLACEBO"])
t.test(active_only$transformed_rate[complete$tissue.use== "MEDIUM"],active_only$transformed_rate[complete$tissue.use == "HIGH"]) # this is the question 1
t.test(placebo_only$transformed_rate[complete$tissue.use== "MEDIUM"],placebo_only$transformed_rate[complete$tissue.use == "HIGH"])  # this is the question 1

t.test(active_only$transformed_rate,placebo_only$transformed_rate) 
#this should have been important and was one of the important factors let the trial go to phase III
# not significant, however transformed_rate is not normal, so other type of non parametric test needs to be used.
wilcox.test(active_only$transformed_rate,placebo_only$transformed_rate, alternative = "greater")
wilcox.test(placebo_only$transformed_rate,active_only$transformed_rate, alternative = "greater")

wilcox.test(active_only$transformed_rate,placebo_only$transformed_rate, alternative = "greater")


t.test(active_only$mucus.viscosity,placebo_only$mucus.viscosity) # same as before 



summary(aov(active_only$transformed_rate ~ active_only$country))
summary(aov(placebo_only$transformed_rate ~ placebo_only$country))
summary(aov(complete$transformed_rate ~ complete$country))
summary(aov(complete$nosebleeds ~ complete$country)) # significant  this is question 2
summary(aov(complete$previous.year ~ complete$country)) # significant # this is question 2
summary(aov(complete$duration ~ complete$country)) # not significant
summary(aov(complete$mucus.viscosity ~ complete$country)) # significant



########################################
### REGRESSION MODELS
########################################


anova(model1,model2,model3)
lmtest::lrtest(model1,model2)

### the ya axis is not normal , so they have to be adjusted

plot(density(rgamma(2,1)))

model1 <- glm(complete$previous.year~complete$mucus.viscosity, family = poisson())
summary(model1)
#broom::tidy(model)
#anova(model)
#AIC(model)
#plot(model)
model2 <- glm(complete$nosebleeds~complete$mucus.viscosity, family = poisson())
summary(model2)
#anova(model)
#AIC(model)
#plot(model)
model3 <- glm(complete$nosebleeds~ (complete$mucus.viscosity + complete$arm +  complete$country + complete$eye.colour + complete$tissue.use ) * complete$duration +  complete$previous.year, family = poisson())
summary(model3)
#anova(model)
#AIC(model)
#plot(model)
model4 <- glm(complete$difference~complete$mucus.viscosity, family = poisson())
summary(model4)
#anova(model)
#AIC(model)
#plot(model)
model5 <- glm(complete$difference~complete$mucus.viscosity*complete$duration, family = poisson())
summary(model5)
#anova(model)
#AIC(model)
#plot(model)
model6 <- glm(complete$transformed_rate~complete$mucus.viscosity +1 )
summary(model6)
#anova(model)
#AIC(model)
#plot(model)
model7 <- glm(complete$transformed_rate~complete$mucus.viscosity + complete$arm +  complete$country + complete$eye.colour + complete$tissue.use)
summary(model7)
#anova(model)
#AIC(model)
#plot(model)
model8 <- glm(complete$transformed_rate~complete$mucus.viscosity + complete$arm +  complete$country + complete$eye.colour + complete$tissue.use + complete$tissue.use * complete$mucus.viscosity)
summary(model8)
#anova(model)
#AIC(model)
#plot(model)
model9 <- glm(complete$transformed_rate~complete$mucus.viscosity + complete$arm +  complete$country + complete$eye.colour + complete$tissue.use + 1)
summary(model9)
#anova(model)
#AIC(model)
#plot(model)
model10 <- glm(complete$transformed_rate~complete$mucus.viscosity + complete$arm +  complete$country + complete$eye.colour + complete$tissue.use + complete$tissue.use * complete$mucus.viscosity)
summary(model10)
#anova(model)
#AIC(model)
#plot(model)
model11 <- lmer(complete$transformed_rate~complete$mucus.viscosity + complete$arm +  (1 | complete$country)  + complete$tissue.use)
summary(model11)
#anova(model)
#AIC(model)
#plot(model)
model12 <- lm(complete$transformed_rate~complete$mucus.viscosity + complete$arm +   complete$country  +  complete$tissue.use)
summary(model12)
#anova(model)
#AIC(model)
#plot(model)
model13 <- glm(active_only$difference+3~active_only$mucus.viscosity, family = poisson())
summary(model13)
#anova(model)
#AIC(model)
#plot(model)
model14 <- glm(active_only$difference+3~active_only$mucus.viscosity*active_only$duration, family = poisson())
summary(model14)
#anova(model)
#AIC(model)
#plot(model)
model15 <- glm(active_only$difference +3  ~ active_only$mucus.viscosity + active_only$duration, family = poisson())
summary(model15)
#anova(model)
#AIC(model)
#plot(model)
model16 <- glm(active_only$transformed_rate ~ active_only$mucus.viscosity + 1 )
summary(model16)
#anova(model)
#AIC(model)
#plot(model)
model17 <- glm(active_only$transformed_rate ~ active_only$mucus.viscosity + active_only$tissue.use)
summary(model17)
#anova(model)
#AIC(model)
#plot(model)
model18 <- glm(active_only$transformed_rate ~ active_only$mucus.viscosity + active_only$tissue.use + active_only$mucus.viscosity * active_only$tissue.use)
summary(model18)
## this one is the best so far
#anova(model)
#AIC(model)
#plot(model)
model19 <- glm(active_only$transformed_rate ~ active_only$mucus.viscosity + active_only$tissue.use + active_only$country)
summary(model19)
#anova(model)
#AIC(model)
#plot(model)
model20 <- glm(active_only$transformed_rate ~ active_only$mucus.viscosity + active_only$tissue.use + active_only$country + active_only$mucus.viscosity * active_only$tissue.use)
summary(model20)
#anova(model)
#AIC(model)
#plot(model)
model21 <- glm(active_only$transformed_rate ~ active_only$mucus.viscosity + active_only$tissue.use + active_only$country + active_only$mucus.viscosity * active_only$tissue.use * active_only$country )
summary(model21)
#anova(model)
#AIC(model)
#plot(model)
model22 <- glm(active_only$transformed_rate~active_only$mucus.viscosity  +  active_only$country + active_only$eye.colour + active_only$tissue.use)
summary(model22)
#anova(model)
#AIC(model)
#plot(model)
model23 <- glm(active_only$transformed_rate~active_only$mucus.viscosity  +  active_only$country + active_only$eye.colour + active_only$tissue.use + active_only$tissue.use * active_only$mucus.viscosity)
summary(model23)
#anova(model)
#AIC(model)
#plot(model)
model24 <- glm(active_only$transformed_rate~active_only$mucus.viscosity  +  active_only$country + active_only$eye.colour + active_only$tissue.use + active_only$tissue.use * active_only$mucus.viscosity * active_only$country)
summary(model24)
#anova(model)
#AIC(model)
#plot(model)

plot(active_only$mucus.viscosity,active_only$transformed_rate)
plot(placebo_only$mucus.viscosity,placebo_only$transformed_rate)





