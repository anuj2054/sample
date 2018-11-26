
library(dplyr)
library(ggplot2)

##################
## DATA CLEANING
####################


complete <- dplyr::full_join(efficacy,randomization)
complete <- dplyr::full_join(complete,subject)

# removing the mucus viscorit that is NA
complete <- complete %>% filter(!(is.na(mucus.viscosity)))

## removing the outlier
complete <- complete %>% filter(!(rate == 146))
