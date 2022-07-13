# Script for survival analysis learning, following the tutorial in:
# https://www.emilyzabor.com/tutorials/survival_analysis_in_r_tutorial.html#what_is_survival_data
# As well as the user guides of the survival package and the following tutorials:
# https://rviews.rstudio.com/2017/09/25/survival-analysis-with-r/
# https://www.openintro.org/download.php?file=survival_analysis_in_R&referrer=/stat/surv.php
# http://www.sthda.com/english/wiki/cox-proportional-hazards-model
# https://stats.stackexchange.com/questions/119790/difference-in-chi-squared-calculated-by-anova-from-cph-and-coxph

# Install and/or load packages ----
NPacks <- c("survival","GGally","ggplot2", "survminer","knitr","emmeans","gtsummary")

#install.packages(c(NPacks,"pacman"))

## Load all the packages and try to install them if they are not (Uses/require pacman package)
pacman::p_load(char=NPacks)

# Load data ----
my.data<-read.csv(file.choose(),row.names=1,check.names=FALSE)
data(lung)
# Lung data set contains time: Survival time in days, status: censoring status 1=censored, 2=dead, sex: Male=1 Female=2
head(lung)
str(lung)

# The survival::Surv function accepts 1/0 (1=event), T/F (T/TRUE = event) or 1/2 (2=event)
# It creates a survival object for use as the response in a model formula. 
## Each entry correspond to one subject's survival time, followed by a + if the subject was censored
my.surv <- Surv(lung$time, lung$status)
my.surv[1:10] #This shows the first ten observations in the survival object

# Kaplan-Meier method for survival curves ----
## It is the most common method to estimate survival times and probabilities.
## It is non-parametric and results in a step function, with a step each time an event occurs.
my.surv <- Surv(lung$time, lung$status)
KM1 <- survfit(my.surv~1, data=lung) # This generates the overall survival curve 
                                   # Note: The 1 indicates that you only want a single survival curve
KM1 # The median survival time indicates the average survival time under non-normality)
   # This corresponds to the time where the survival probability is 0.5
names(KM1)
plot(KM1, xlab="Days",ylab="Overall survival probability")
GGally::ggsurv(KM1, back.white=T, xlab="Days",ylab="Overall survival probability")
survminer::ggsurvplot(KM1, xlab="Days",ylab="Overall survival probability")

## It is possible to estimate the x-time survival (probability of surviving to x-time)
summary(KM1,times = 365.25) # One-year survival estimation. In this case is 0.409, so 40.9%

## It is possible to compare survival times between groups
KM2<-survfit(my.surv~sex, data = lung) # One curve per sex group to compare between sexes
KM2
plot(KM2, xlab="Days",ylab="Survival probability for sexes")
surv.p <- ggsurv(KM2, back.white=T, xlab="Days",ylab="Survival probability for sexes")
surv.p
surv.p + ggplot2::guides(linetype=F) + 
  ggplot2::scale_color_discrete(name='sex',breaks=c(1,2),labels=c('Male','Female'))
surv.p <- ggsurvplot(KM2,conf.int = T, xlab="Days",ylab="Survival probability for sexes") 
surv.p$plot + scale_colour_viridis_d() + scale_fill_viridis_d()

# Between-group significance log-rank test
sd<-survdiff(my.surv~sex, data = lung)
sd

# Cox proportional hazards model ----
## It is a semi-parametric model that can be used to fit univariable and multivariable regression models with survival outcomes
## Assumes non-informative censoring and
## proportional hazards (The effect of the covariads on the survival is lineal)
CPH<-coxph(my.surv~sex,data=lung)
CPH
summary(CPH) # Shows a summary of the regression model
             # The exp(coef) is the HR
## The quantity of interest from a Cox regression model is a hazard ratio (HR).
## The HR represents the ratio of hazards between two groups at any particular point in time.
## The HR is interpreted as the instantaneous rate of occurrence of the event of interest in those who are still at risk for the event.
## It is not a risk, though it is commonly interpreted as such.
## A HR < 1 indicates reduced hazard of death whereas a HR > 1 indicates an increased hazard of death.
## The direction of the HR is interpreted with respect the first group in the data set.
## For the lung data example, the first subject is male (1), thus the 0.588 implies that females die less than males
## 0.588 times as many females are dying as males at any given time.

## For multivariate analysis it's necessary to test the proportional hazards assumption
my.surv<-Surv(lung$time, lung$status)
CPH2<-coxph(my.surv~sex+age#+strata(inst)
            ,data=lung)
CPH2
summary(CPH2)
anova(CPH2) #Not understood yet ?anova.coxph()
#emmeans(CPH2, list(pairwise ~sex), adjust = "tukey") #NOT WORKIGN YET
CZ<-cox.zph(CPH2) # Tests proportional hazards
CZ # A significant p-value indicates that the proportional hazards assumption is violated
plot(CZ) # Plots of the Schoenfeld residuals (produce two plots)
         # Deviation from a zero-slope line is evidence that the proportional hazards assumption is violated

CPH3<-coxph(my.surv~sex*age#+strata(inst)   # The * indicates we also evaluate the interaction of the variables
            ,data=lung)
CPH3
summary(CPH3)
CZ2<-cox.zph(CPH3) # Tests proportional hazards
CZ2 # A significant p-value indicates that the proportional hazards assumption is violated
plot(CZ2) # Plots of the Schoenfeld residuals (produce two plots)
# Deviation from a zero-slope line is evidence that the proportional hazards assumption is violated

anova(CPH2,CPH3) # Compares both models #Not understood yet

## Let's obtain and plot the survival curves
Cox.p<-survminer::ggadjustedcurves(CPH,variable="sex",method='conditional',xlab="Days",ylab="Survival")
Cox.p
Cox.p + ggplot2::guides(linetype=F) + 
  ggplot2::scale_color_discrete(name='sex',breaks=c(1,2),labels=c('Male','Female'))

Cox.p2<-survminer::ggadjustedcurves(CPH2,variable="sex"
                                    ,method='conditional',xlab="Days",ylab="Survival")
Cox.p2
Cox.p2 + ggplot2::guides(linetype=F) + 
  ggplot2::scale_color_discrete(name='sex',breaks=c(1,2),labels=c('Male','Female'))

## A different approach that allows ploting the confidence interval
### Consider that, we want to assess the impact of the sex on the estimated survival probability. 
### In this case, we construct a new data frame with two rows, one for each value of sex;
### the other covariates are fixed to their average values (if they are continuous variables) 
### or to their lowest level (if they are discrete variables). 
### For a dummy covariate, the average value is the proportion coded 1 in the data set.
### This data frame is passed to survfit() via the newdata argument:
sex_df <- with(lung, data.frame(sex = c(1, 2)
                                ,age = rep(mean(age, na.rm = TRUE), 2)
                                ,ph.ecog = c(1, 1)
                                )
               )
sex_df
fit <- survfit(CPH2, newdata = sex_df)
ggsurvplot(fit, data=sex_df, conf.int = TRUE, legend.labs=c("Sex=Male", "Sex=Female"),
                   ggtheme = theme_minimal())
