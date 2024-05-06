# clear memory if necessary
rm(list = ls())
# load required packages
library(lme4)
library(lmerTest)
library(MuMIn)
library(tidyverse)
library(multcomp)
library(ggplot2)
library(effects)
library(ggthemes)
library(LMERConvenienceFunctions)
library(AER)
library(emmeans)
library(data.table)
options(scipen=999)

#### notes ####
# new coding scheme after LL rejection 
# for both production and priming sessions, code 6 and higher as a 1, everything else as a 0

# SSLA revision
# redid priming model to include all interactions
# redid production model as a single model, then used emmeans package to do post hocs
# made a new figure for the production data


####  priming session data ####

a <- read.csv("~/Dropbox/R/stranded preposition priming sessions.csv", header =TRUE) # load the data
a$subject<-as.factor(a$subject) # turn subjects into factors
str(a)

#### descriptives for logistic regression approach priming data ####

# turn DV into categorical
a$binary_score<-as.factor(a$binary_score)


# split into priming sessions
s1<-subset(a,Session=="session1")
s2<-subset(a,Session=="session2")

# slit data based on prime type for session 1
s1Prime<-subset(s1,Type=="prime")
s1NP<-subset(s1,Type=="non-prime")

# calculate frequency counts for session 1
tapply(s1Prime$binary_score,list(s1Prime$Modality),summary)
tapply(s1NP$binary_score,list(s1NP$Modality),summary)

# split data based on prime type for session 2
s2Prime<-subset(s2,Type=="prime")
s2NP<-subset(s2,Type=="non-prime")

# calculate frequency counts for session 2
tapply(s2Prime$binary_score,list(s2Prime$Modality),summary)
tapply(s2NP$binary_score,list(s2NP$Modality),summary)


#### check predictors in priming data for mulitcollinearity using VIF function #####

# make a df of the predictors only

preds<-a[c("Cloze","WMC","prod_pre","rec_pre_1_2")]

# run the VIF Function with a threshold of three
check.vif<-vif_func(in_frame=preds,thresh=3,trace=T)
# good to go there, max VIF is 1.41

#### check predictors in priming data for separation ####
ftable(xtabs(~binary_score+Modality,data=a))
ftable(xtabs(~binary_score+Session,data=a))
ftable(xtabs(~binary_score+Type,data=a))
ftable(xtabs(~binary_score+prod_pre,data=a))
ftable(xtabs(~binary_score+rec_pre_1_2,data=a))
ftable(xtabs(~WMC+binary_score,data=a))
ftable(xtabs(~Cloze+binary_score,data=a))
ftable(xtabs(~subject+binary_score,data=a))
ftable(xtabs(~verb+binary_score,data=a))



#### logistic regression (GLMM) of priming data #### 

# let's get ready to build the model
# center the numerical predictor variables (always good practice)
a$Cloze<-scale(a$Cloze,center=TRUE) # cloze score
a$typingSpeed<-scale(a$typingSpeed,center=TRUE) # typing speed (I don't think we need this)
a$WMC<-scale(a$WMC,center=TRUE) # working memory capacity
# a$rec_pre<-scale(a$rec_pre,center=TRUE) # score on the receptive pre test (GJT???) don't seem to use this one
a$prod_pre<-scale(a$prod_pre,center=TRUE) # score on the production pre test
a$rec_pre_1_2<-scale(a$rec_pre_1_2,center=TRUE) # some other receptive pre test score coded differently
a$trial_order<-scale(a$trial_order,center=TRUE) # captures the sequence of trials during session to control for task effects, etc

# build a model that contains only the main effects
# we increase the number of functions as well as the optimizer to avoid non-convergence
priming.model<-glmer(binary_score~Modality+trial_order+Session+Type+Cloze+WMC+prod_pre+rec_pre_1_2+(1+Type|subject)+(1|verb),a,family="binomial"(link="logit"),glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun = 100000)))
summary(priming.model)
test <- tidy(priming.model)
test

# check model assumptions and gather effect size
mcp.fnc(priming.model)
plot(priming.model)
qqPlot(resid(priming.model))
r.squaredGLMM(priming.model)

# calculate odds ratio for the main effects only model
# calculate confidence intervals using the "Wald" method
cc<-confint(priming.model,parm="beta_",level=0.90,method="Wald")
# save the confidence interval and coefficients to a table
ctab <- cbind(est=fixef(priming.model),cc)
# exponentiate the values to get odds rations
rtab <- exp(ctab)
# print them
print(rtab,digits=3)

# now a model with only significant interactions (but we won't be reporting this)
# priming.model.sig.interactions<-glmer(binary_score~Modality*Type+trial_order*Type+Session*Type+Cloze+WMC+prod_pre*Type+rec_pre_1_2+(1+Type|subject)+(1|verb),a,family="binomial"(link="logit"),glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun = 100000)))
# summary(priming.model.sig.interactions)
# r.squaredGLMM(priming.model.sig.interactions)

# I actually think the interaction b/w modality and WMC is not needed based on the research questions!
priming.model.full.interactions.m<-glmer(binary_score~Modality*Type+trial_order*Type+Session*Type+Cloze*Type+WMC*Type+prod_pre*Type+rec_pre_1_2*Type+(1+Type|subject)+(1|verb),a,family="binomial"(link="logit"),glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun = 100000)))
summary(priming.model.full.interactions.m)
r.squaredGLMM(priming.model.full.interactions.m)

# calculate odds ratio
# calculate confidence intervals using the "Wald" method
cc<-confint(priming.model.full.interactions.m,parm="beta_",level=0.80,method="Wald")
# save the confidence interval and coefficients to a table
ctab <- cbind(est=fixef(priming.model.full.interactions.m),cc)
# exponentiate the values to get odds rations
rtab <- exp(ctab)
# print them
print(rtab,digits=3)

# run post hocs analysis (not really necessary right?)
# within.groups.priming<-emmeans(priming.model.full.interactions.m,"Type",by="Modality",type="response")
# pairs(within.groups.priming)

# between.groups.priming<-emmeans(priming.model.full.interactions.m,"Modality",by="Type",type="response")
# pairs(between.groups.priming)
# emmip(priming.model.full.interactions.m,Modality~Type,type="response")


# I wonder if we ever need to provide the LME tables . . . 
# now let's plot significant interactions we want to display
# interaction between modality and prime type

# first use emmeans to plot the interaction
emmip(priming.model.full.interactions.m, Modality~Type, type="response")

# now make a ggplot object of the same interaction
# save the specific effect to a variable
ef<-effect("Modality*Type",priming.model.full.interactions.m)

# simple summary and plot
summary(ef)
plot(ef)

# but we want to make is prettier
# first convert the saved effect it into a data frame
dd<-as.data.frame(ef)

# now create ggplot object of the modality*type interaction
ggplot(dd,aes(Type,fit,linetype=Modality))+ #plot with type on x, fit(coef) on y, and linetype by modality
  geom_line(aes(group=Modality))+ #add a linetype with modality as the grouping aesthetic
  theme_classic(base_size=14)+ # use a clean theme
  geom_point(aes(group=Modality))+ # add points to the lines
  # theme_gdocs(base_size=14)+ #use a nice google theme with a larger font
  labs(x="",y="Primed Production (log odds)")+ #remove x axis label and change name of y label
  scale_x_discrete(expand=c(0.3,0),labels=c("Filler Trials","Prime Trials"))+ #reduce space between ticks on x axis and rename the labels
  theme(axis.title.y=element_text(face="bold",size=14))+ #change font and size of the y axis
  theme(axis.title.x=element_text(face="bold",size=14))+ #change font and size of the x axis
  scale_linetype_manual(values=c(1,2),name="Modality",guide=guide_legend(reverse=TRUE))+ #flip the ordering of the legend
  theme(legend.justification=c(-.5,.75),legend.position=c(.2,.9))+
  ylim(0,.5)

#### check "dispersion" with a function ####
# may not actually have to do this 
# (https://stats.stackexchange.com/questions/92156/how-to-handle-underdispersion-in-glmm-binomial-outcome-variable)

overdisp_fun <- function(model) {
  ## number of variance parameters in 
  ##   an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

# run the function on main effects and interaction models

overdisp_fun(priming.model)
overdisp_fun(priming.model.backfit)

# no idea what these mean, but apparently a similar approach
aods3::gof(priming.model)
aods3::gof(priming.model.backfit)


#### PRODUCTION DATA ####
# clear memory
rm(list = ls())

b <- read.csv("~/Dropbox/R/stranded preposition production sessions.csv", header =TRUE)
b$subject <- as.factor(b$subject)
str(b)


b$test_order <- factor(b$test_order,levels = c("a_pre", "post1", "post2"), labels =c("Pretest","Immediate Post", "Delayed Post"))
summary(b$test_order)
b$test_order <- relevel(b$test_order,ref="Pretest")

# check their default contrasts, which is treatment/dummy coding
contrasts(b$test_order)


dummy <- glmer(binary_score ~ test_order + (1|subject), data = b, family="binomial"(link="logit"))
summary(dummy)


contrasts(b$test_order) <- c(1, -.5, -.5)

sum_cont <- glmer(binary_score ~ test_order + (1|subject), data = b, family="binomial"(link="logit"))
summary(sum_cont)

# simple coding

original_contrasts <- contr.treatment(3)
my.coding <- matrix(rep(1/3, 6), ncol = 2)
original_contrasts
my.coding

my.simple <- original_contrasts - my.coding
my.simple

contrasts(b$test_order) <- my.simple

simple_contrasts <- glmer(binary_score ~ test_order + (1|subject), data = b, family="binomial"(link="logit"))
summary(simple_contrasts)


means <- b %>%
  group_by(test_order) %>%
  summarise(mean_score = mean(binary_score))

#  hhh

#### check mulitcollinearity in production data ####

# make a df of the predictors

prod_preds<-b[c("priming_amount","trial_order","cloze","WMC2","rec_pre_1_2","group","modality")]
car::vif(lm(binary_score~priming_amount+trial_order+cloze+WMC2+rec_pre_1_2+group+modality,data=b))

# run the VIF Function with a threshold of three
check.vif<-vif_func(in_frame=prod_preds,thresh=3,trace=T)
# good to go there, max VIF is 1.33


#### check separation in production data ####

ftable(xtabs(~binary_score+modality ,data=b))
ftable(xtabs(~binary_score+group ,data=b))
ftable(xtabs(~binary_score+condition+test_order,data=b))
ftable(xtabs(~binary_score+WMC2 ,data=b))
ftable(xtabs(~binary_score+test_order ,data=b))
ftable(xtabs(~binary_score+rec_pre ,data=b))
ftable(xtabs(~subject+binary_score ,data=b))
ftable(xtabs(~priming_amount+condition ,data=b))

ftable(xtabs(~priming_amount+test_order ,data=b))


#### descriptive stats for GLM ####
b$binary_score<-as.factor(b$binary_score)

# replace all my previous hard work with one simple command, lol
ftable(xtabs(~binary_score+condition+test_order,data=b))

#### production logistic model using a single model ####
# center our variables first

# have to take primimg amount out as it is MC with control group
b$cloze<-scale(b$cloze,center=TRUE)
b$WMC2<-scale(b$WMC2,center=TRUE)
b$rec_pre_1_2<-scale(b$rec_pre_1_2,center=TRUE)
b$trial_order<-scale(b$trial_order,center=TRUE)
# b$priming_amount<-scale(b$priming_amount,center=TRUE) # can't use in single model method

# create a main effects model
production.main.effects.model<-glmer(binary_score~group+modality+trial_order+test_order+rec_pre_1_2+cloze+WMC2+(1+test_order|subject)+(1|verb),b,family="binomial"(link="logit"),glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun = 100000)))

summary(production.main.effects.model)

# calculate odds ratio
# calculate confidence intervals using the "Wald" method
cc<-confint(production.main.effects.model,parm="beta_",level=0.90,method="Wald")
# save the confidence interval and coefficients to a table
ctab <- cbind(est=fixef(production.main.effects.model),cc)
# exponentiate the values to get odds rations
rtab <- exp(ctab)
# print them
print(rtab,digits=3)
r.squaredGLMM(production.main.effects.model)

# fit a model with interaction between group, modality, and test order
# also include interactions between individual differences and test order here. 
production.interaction.model<-glmer(binary_score~group*modality*test_order+trial_order+rec_pre_1_2+cloze+rec_pre_1_2+cloze+WMC2+(1+test_order|subject)+(1|verb),b,family="binomial"(link="logit"),glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun = 100000)))
summary(production.interaction.model)

# calculate odds ratio
# calculate confidence intervals using the "Wald" method
cc<-confint(production.interaction.model,parm="beta_",level=0.80,method="Wald")
# save the confidence interval and coefficients to a table
ctab <- cbind(est=fixef(production.interaction.model),cc)
# exponentiate the values to get odds rations
rtab <- exp(ctab)
# print them
print(rtab,digits=3)
r.squaredGLMM(production.interaction.model)


# run post hocs, first check each groups' performance over time
# gives estimate
compare.within.groups<-emmeans(production.interaction.model,c("test_order"),by=c("group","modality"))

# gives odds ratios
compare.within.groups<-emmeans(production.interaction.model,c("test_order"),by=c("group","modality"))
pairs(compare.within.groups,reverse=TRUE)
confint(pairs(compare.within.groups,type="response",reverse=TRUE),level=.9)


# now compare groups at each of the time points
# gives estimates
compare.between.groups<-emmeans(production.interaction.model,c("group","modality"),by=c("test_order"))

compare.between.groups<-emmeans(production.interaction.model,c("group","modality"),by=c("test_order"),type="response")

pairs(compare.between.groups,reverse=TRUE)
confint(pairs(compare.between.groups),reverse=TRUE,level=.90)


#### production data plots ####

# construct a plot using emmeans package
help(emmip)
emmip(production.interaction.model,group~modality | test_order,type="response")


# now create a custom ggplot object of the same interaction
ef<-effect("group*modality*test_order",production.interaction.model)
plot(ef)

ef<-as.data.frame(ef)
ef

summary(ef$test_order)
ef$test_order<-relevel(ef$test_order,ref="Pretest")
ef$test_order<-factor(ef$test_order,levels=c("Pretest","Immediate Post","Delayed Post"))
ef$modality<-revalue(ef$modality,c("CMC"= "SCMC"))
# ef$modality<-factor(ef$modality,levels=c("Pre-Test","Immediate Post","Delayed Post"),ordered=TRUE)

ggplot(ef,aes(modality,fit,shape=group,linetype=group))+
  facet_grid(.~test_order)+
  geom_point(aes(shape=group),size=3)+
  geom_line(aes(group=group))+
  theme_base(base_size=12)+
  scale_shape_manual(values=c(1,2),name="",labels=c("Control","Priming"))+
  scale_linetype_manual(values=c(1,2),name="",labels=c("Control","Priming"))+
  theme(axis.title.y=element_text(face="bold",size=12))+
  theme(legend.title=element_text(color="black",size=12,face="bold"))+
  theme(legend.text=element_text(face="bold", size=12))+
  labs(y="Stranded Preposition Production (log odds)",x="")+
  ylim(0,1)+
  theme(legend.justification=c(-.25,1),legend.position=c(0,.9))+
  theme(strip.text.x = element_text(size = 12,face="bold"))




# fit a model testing just the exp groups to examing effects of priming amount and other individual differences####
rm(list = ls())

# reload the data to fix centering and subsets
b<-read.csv("~/Dropbox/R/stranded preposition production sessions.csv", header =TRUE)
b$subject<-as.factor(b$subject)

# change it to a data table
exp.b<-setDT(b)

# select just the experimental group
exp.b<-b[group=="exp"]

tapply(exp.b$priming_amount,exp.b$modality,mean)
tapply(exp.b$priming_amount,exp.b$modality,sd)

# center predictors
exp.b$cloze<-scale(exp.b$cloze,center=TRUE)
exp.b$WMC2<-scale(exp.b$WMC2,center=TRUE)
exp.b$rec_pre_1_2<-scale(exp.b$rec_pre_1_2,center=TRUE)
exp.b$trial_order<-scale(exp.b$trial_order,center=TRUE)
exp.b$priming_amount<-scale(exp.b$priming_amount,center=TRUE)



# main effects only
exp.prod.only.main<-glmer(binary_score~cloze+rec_pre_1_2+priming_amount+WMC2+modality+(1+test_order|subject)+(1|verb),exp.b,family="binomial"(link="logit"),glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun = 100000)))
summary(exp.prod.only.main)

# calculate odds ratio
# calculate confidence intervals using the "Wald" method
cc<-confint(exp.prod.only.main,parm="beta_",level=0.80,method="Wald")
# save the confidence interval and coefficients to a table
ctab <- cbind(est=fixef(exp.prod.only.main),cc)
# exponentiate the values to get odds rations
rtab <- exp(ctab)
# print them
print(rtab,digits=3)
r.squaredGLMM(exp.prod.only.main)

# fit a model with interactions between test order and individual differences variables
exp.prod.only.test.order<-glmer(binary_score~cloze*test_order+rec_pre_1_2*test_order+priming_amount*test_order+WMC2*test_order+modality+(1+test_order|subject)+(1|verb),exp.b,family="binomial"(link="logit"),glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun = 100000)))
summary(exp.prod.only.test.order)

# calculate odds ratio
# calculate confidence intervals using the "Wald" method
cc<-confint(exp.prod.only.test.order,parm="beta_",level=0.80,method="Wald")
# save the confidence interval and coefficients to a table
ctab <- cbind(est=fixef(exp.prod.only.test.order),cc)
# exponentiate the values to get odds rations
rtab <- exp(ctab)
# print them
print(rtab,digits=3)
r.squaredGLMM(exp.prod.only.test.order)


# need a post hoc for comparison between immediate and delayed

my_comparison<-rbind(
  "immediate vs. delayed posttest"=c(0,0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0)
  
)

last.ph<-summary(glht(exp.prod.only.test.order, my_comparison), test = adjusted("none"))
confint(last.ph,level=.9)
exp(-0.338)


priming.effect<-effect("test_order*priming_amount",exp.prod.only.test.order)
plot(priming.effect)

emmip(exp.prod.only.test.order,priming_amount~test_order)
help(emmip)
