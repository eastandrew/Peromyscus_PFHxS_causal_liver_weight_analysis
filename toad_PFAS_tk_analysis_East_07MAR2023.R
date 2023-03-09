######################## toad uptake kinetics ######################################


library(readr)
toad <- read_csv("P:/Tox/HEF/TOAD PFAS BIOAccumulation/stats/stackedToaddata_east_07MAR2023.csv")

toad$PFAS <- factor(toad$PFAS)
toad <- toad %>%
  mutate(
    period = case_when(
      day==0 ~ "background",
      day%in%c(1:28) ~ "uptake",
      day%in%c(28:56) ~ "elimination"
    )
  )
toad$period <- factor(toad$period, levels=c("background","elimination","uptake"))

library(tidyverse)

toad <- toad %>%
  pivot_wider(names_from=tissue, values_from=concentration_ug_per_kg) %>%
  mutate(liver_div_wholebody = liver/wholebody) %>%
  pivot_longer(cols=liver:wholebody, values_to="concentration_ug_per_kg", names_to="tissue")


ggplot(toad, aes(x=day, y=concentration_ug_per_kg)) +
  geom_point() +
  geom_smooth(method="loess", se=F, span=0.75) +
  facet_wrap(~fct_reorder(PFAS,concentration_ug_per_kg, .desc=T)*tissue, scales="free_y") +
  theme_classic()

ggplot(toad, aes(x=day, y=concentration_ug_per_kg)) +
  geom_point() +
  geom_smooth(method="loess", se=F, span=0.75) +
  facet_wrap(~fct_reorder(PFAS,concentration_ug_per_kg, .desc=T)*tissue, scales="free_y") +
  theme_classic() +
  scale_y_log10()

ggplot(toad, aes(x=day, y=liver_div_wholebody)) +
  geom_point() +
  geom_smooth(method="loess", se=F, span=0.75) +
  facet_wrap(~fct_reorder(PFAS,liver_div_wholebody, .desc=T), scales="free_y") +
  theme_classic() +
  scale_y_log10()
  

toadsummPFOSmed <- toad %>%
  filter(PFAS=="PFOS") %>%
  group_by(tissue,day) %>%
  summarize(medianPFOS = median(concentration_ug_per_kg, na.rm=T))

toad <- merge(toad, toadsummPFOSmed, by=c("tissue", "day"))
toad$relative_to_PFOSmed <- toad$concentration_ug_per_kg/toad$medianPFOS

toadsummday0med <- toad %>%
  filter(day==0) %>%
  group_by(tissue,PFAS) %>%
  summarize(medianday0 = median(concentration_ug_per_kg, na.rm=T))

toad <- merge(toad, toadsummday0med, by=c("PFAS", "tissue"))
toad$relative_to_day0 <- toad$concentration_ug_per_kg/toad$medianday0


ggplot(toad, aes(x=day, y=relative_to_PFOSmed)) +
  geom_point() +
  geom_smooth(method="loess", se=F, span=0.75) +
  facet_wrap(~fct_reorder(PFAS,relative_to_PFOSmed, .desc=T)*tissue, scales="free_y") +
  theme_classic()

ggplot(toad, aes(x=day, y=relative_to_day0)) +
  geom_point() +
  geom_smooth(method="loess", se=F, span=0.75) +
  facet_wrap(~fct_reorder(PFAS,relative_to_day0, .desc=T)*tissue, scales="free_y") +
  theme_classic()


pfosliver <- toad %>%
  filter(tissue=="liver"&PFAS=="PFOS") %>%
  droplevels()

pfosliversumm <- pfosliver %>%
  group_by(day) %>%
  summarize(mean=mean(concentration_ug_per_kg, na.rm=T),
            median=median(concentration_ug_per_kg, na.rm=T),
            sd=sd(concentration_ug_per_kg,na.rm=T))
pfosliversumm

















# need toad mass data
# need worm concentrations
# tween these two can solve for daily dose


# PFOS worm conc is 697 ug/kg
# toads fed 0.2g MWF, 0.0002kg*3=0.0006 weekly kg of worm, 0.0006kg*697ug/kg=0.4182ug weekly dose
D <- 0.4182/7

pfoslivup <- pfosliver %>% filter(day<=28) %>%droplevels()
pfoslivupmean <- pfoslivup %>% group_by(day)%>%summarize(mean=mean(concentration_ug_per_kg))

C = ((D*k01)/(V*(k01-k10)))*(e^-k10*t-e^-k01*t)


uplmall <- lm(log(concentration_ug_per_kg)~day, data=pfoslivup)
summary(uplmall)
plot(log(mean)~day, data=pfoslivupmean)
abline(uplmall)
# uptake (k1?) is 0.06

library(minpack.lm)

nlspfosmeanmodel <- nlsLM(mean~((0.4182*k01)/(V*(k01-k10)))*(exp(-k10*day)-exp(-k01*day)), data=pfoslivupmean, start=list(k01=0.06,k10=0.007, V=0.0007), trace=T, control=list(minFactor=0.000001))
summary(nlspfosmeanmodel)
preds <- predict(nlspfosmeanmodel, newdata=data.frame(day=seq(from=0, to=56, by=0.1)))
plot(mean~day, data=pfoslivupmean, ylim=c(0,2000), xlim=c(0,56))
points(preds~seq(from=0, to=56, by=0.1), type="l")



nlspfoswholemeanmodel <- nlsLM(mean~((0.4182*k01)/(V*(k01-k10)))*(exp(-k10*day)-exp(-k01*day)), data=pfosliversumm, start=list(k01=0.06,k10=0.007, V=0.0007), trace=T, control=list(minFactor=0.01))
summary(nlspfoswholemeanmodel)
predsw <- predict(nlspfoswholemeanmodel, newdata=data.frame(day=seq(from=0, to=56, by=0.1)))
plot(mean~day, data=pfosliversumm, ylim=c(0,2000), xlim=c(0,56))
points(predsw~seq(from=0, to=56, by=0.1), type="l")


# with intake  (not right, need to use two models to connect)
pfosliversumm$intake <- c(rep(0.4182,4),rep(0,5))
nlspfoswholemeanmodelin <- nlsLM(mean~((intake*k01)/(V*(k01-k10)))*(exp(-k10*day)-exp(-k01*day)), data=pfosliversumm, start=list(k01=0.06,k10=0.007, V=0.0007), trace=T, control=list(minFactor=0.01))
summary(nlspfoswholemeanmodelin)
newdata <- data.frame(day=seq(from=0, to=56, by=0.1)) %>% mutate(intake=case_when(
  day<28 ~ 0.4182,
  T ~ 0
))
predswin <- predict(nlspfoswholemeanmodelin, newdata=newdata)
plot(mean~day, data=pfosliversumm, ylim=c(0,2000), xlim=c(0,56))
points(predswin~seq(from=0, to=56, by=0.1), type="l")




nlspfoswholemodel <- nlsLM(concentration_ug_per_kg~259+((0.4182/7*k01)/(V*(k01-k10)))*(exp(-k10*day)-exp(-k01*day)), data=pfosliver, start=list(k01=0.113,k10=0.00007, V=0.0007), trace=T)
summary(nlspfoswholemodel)
predsall <- predict(nlspfoswholemodel, newdata=data.frame(day=seq(from=0, to=56, by=0.1)))
plot(concentration_ug_per_kg~day, data=pfosliver, ylim=c(0,2000), xlim=c(0,56), pch=16, cex=0.5)
points(predsall~seq(from=0, to=56, by=0.1), type="l")


nlspfoswholemodel2 <- nlsLM(concentration_ug_per_kg~259+((0.4182*0.06)/(V*(0.06-0.007)))*(exp(-0.007*day)-exp(-0.06*day)), data=pfosliver, start=list(V=0.0007), trace=T)
summary(nlspfoswholemodel2)
predsall2 <- predict(nlspfoswholemodel2, newdata=data.frame(day=seq(from=0, to=56, by=0.1)))
plot(concentration_ug_per_kg~day, data=pfosliver, ylim=c(0,2000), xlim=c(0,56), pch=16, cex=0.5)
points(predsall2~seq(from=0, to=56, by=0.1), type="l")

nlspfosmeanmodel2 <- nlsLM(median~259+((0.4182*0.06)/(V*(0.06-0.007)))*(exp(-0.007*day)-exp(-0.06*day)), data=pfosliversumm, start=list(V=0.0007), trace=T)
summary(nlspfosmeanmodel2)
predsall2b <- predict(nlspfosmeanmodel2, newdata=data.frame(day=seq(from=0, to=56, by=0.1)))
plot(median~day, data=pfosliversumm, ylim=c(0,2000), xlim=c(0,56), pch=16, cex=0.5)
points(predsall2b~seq(from=0, to=56, by=0.1), type="l")


library(deSolve)
model <- function(t, C, parms){
  with(as.list(c(C, parms)),{
    dC <- (ku*697)-(ke*697)*C
    
    list(c(dC))
  })
}

modele <- function(t, C, parms){
  with(as.list(c(C, parms)),{
    dC <- -ke*C
    
    list(c(dC))
  })
}





state <- c(C=250)
params <- c(ku=0.12, ke=0.00007)

times <- seq(0,28, by=1)

out <- ode(y=state, times=times, func=model, parms=params)
plot(C~time, data=out, log="", type="l")

statee <- c(C=max(out))
paramse <- c(ke=0.007)
timese <- seq(28,56, by=1)
oute <- ode(y=statee, times=timese, func=modele, parms=paramse)

plot(mean~day, data=pfosliversumm, ylim=c(0,2000), xlim=c(0,56), type="b", col="red", pch=16)
points(C~time, data=out, type="l")
points(C~time, data=oute, type="l")
points(concentration_ug_per_kg~day, data=pfosliver, pch=16, cex=0.5)
abline(v=28)




ssq=function(parms){
  # inital concentration
  cinit=c(C=250)
  # time points for which conc is reported
  # include the points where data is available
  t=c(seq(0,28,1),pfoslivupmean$day)
  t=sort(unique(t))
  # parameters from the parameter estimation routine
  ku=parms[1]
  ke=parms[1]
  # solve ODE for a given set of parameters
  out=ode(y=cinit,times=t,func=model,parms=c(ku=ku, ke=ke))
  # Filter data that contains time points where data is available
  outdf=data.frame(out)
  outdf=outdf[outdf$time %in% pfoslivupmean$day, ]
  # Evaluate predicted vs experimental residual
  preddf=outdf
  expdf=pfoslivupmean
  ssqres=preddf$C-expdf$mean
  # return predicted vs experimental residual
  #return(sum(ssqres^2))
  return(ssqres)
  
}


library(minpack.lm)
# parameter fitting using levenberg marquart algorithm
# initial guess for parameters
parms<-c(0.113,0.000007)
ssq(parms)

# fitting
fitval<-nls.lm(par=parms,fn=ssq)

fitval
summary(fitval)

opt <- optim(par=parms, fn=ssq, method="Nelder-Mead")
opt






pfosliverelim <- pfosliver %>% filter(day>=28)%>%droplevels()
pfosliverelimsumm <- pfosliverelim %>% group_by(day)%>%summarize(mean=mean(concentration_ug_per_kg, na.rm=T))
elimlm <- lm(log(mean)~day, data=pfosliverelimsumm)
summary(elimlm)
elimlmall <- lm(log(concentration_ug_per_kg)~day, data=pfosliverelim)
summary(elimlmall)
plot(log(mean)~day, data=pfosliverelimsumm)
abline(elimlmall)
# k2 is then -0.007
k2 <- -0.007
C0d <- mean(pfosliver$concentration_ug_per_kg[pfosliver$day==28])
I <- 0.002/0.020 #(0.002kg food/0.02kg toad)
Cfood <- 697
alpha <- ((C0d*k2)/(I*Cfood))*(1/(1-exp(-k2*0:56)))

plot(alpha~c(0:56), xlim=c(0,56), type="b", log="y")


pfoslivupmean$guess <- ((0.4182*0.309)/(0.0006*(0.309--0.025)))*(exp(--0.025*pfoslivupmean$day)-exp(-0.309*pfoslivupmean$day))
plot(mean~day, data=pfoslivupmean, ylim=c(0,2000))
points(guess~day, data=pfoslivupmean, type="b")


nlspfosmodel <- nls(concentration_ug_per_kg~((0.4182*k01)/(2.5*(k01-k10)))*(exp(-k10*day)-exp(-k01*day)), data=pfosliver, start=list(k01=0.02,k10=0.008), trace=T, control=list(minFactor=0.0000000001, warnOnly=T))

pfosliver$guesstimate <- ((0.4*0.9)/(0.25*(0.9-0.0008)))*(exp(-0.0008*pfosliver$day)-exp(-0.9*pfosliver$day))
plot(guesstimate~day, data=pfosliver)


