
library(tidyverse)
library(gridExtra)
library(ggbeeswarm)
library(mediation)

#liv2 <- data.frame(read_csv("C:/Users/andrew.east2/Desktop/buckmalelivers.csv"))
liv2 <- data.frame(read_csv("pfhxslivers.csv", col_types = cols(liver = col_double()))) %>%
  filter(sex=="female") %>%
  filter(!is.na(liver)) %>%
  droplevels() %>%
  ungroup()




## select colours, theme, and panel letter font for plot
cols <- RColorBrewer::brewer.pal(6, "Blues")
theme_set(theme_classic(base_size = 12))
letter <- theme(plot.title = element_text(face="bold"))




liv2 %>%
  mutate(
    ## convert Dose to a factor
    fac.dose=factor(dose),
    ## ratio adjustment
    ratio=liver/bodyweight,
    ## ANCOVA adjustment
    adj.liver=resid(lm(liver ~ bodyweight)) + mean(liver)) ->
  liv2


## calculate group means for plotting
liv2 %>%
  group_by(fac.dose) %>% 
  summarise(mean.liver=mean(liver), 
            mean.body=mean(bodyweight),
            mean.ratio=mean(ratio)) ->
  liv2.means



############# no axis control
png("sixplots.png", height=6, width=9, res=300, units="in")
grid.arrange(
  ggplot(liv2, aes(y=liver, x=fac.dose, group=fac.dose)) +
    geom_quasirandom(shape=21, color="royalblue", fill=cols[3], cex=1.2, width=0.1) +
    geom_segment(aes(x = as.numeric(fac.dose)-0.2, y = mean.liver,
                     xend = as.numeric(fac.dose) + 0.2, yend = mean.liver),
                 lwd=1.2, col="firebrick", data=liv2.means) + 
    ylab("Liver weight (g)") + xlab("Dose") + ggtitle("A") +
    letter, 
  
  ggplot(liv2, aes(y=bodyweight, x=fac.dose, group=fac.dose)) +
    geom_quasirandom(shape=21, color="royalblue", fill=cols[3], cex=1.2, width=0.1) +
    geom_segment(aes(x = as.numeric(fac.dose)-0.2, y = mean.body,
                     xend = as.numeric(fac.dose) + 0.2, yend = mean.body),
                 lwd=1.2, col="firebrick", data=liv2.means) + 
    ylab("Body weight (g)") + xlab("Dose") +ggtitle("B") +
    letter, 
  
  ggplot(liv2, aes(y=ratio, x=fac.dose, group=fac.dose)) +
    geom_quasirandom(shape=21, color="royalblue", fill=cols[3], cex=1.2, width=0.1) +
    geom_segment(aes(x = as.numeric(fac.dose)-0.2, y = mean.ratio,
                     xend = as.numeric(fac.dose) + 0.2, yend = mean.ratio),
                 lwd=1.2, col="firebrick", data=liv2.means) + 
    ylab("Liver/Body weight ratio") + xlab("Dose") + ggtitle("C") +
    letter, 
  
  ggplot(liv2, aes(y=liver, x=bodyweight, fill=fac.dose)) + 
    geom_point(shape=21, color="royalblue", size=2) +
    scale_fill_manual(name = "Dose", values=cols) +
    ylab("Liver weight (g)") +  xlab("Body weight (g)") +
    ggtitle("D") +
    letter + theme(legend.position = "none"), 
  
  ggplot(liv2, aes(y=ratio, x=bodyweight, fill=fac.dose)) + 
    geom_point(shape=21, color="royalblue", size=2) +
    scale_fill_manual(name = "Dose", values=cols) +
    ylab("Liver/Body weight ratio") +  xlab("Body weight (g)") + ggtitle("E") +
    letter + theme(legend.position = "none"), 
  
  ggplot(liv2, aes(y=adj.liver, x=bodyweight, fill=fac.dose)) + 
    geom_point(shape=21, color="royalblue", size=2) +
    scale_fill_manual(name = "Dose", values=cols) +
    ylab("ANCOVA adjusted liver weight (g)") +
    xlab("Body weight (g)") + ggtitle("F") +
    letter, 
  
  nrow=2, ncol=3)
dev.off()
####### end

## check how adjustments work
cor.test(~liver + bodyweight, data=liv2)
cor.test(~ratio + bodyweight, data=liv2)
cor.test(~adj.liver + bodyweight, data=liv2)


## analysis of liverweight
liv2.mod <- lm(liver ~ fac.dose, data=liv2)
summary(liv2.mod)
par(mfrow=c(2, 2))  
plot(liv2.mod)

## analysis of bodyweight
bw2.mod <- lm(bodyweight ~ fac.dose, data=liv2)
summary(bw2.mod)
par(mfrow=c(2, 2))  
plot(bw2.mod)

## analysis of ratios
rel2.mod <- lm(ratio ~ fac.dose, data=liv2)
summary(rel2.mod)
par(mfrow=c(2, 2))
plot(rel2.mod)

## calculate 95% CI
cis <- confint(rel2.mod)[-1, ]

par(mfrow=c(1,1))
plot(c(1:4) ~ coef(rel2.mod)[-1], xlim=c(0.0035,0.025))
segments(cis[, 1], 1:4, cis[, 2], 1:4)
abline(v=0, lty=2)













# meditation stuff



liv3 <- subset(liv2, fac.dose%in%c("0","14"))
liv3 <- droplevels(liv3)
liv3$dose10 <- ifelse(liv3$dose==0,0,1)
str(liv3)


## results
coef(summary(lm(liver ~ dose10, data=liv3)))
coef(summary(lm(liver ~ bodyweight + dose10, data=liv3)))
coef(summary(lm(bodyweight ~ dose10, data=liv3)))
coef(summary(lm(ratio ~ dose10, data=liv3)))



## recode for nicer plots
liv3 <- liv3 %>%
  mutate(cols=recode(dose10, `0`=cols[1], `1`=cols[4]),
         pch=recode(dose10, `0`=1, `1`=17),
         X=factor(recode(dose10, `0`="Control", `1`="High")))


#pdf("../figs/Sim.pdf", height=7, width=7)
#png("../figs/Sim.png", height=7, width=7, res=300, units="in")
grid.arrange(
  ggplot(liv3, aes(y=liver, x=X, group=X)) +
    geom_quasirandom(shape=21, color="royalblue", fill=liv3$cols, cex=1.95, width=0.1) +
    ylab("Liver weight (g)") + xlab("Dose") +  ggtitle("A") +
    letter,
  
  ggplot(liv3, aes(y=bodyweight, x=X, group=X)) +
    geom_quasirandom(shape=21, color="royalblue", fill=liv3$cols, cex=1.95, width=0.1) +
    ylab("Body weight (g)") + xlab("Dose") + 
    ggtitle("B") +
    letter,
  
  ggplot(liv3, aes(y=ratio, x=X, group=X)) +
    geom_quasirandom(shape=21, color="royalblue", fill=liv3$cols, cex=1.95, width=0.1) +
    ylab("Liver/Body weight") + xlab("Dose") + 
    ggtitle("C") +
    letter,
  
  ggplot(liv3, aes(y=liver, x=bodyweight, group=X)) +
    geom_point(shape=21, color="royalblue", fill=liv3$cols, cex=1.95) +
    ylab("Liver weight (g)") + xlab("Body weight (g)") +
    geom_smooth(method = "lm", se = FALSE, lwd=0.5) + 
    ggtitle("D") + 
    letter,
  
  nrow=2, ncol=2)
#dev.off()


## mediation model
m.mod <- lm(bodyweight ~ X, data=liv3)

## outcome model
y.mod <- lm(liver ~  bodyweight + X, data=liv3)

## fit model
med.mod <- mediate(m.mod, y.mod,
                   treat = "X", mediator = "bodyweight",
                   control.value = "Control", treat.value = "High",
                   robustSE = TRUE)

summary(med.mod)


par(las=1,
    mar=c(4.2, 6, 1, 1))
plot(med.mod, xlab="Change in liver weight (g)", yaxt="n")
axis(2, at=1:3, labels = c("Total\neffect", "Direct\neffect", "Mediated via\nBody weight"))
#dev.off()


##############################################################################################################################
# as of 06MAR2023, have not checked any of the treatment specific mediation stuff beyond high, belowNOT RUN ##################
##############################################################################################################################


liv3 <- subset(liv2, fac.dose%in%c("0","1000"))
liv3 <- droplevels(liv3)
liv3$dose10 <- ifelse(liv3$dose==0,0,1)
str(liv3)


## results
coef(summary(lm(liver ~ dose10, data=liv3)))
coef(summary(lm(liver ~ bodyweight + dose10, data=liv3)))
coef(summary(lm(bodyweight ~ dose10, data=liv3)))
coef(summary(lm(ratio ~ dose10, data=liv3)))



## recode for nicer plots
liv3 <- liv3 %>%
  mutate(cols=recode(dose10, `0`=cols[1], `1`=cols[4]),
         pch=recode(dose10, `0`=1, `1`=17),
         X=factor(recode(dose10, `0`="Control", `1`="MediumHigh")))


#pdf("../figs/Sim.pdf", height=7, width=7)
#png("../figs/Sim.png", height=7, width=7, res=300, units="in")
grid.arrange(
  ggplot(liv3, aes(y=liver, x=X, group=X)) +
    geom_quasirandom(shape=21, color="royalblue", fill=liv3$cols, cex=1.95, width=0.1) +
    ylab("Liver weight (g)") + xlab("Dose") + ylim(1,3) + ggtitle("A") +
    letter,
  
  ggplot(liv3, aes(y=bodyweight, x=X, group=X)) +
    geom_quasirandom(shape=21, color="royalblue", fill=liv3$cols, cex=1.95, width=0.1) +
    ylab("Body weight (g)") + xlab("Dose") + ylim(30,50) +
    ggtitle("B") +
    letter,
  
  ggplot(liv3, aes(y=ratio, x=X, group=X)) +
    geom_quasirandom(shape=21, color="royalblue", fill=liv3$cols, cex=1.95, width=0.1) +
    ylab("Liver/Body weight") + xlab("Dose") + #ylim(7, 15)
    ggtitle("C") +
    letter,
  
  ggplot(liv3, aes(y=liver, x=bodyweight, group=X)) +
    geom_point(shape=21, color="royalblue", fill=liv3$cols, cex=1.95) +
    ylab("Liver weight (g)") + xlab("Body weight (g)") + ylim(1,3) + xlim(30,50) +
    geom_smooth(method = "lm", se = FALSE, lwd=0.5) + 
    ggtitle("D") + 
    letter,
  
  nrow=2, ncol=2)
#dev.off()


## mediation model
m.mod <- lm(bodyweight ~ X, data=liv3)

## outcome model
y.mod <- lm(liver ~  bodyweight + X, data=liv3)

## fit model
med.mod <- mediate(m.mod, y.mod,
                   treat = "X", mediator = "bodyweight",
                   control.value = "Control", treat.value = "MediumHigh",
                   robustSE = TRUE)

summary(med.mod)


par(las=1,
    mar=c(4.2, 6, 1, 1))
plot(med.mod, xlab="Change in liver weight (g)", yaxt="n")
axis(2, at=1:3, labels = c("Total\neffect", "Direct\neffect", "Mediated via\nBody weight"))
#dev.off()





liv3 <- subset(liv2, fac.dose%in%c("0","500"))
liv3 <- droplevels(liv3)
liv3$dose10 <- ifelse(liv3$dose==0,0,1)
str(liv3)


## results
coef(summary(lm(liver ~ dose10, data=liv3)))
coef(summary(lm(liver ~ bodyweight + dose10, data=liv3)))
coef(summary(lm(bodyweight ~ dose10, data=liv3)))
coef(summary(lm(ratio ~ dose10, data=liv3)))



## recode for nicer plots
liv3 <- liv3 %>%
  mutate(cols=recode(dose10, `0`=cols[1], `1`=cols[4]),
         pch=recode(dose10, `0`=1, `1`=17),
         X=factor(recode(dose10, `0`="Control", `1`="Medium")))


#pdf("../figs/Sim.pdf", height=7, width=7)
#png("../figs/Sim.png", height=7, width=7, res=300, units="in")
grid.arrange(
  ggplot(liv3, aes(y=liver, x=X, group=X)) +
    geom_quasirandom(shape=21, color="royalblue", fill=liv3$cols, cex=1.95, width=0.1) +
    ylab("Liver weight (g)") + xlab("Dose") + ylim(1,3) + ggtitle("A") +
    letter,
  
  ggplot(liv3, aes(y=bodyweight, x=X, group=X)) +
    geom_quasirandom(shape=21, color="royalblue", fill=liv3$cols, cex=1.95, width=0.1) +
    ylab("Body weight (g)") + xlab("Dose") + ylim(28,50) +
    ggtitle("B") +
    letter,
  
  ggplot(liv3, aes(y=ratio, x=X, group=X)) +
    geom_quasirandom(shape=21, color="royalblue", fill=liv3$cols, cex=1.95, width=0.1) +
    ylab("Liver/Body weight") + xlab("Dose") + #ylim(7, 15)
    ggtitle("C") +
    letter,
  
  ggplot(liv3, aes(y=liver, x=bodyweight, group=X)) +
    geom_point(shape=21, color="royalblue", fill=liv3$cols, cex=1.95) +
    ylab("Liver weight (g)") + xlab("Body weight (g)") + ylim(1,3) + xlim(28,50) +
    geom_smooth(method = "lm", se = FALSE, lwd=0.5) + 
    ggtitle("D") + 
    letter,
  
  nrow=2, ncol=2)
#dev.off()


## mediation model
m.mod <- lm(bodyweight ~ X, data=liv3)

## outcome model
y.mod <- lm(liver ~  bodyweight + X, data=liv3)

## fit model
med.mod <- mediate(m.mod, y.mod,
                   treat = "X", mediator = "bodyweight",
                   control.value = "Control", treat.value = "Medium",
                   robustSE = TRUE)

summary(med.mod)


par(las=1,
    mar=c(4.2, 6, 1, 1))
plot(med.mod, xlab="Change in liver weight (g)", yaxt="n")
axis(2, at=1:3, labels = c("Total\neffect", "Direct\neffect", "Mediated via\nBody weight"))
#dev.off()





liv3 <- subset(liv2, fac.dose%in%c("0","250"))
liv3 <- droplevels(liv3)
liv3$dose10 <- ifelse(liv3$dose==0,0,1)
str(liv3)


## results
coef(summary(lm(liver ~ dose10, data=liv3)))
coef(summary(lm(liver ~ bodyweight + dose10, data=liv3)))
coef(summary(lm(bodyweight ~ dose10, data=liv3)))
coef(summary(lm(ratio ~ dose10, data=liv3)))



## recode for nicer plots
liv3 <- liv3 %>%
  mutate(cols=recode(dose10, `0`=cols[1], `1`=cols[4]),
         pch=recode(dose10, `0`=1, `1`=17),
         X=factor(recode(dose10, `0`="Control", `1`="Low")))


#pdf("../figs/Sim.pdf", height=7, width=7)
#png("../figs/Sim.png", height=7, width=7, res=300, units="in")
grid.arrange(
  ggplot(liv3, aes(y=liver, x=X, group=X)) +
    geom_quasirandom(shape=21, color="royalblue", fill=liv3$cols, cex=1.95, width=0.1) +
    ylab("Liver weight (g)") + xlab("Dose") + ylim(1,3) + ggtitle("A") +
    letter,
  
  ggplot(liv3, aes(y=bodyweight, x=X, group=X)) +
    geom_quasirandom(shape=21, color="royalblue", fill=liv3$cols, cex=1.95, width=0.1) +
    ylab("Body weight (g)") + xlab("Dose") + ylim(28,50) +
    ggtitle("B") +
    letter,
  
  ggplot(liv3, aes(y=ratio, x=X, group=X)) +
    geom_quasirandom(shape=21, color="royalblue", fill=liv3$cols, cex=1.95, width=0.1) +
    ylab("Liver/Body weight") + xlab("Dose") + #ylim(7, 15)
    ggtitle("C") +
    letter,
  
  ggplot(liv3, aes(y=liver, x=bodyweight, group=X)) +
    geom_point(shape=21, color="royalblue", fill=liv3$cols, cex=1.95) +
    ylab("Liver weight (g)") + xlab("Body weight (g)") + ylim(1,3) + xlim(28,50) +
    geom_smooth(method = "lm", se = FALSE, lwd=0.5) + 
    ggtitle("D") + 
    letter,
  
  nrow=2, ncol=2)
#dev.off()


## mediation model
m.mod <- lm(bodyweight ~ X, data=liv3)

## outcome model
y.mod <- lm(liver ~  bodyweight + X, data=liv3)

## fit model
med.mod <- mediate(m.mod, y.mod,
                   treat = "X", mediator = "bodyweight",
                   control.value = "Control", treat.value = "Low",
                   robustSE = TRUE)

summary(med.mod)


par(las=1,
    mar=c(4.2, 6, 1, 1))
plot(med.mod, xlab="Change in liver weight (g)", yaxt="n")
axis(2, at=1:3, labels = c("Total\neffect", "Direct\neffect", "Mediated via\nBody weight"))
#dev.off()



liv3 <- subset(liv2, fac.dose%in%c("0","125"))
liv3 <- droplevels(liv3)
liv3$dose10 <- ifelse(liv3$dose==0,0,1)
str(liv3)


## results
coef(summary(lm(liver ~ dose10, data=liv3)))
coef(summary(lm(liver ~ bodyweight + dose10, data=liv3)))
coef(summary(lm(bodyweight ~ dose10, data=liv3)))
coef(summary(lm(ratio ~ dose10, data=liv3)))



## recode for nicer plots
liv3 <- liv3 %>%
  mutate(cols=recode(dose10, `0`=cols[1], `1`=cols[4]),
         pch=recode(dose10, `0`=1, `1`=17),
         X=factor(recode(dose10, `0`="Control", `1`="LowLow")))


#pdf("../figs/Sim.pdf", height=7, width=7)
#png("../figs/Sim.png", height=7, width=7, res=300, units="in")
grid.arrange(
  ggplot(liv3, aes(y=liver, x=X, group=X)) +
    geom_quasirandom(shape=21, color="royalblue", fill=liv3$cols, cex=1.95, width=0.1) +
    ylab("Liver weight (g)") + xlab("Dose") + ylim(1,3) + ggtitle("A") +
    letter,
  
  ggplot(liv3, aes(y=bodyweight, x=X, group=X)) +
    geom_quasirandom(shape=21, color="royalblue", fill=liv3$cols, cex=1.95, width=0.1) +
    ylab("Body weight (g)") + xlab("Dose") + ylim(28,50) +
    ggtitle("B") +
    letter,
  
  ggplot(liv3, aes(y=ratio, x=X, group=X)) +
    geom_quasirandom(shape=21, color="royalblue", fill=liv3$cols, cex=1.95, width=0.1) +
    ylab("Liver/Body weight") + xlab("Dose") + #ylim(7, 15)
    ggtitle("C") +
    letter,
  
  ggplot(liv3, aes(y=liver, x=bodyweight, group=X)) +
    geom_point(shape=21, color="royalblue", fill=liv3$cols, cex=1.95) +
    ylab("Liver weight (g)") + xlab("Body weight (g)") + ylim(1,3) + xlim(28,50) +
    geom_smooth(method = "lm", se = FALSE, lwd=0.5) + 
    ggtitle("D") + 
    letter,
  
  nrow=2, ncol=2)
#dev.off()


## mediation model
m.mod <- lm(bodyweight ~ X, data=liv3)

## outcome model
y.mod <- lm(liver ~  bodyweight + X, data=liv3)

## fit model
med.mod <- mediate(m.mod, y.mod,
                   treat = "X", mediator = "bodyweight",
                   control.value = "Control", treat.value = "LowLow",
                   robustSE = TRUE)

summary(med.mod)


par(las=1,
    mar=c(4.2, 6, 1, 1))
plot(med.mod, xlab="Change in liver weight (g)", yaxt="n")
axis(2, at=1:3, labels = c("Total\neffect", "Direct\neffect", "Mediated via\nBody weight"))
#dev.off()

##############################################################################################################################
# as of 06MAR2023, have not checked any of the treatment specific mediation stuff beyond high, above/between NOT RUN ##################
##############################################################################################################################










# Simulation/Bayesian shenanigans
library(tidyverse)
library(mcmcplots)
library(rstan)
library(rstanarm)
library(gridExtra)
library(extrafont)

source("functions.R")

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


#liv2 <- data.frame(read_csv("C:/Users/andrew.east2/Desktop/buckmalelivers.csv"))

liv2 %>%
  mutate(
    ## convert Dose to a factor
    fac.dose = factor(dose),
    ## calc liver-free body weight
    bw = bodyweight - liver, 
    ## center body and liver weight (helps to set the prior)
    c.BodyWt = bw - mean(bw),
    c.LiverWt = liver - mean(liver)) -> 
  liv2


## get an idea of values to help set priors
tapply(liv2$c.BodyWt, liv2$fac.dose, mean)
mean(tapply(liv2$c.BodyWt, liv2$fac.dose, sd)) # sigma
sd(tapply(liv2$c.BodyWt, liv2$fac.dose, mean)) # grp_sigma

tapply(liv2$c.LiverWt, liv2$fac.dose, mean)
mean(tapply(liv2$c.LiverWt, liv2$fac.dose, sd)) # sigma
sd(tapply(liv2$c.LiverWt, liv2$fac.dose, mean)) # grp_sigma




## Hierarchical model with unequal variances
compiled.model <- stan_model("multilev_het_var_mouse.stan")

## data to feed into model
d.hier <- with(liv2,
               list(
                 N = nrow(liv2),
                 P = length(unique(fac.dose)),
                 grp = as.integer(as.numeric(fac.dose)),
                 bw = c.BodyWt,
                 lw = c.LiverWt,
                 prior_only = 0)
)

## draw samples from posterior
fit.hier <- sampling(compiled.model, data=d.hier, iter=10000, chains=3,
                     seed=123)

check_hmc_diagnostics(fit.hier)


## get draws from posterior
post.hier <- rstan::extract(fit.hier)

## calculate means and 95% HPDI for plotting
post.hier$beta %>%
  calc.stats() ->
  hier.bw.result

post.hier$gamma %>%
  calc.stats() ->
  hier.lw.result







## calculate effects for each dose
D1 <- calc.effects(post.hier$theta[, 2] - post.hier$theta[, 1],
                   post.hier$beta[, 2] - post.hier$beta[, 1],
                   post.hier$alpha)
D2 <- calc.effects(post.hier$theta[, 3] - post.hier$theta[, 1],
                   post.hier$beta[, 3] - post.hier$beta[, 1],
                   post.hier$alpha)
D3 <- calc.effects(post.hier$theta[, 4] - post.hier$theta[, 1],
                   post.hier$beta[, 4] - post.hier$beta[, 1],
                   post.hier$alpha)
D4 <- calc.effects(post.hier$theta[, 5] - post.hier$theta[, 1],
                   post.hier$beta[, 5] - post.hier$beta[, 1],
                   post.hier$alpha)



## density plots for the direct effects
dens.D1 <- density(D1$Direct, adjust=1)
dens.D2 <- density(D2$Direct, adjust=1)
dens.D3 <- density(D3$Direct, adjust=1)
dens.D4 <- density(D4$Direct, adjust=1)




#png("posteriors.png", height=8, width=6, res=300, units="in") 
par(mfrow=c(4,2),
    lend="square",
    mar=c(3,5,2,0.1),
    oma=c(0,0,2,0),
    font.axis=1)

caterplot(D1, style="plain", col="royalblue", pch=16, reorder=FALSE, 
          quantiles=list(outer=c(0.025,0.975),inner=c(0.25,0.75)),
          val.lim=c(-0.05, 0.5), cex.labels=1, lwd=c(1,3), bty="L")
abline(v=c(0,0.138), lty=c(2,2), col=c("black","darkgray"))
points(c(3:1) ~ apply(D1, 2, mean), pch=23, bg="white", col="royalblue")
mtext("1.6 mg/kg-d PFHxS", side=3, line=0, cex=1, adj=0)
mtext("Effect decomposition", side=3, line=2, cex=1, font=2)


plot(dens.D1, lwd=2, col="royalblue", ylab="", bty="L",
     xlab="", yaxt="n", main="", xlim=c(-0.05, 0.65), ylim=c(0, 17))
polygon(dens.D1$x, dens.D1$y, col=rgb(0, 0, 1, 0.25), border=NA)
abline(v=c(0,0.138), lty=c(2,2), col=c("black","darkgray"))
#mtext("Dose = 1.6", side=3, line=0, cex=1, adj=0)
legend("topright", bty = "n",
       legend=paste0("P(Eff<0.138) = ", round(mean(D1$Direct < 0.138), 2)))
mtext("Direct effect", side=3, line=2, cex=1, font=2)

caterplot(D2, style="plain", col="royalblue", pch=19, reorder=FALSE, 
          quantiles=list(outer=c(0.025,0.975),inner=c(0.25,0.75)),
          val.lim=c(-0.05, 0.5), cex.labels=1, lwd=c(1,3), bty="L")
abline(v=c(0,0.138), lty=c(2,2), col=c("black","darkgray"))
points(c(3:1) ~ apply(D2, 2, mean), pch=23, bg="white", col="royalblue")
mtext("3.5 mg/kg-d PFHxS", side=3, line=0, cex=1, adj=0)

plot(dens.D2, lwd=2, col="royalblue", ylab="", bty="L",
     xlab="", yaxt="n", main="", xlim=c(-0.05, 0.65), ylim=c(0,17))
polygon(dens.D2$x, dens.D2$y, col=rgb(0, 0, 1, 0.25), border=NA)
abline(v=c(0,0.138), lty=c(2,2), col=c("black","darkgray"))
#mtext("Dose = 3.5", side=3, line=0, cex=1, adj=0)
legend("topright", bty = "n",
       legend=paste0("P(Eff<0.138) = ", round(mean(D2$Direct < 0.138), 2)))


caterplot(D3, style="plain", col="royalblue", pch=19, reorder=FALSE, 
          quantiles=list(outer=c(0.025,0.975),inner=c(0.25,0.75)),
          val.lim=c(-0.05, 0.5), cex.labels=1, lwd=c(1,3), bty="L")
abline(v=c(0,0.138), lty=c(2,2), col=c("black","darkgray"))
points(c(3:1) ~ apply(D3, 2, mean), pch=23, bg="white", col="royalblue")
mtext("7 mg/kg-d PFHxS", side=3, line=0, cex=1, adj=0)

plot(dens.D3, lwd=2, col="royalblue", ylab="", bty="L",
     xlab="", yaxt="n", main="", xlim=c(-0.05, 0.65), ylim=c(0, 17))
polygon(dens.D3$x, dens.D3$y, col=rgb(0, 0, 1, 0.25), border=NA)
abline(v=c(0,0.138), lty=c(2,2), col=c("black","darkgray"))
#mtext("Dose = 7", side=3, line=0, cex=1, adj=0)
legend("topright", bty = "n",
       legend=paste0("P(Eff<0.138) = ", round(mean(D3$Direct < 0.138), 2)))


caterplot(D4, style="plain", col="royalblue", pch=19, reorder=FALSE, 
          quantiles=list(outer=c(0.025,0.975),inner=c(0.25,0.75)),
          val.lim=c(-0.05, 0.5), cex.labels=1, lwd=c(1,3), bty="L")
abline(v=c(0,0.138), lty=c(2,2), col=c("black","darkgray"))
points(c(3:1) ~ apply(D4, 2, mean), pch=23, bg="white", col="royalblue")
mtext("14 mg/kg-d PFHxS", side=3, line=0, cex=1, adj=0)

plot(dens.D4, lwd=2, col="royalblue", ylab="", bty="L",
     xlab="", yaxt="n", main="", xlim=c(-0.05, 0.65), ylim=c(0, 17))
polygon(dens.D4$x, dens.D4$y, col=rgb(0, 0, 1, 0.25), border=NA)
abline(v=c(0,0.138), lty=c(2,2), col=c("black","darkgray"))
#mtext("Dose = 14 mg/kg-d", side=3, line=0, cex=1, adj=0)
legend("topright", bty = "n",
       legend=paste0("P(Eff<0.138) = ", round(mean(D4$Direct < 0.138), 2)))




#dev.off()





## prob of safety
tapply(liv2$liver, liv2$fac.dose, mean)

## interested in a 20% effect
(0.69 * 0.2) # +/- 0.138 g

## proportion of posterior between +/- 20%
mean(D1$Direct > -0.138 & D1$Direct < 0.138)
mean(D2$Direct > -0.138 & D2$Direct < 0.138)
mean(D3$Direct > -0.138 & D3$Direct < 0.138)
mean(D4$Direct > -0.138 & D4$Direct < 0.138)



## proportion of posterior below and above threshold
mean(D1$Direct < -0.138)
mean(D1$Direct > 0.138)
1-sum(mean(D1$Direct < -0.138),mean(D1$Direct > 0.138))

mean(D2$Direct < -0.138)
mean(D2$Direct > 0.138)
1-sum(mean(D2$Direct < -0.138),mean(D2$Direct > 0.138))

mean(D3$Direct < -0.138)
mean(D3$Direct > 0.138)
1-sum(mean(D3$Direct < -0.138),mean(D3$Direct > 0.138))

mean(D4$Direct < -0.138)
mean(D4$Direct > 0.138)


#png("ROPE_at_1.6.png", height=4, width=4, res=300, units="in") 
par(mfrow=c(1,1),
    mar=c(4, 1, 2, 0.1))

plot(dens.D1, lwd=2, col="royalblue", ylab="", bty="L",
     xlab="", yaxt="n", main="", xlim=c(-0.2, 0.65), ylim=c(0, 20))
polygon(dens.D1$x, dens.D1$y, col=rgb(0, 0, 1, 0.25), border=NA)
abline(v=c(-0.138, 0.138), lty=2)
mtext("Dose = 1.6", side=3, line=0, cex=1, adj=0)
mtext("Change in liver weight (g)", side=1, line=2.5)

text(y=19, x=0, label="0.56", font=2)
text(y=19, x=0.138+0.05, label="0.44", font=2)
text(y=19, x= -0.138-0.05, label="0.00", font=2)

#dev.off()

#png("ROPE_at_3.5.png", height=4, width=4, res=300, units="in") 
par(mfrow=c(1,1),
    mar=c(4, 1, 2, 0.1))

plot(dens.D2, lwd=2, col="royalblue", ylab="", bty="L",
     xlab="", yaxt="n", main="", xlim=c(-0.2, 0.65), ylim=c(0, 20))
polygon(dens.D2$x, dens.D2$y, col=rgb(0, 0, 1, 0.25), border=NA)
abline(v=c(-0.138, 0.138), lty=2)
mtext("Dose = 3.5", side=3, line=0, cex=1, adj=0)
mtext("Change in liver weight (g)", side=1, line=2.5)

text(y=19, x=0, label="0.02", font=2)
text(y=19, x=0.138+0.05, label="0.98", font=2)
text(y=19, x= -0.138-0.05, label="0.00", font=2)

#dev.off()

#png("../figs/ROPE.png", height=4, width=4, res=300, units="in") 
par(mfrow=c(1,1),
    mar=c(4, 1, 2, 0.1))

plot(dens.D3, lwd=2, col="royalblue", ylab="", bty="L",
     xlab="", yaxt="n", main="", xlim=c(-0.2, 0.65), ylim=c(0, 20))
polygon(dens.D3$x, dens.D3$y, col=rgb(0, 0, 1, 0.25), border=NA)
abline(v=c(-0.138, 0.138), lty=2)
mtext("Dose = 7", side=3, line=0, cex=1, adj=0)
mtext("Change in liver weight (g)", side=1, line=2.5)

text(y=19, x=0, label="0.00", font=2)
text(y=19, x=0.138+0.05, label="1.00", font=2)
text(y=19, x= -0.138-0.05, label="0.00", font=2)

#dev.off()




## prob that effect is > 0.35 for each dose
mean(D1$Direct > 0.138)
mean(D2$Direct > 0.138)
mean(D3$Direct > 0.138)
mean(D4$Direct > 0.138)

probs <- c(mean(D1$Direct > 0.138),
           mean(D2$Direct > 0.138),
           mean(D3$Direct > 0.138),
           mean(D4$Direct > 0.138))
x <- c(1.6,3.5,7,14)

png("dose_response_area_inside_20percent_effect.png", height=4, width=4, res=300, units="in") 
par(mai=c(1,1,0.1,0.1))
plot(probs~x, log="x", type="b", ylim=c(0,1), ylab="Probability of > 0.138 g change in liver weight", xlab="mg/kg d PFNA, female Peromyscus mice\nlog scaled x-axis", bty="n")
abline(h=c(0,1),lty=2, col="darkgray")
dev.off()


# prob that at least one dose has an effect > 0.138
mean(D1$Direct > 0.138 |
       D2$Direct > 0.138 |
       D3$Direct > 0.138 |
       D4$Direct > 0.138)
