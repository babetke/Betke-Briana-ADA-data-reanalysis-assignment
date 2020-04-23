# Data reanalysis proj. of Dan Becker's livestock abundance and vampire bat paper
# install the lme4 package
install.packages("lme4")
install.packages('lmerTest')
library(lme4)
library(lmerTest)# for da models
library(tidyverse) # for da tidy data
library(cowplot) # for da neato graphs
theme_set(theme_cowplot()) # so you dont have to remove grids from every graph...

# read in the dataset from the data file in the rproj, think about how it would work if 
# it were a link
d <- read_csv("data/becker et al_philtransbdata.csv",col_names = TRUE)
head(d)

# lets try recreating the bat isotope graph, figure S1. I do not have the isotopes for the 
# prey unfortunately so this is only the bat isotopes.

# also need to read in the wildlife data
d2 <- read_csv("data/prey_isotopes.csv")
head(d2, 3)

# Subsetting the sites
LR <- filter(d, site %in% c("LR1","LR2","LR3","LR4"))
head(LR)

AM <- filter(d, site %in% c("AM1","AM2","AM3","CA1"))
head(AM)

OW <- filter(d, site %in% c("OW1","OW2"))
OW

LR_prey <- filter(d2, site %in% c("LR1","LR2","LR3","LR4","Rio Nanay (near Iquitos)"))
AM_prey <- filter(d2, site %in% c("Puerto Pakuy (near AM1–3)","CA1"))
OW_prey <- filter(d2, site %in% c("OW2", "Pacbitun (near OW1 & OW2)"))

d2$class <- as.factor(d2$class)

# Creating the plots
pAM <- ggplot(NULL, aes(y = dn15, x = dc13)) + 
    geom_point(AM_prey, mapping = aes(shape = class)) +
    scale_shape_manual(values=c(18,15)) +
    geom_point(AM, mapping = aes(color = site)) +
    scale_color_manual(values = c("darkorchid4", "darkorchid3", "darkorchid1", "maroon3")) +
    theme(axis.line.x = element_blank(),
          axis.line.y = element_blank(),
          axis.title.y = element_blank(),
          legend.position = c(0.8, 0.8), 
          legend.text = element_text(size = 5),
          legend.title=element_text(size=5),
          panel.background = element_rect(color = 'black', size = 1)) +
    xlim(-30, -5) +
    ylim(4, 18)

pLR <- ggplot(NULL, aes(y = dn15, x = dc13)) + 
  geom_point(LR_prey, mapping = aes(shape = class)) +
  scale_shape_manual(values=c(18,15,24)) + 
  geom_point(LR, mapping = aes(color = site)) +
  scale_color_manual(values = c("forestgreen", "olivedrab", "olivedrab2", "olivedrab1")) +
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        legend.position = c(0.8, 0.8), 
        legend.text = element_text(size = 5),
        legend.title=element_text(size=5),
        panel.background = element_rect(color = 'black', size = 1)) +
  xlim(-30, -5) +
  ylim(4, 18)

pOW <- ggplot(NULL, aes(y = dn15, x = dc13)) + 
  geom_point(OW_prey, mapping = aes(shape = class)) +
  scale_shape_manual(values=c(18,15,24)) + 
  geom_point(OW, mapping = aes(color = site)) + 
  scale_color_manual(values = c("steelblue1", "steelblue4")) +
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = c(0.8, 0.8), legend.text = element_text(size = 5),
        legend.title=element_text(size=5),
        panel.background = element_rect(color = 'black', size = 1)) +
        xlim(-30, -5) +
        ylim(4, 18)

# then make a multiplot
plot_grid(pLR, pAM, pOW, labels = c('A', 'B', 'C'), nrow = 1, rel_widths = c(1.20,1,1))


# lmer
d$year <- as.factor(d$year)
d$livestock <- as.factor(d$livestock)
d$site <- as.factor(d$site)
d$ID <- as.factor(d$ID)

global_mod <- lmer(pca1)

mod_isodist <- lmer(idistance ~ livestock + year + (1|ID/site), data = d)
summary(mod_isodist)
Anova(mod_isodist)

iso_lm <- lm(idistance ~ livestock + year, data=d)
predict(iso_lm, newdata = nd, interval = 'confidence')

ranova(mod_isodist)
summary(aov(idistance ~ livestock, data = d))

nd <- data.frame(livestock = "low", year = "2013")
nd2 <- data.frame(livestock = "high", year = "2013")

predict(mod_isodist,newdata= nd,re.form=NA,
        interval = 'confidence')
predict(mod_isodist, newdata = nd, re.form = NA)


mod_isodist2 <- lmer(idistance ~ livestock + year + (1|ID), data = d)
PI <- predictInterval(merMod = mod_isodist2, newdata = nd, n.sims = 999)


d.no_na <- 


na.action = "na.fail"
mod <- lmer(pca1 ~ factor(rep) + year + (1|ID/site), data = d)
summary(mod)
Anova(mod)

install.packages("MuMIn")
library("MuMIn")

d.mod.build <- drop_na(d) %>%
  select(ID, site, livestock, idistance, year, age, sex, rep, pca1)
d.mod.build

mod <- lm(pca1 ~ livestock + idistance + year + age + sex + rep + livestock*sex + rep*sex + livestock*rep 
          + idistance*rep + idistance*sex, data = d.mod.build)
summary(mod)
combo3 <- dredge(mod, m.lim = c(1,4), beta = "partial.sd")
coefTable(combo3)
print(combo3)
model.avg(combo3, beta = "partial.sd")
summary(model.avg(combo3, subset = cumsum(weight) <= .95))
summary(model.avg(combo3, beta = "partial.sd"))

mod_avg <- summary(model.avg(combo3, subset = cumsum(weight) <= .95))






mod_lmer <- lmer(pca1 ~ livestock + idistance + year + age + sex + rep + livestock*sex + rep*sex + livestock*rep + idistance*rep 
                 + idistance*sex + (1|ID), data = d.mod.build)

mod_lmer2 <- lmer(pca1 ~ livestock + idistance + year + age + sex + rep + livestock*sex + rep*sex + livestock*rep + idistance*rep 
                 + idistance*sex + (1|site/ID), data = d) 


nit <- lm(pca1 ~ livestock + rep, data = d.mod.build)
summary(nit)

lmer(pca1 ~ livestock + rep + (1|ID/site), data = d.mod.build)





