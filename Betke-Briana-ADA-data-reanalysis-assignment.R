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

library("MuMIn")

d.mod.build <- drop_na(d) %>%
  select(ID, site, livestock, idistance, year, age, sex, rep, pca1)
d.mod.build

mod <- lm(pca1 ~ livestock + idistance + year + age + sex + rep + livestock*sex + rep*sex + livestock*rep 
          + idistance*rep + idistance*sex, data = d.mod.build)
summary(mod)
na.action = "na.fail"
combo3 <- dredge(mod, m.lim = c(1,4), beta = "partial.sd")
coefTable(combo3) # view the coefficient table
combo3
model.avg(combo3)
mod_avg <- summary(model.avg(combo3), subset = cumsum(weight) <= .95) # summary of the model
df1 <- as.data.frame(mod_avg$coefmat.full) # save the full model coefficients
CI <- as.data.frame(confint(mod_avg, full=T)) # get the 95% CI

total <- cbind(df1, CI) # bind them together

total <- total %>% # make the row names into a column called coefficent
  rownames_to_column("coefficient")

# make the plot!
ggplot(data=total[2:12,], aes(x = coefficient, y = Estimate)) + #again, excluding intercept because estimates so much larger
  geom_hline(yintercept = 0, color = "black", linetype="dashed")+ #add dashed line at zero
  geom_errorbar(aes(ymin=`2.5 %`, ymax=`97.5 %`), color="grey", #adj SE
                width=0, lwd=10) +
  #coord_flip()+ # flipping x and y axes
  geom_point(size=20, shape = "-")+ ylab("Coefficient")

# Working towards creating panel B
summary(lm(pca1 ~ livestock + rep, data=d.mod.build))
mod_best <- lm(pca1 ~ livestock + rep, data=d.mod.build)
nd <- data.frame(livestock = "low", rep = "Y")
nd2 <- data.frame(livestock = "low", rep = "N")
nd3 <- data.frame(livestock = "high", rep = "Y")
nd4 <- data.frame(livestock = "high", rep = "N")

predict(mod_best,newdata= nd, interval = 'confidence')
predict(mod_best,newdata= nd2, interval = 'confidence')
predict(mod_best,newdata= nd3, interval = 'confidence')
predict(mod_best,newdata= nd4, interval = 'confidence')






