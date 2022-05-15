#######################
##  Analysis of effects of schistocephalus introductions into copepods habitats
## 
#####################

library(brms)
rm(list=ls(all=TRUE))

Data = read.csv("Data/Hatching.csv")

Data$batch = factor(Data$batch)

head(Data)


Bm1 <- brm(hachted | trials(total)  ~  Population * batch, family = binomial, data = Data,
           iter = 6000, warmup = 4000, chains =4, cores = 4,
           control = list(adapt_delta = 0.95, max_treedepth = 15))
summary(Bm1)

png("Plots/ModelPredictions_hatching.png")
conditional_effects(Bm1, effects = "Population:batch")
graphics.off()



Bm2 <- brm(hachted | trials(total)  ~  family * batch, family = binomial, data = Data,
           iter = 6000, warmup = 4000, chains =4, cores = 4,
           control = list(adapt_delta = 0.95, max_treedepth = 15))
summary(Bm2)

png("Plots/ModelPredictions_hatching_family.png")
conditional_effects(Bm2, effects = "family:batch")
graphics.off()
