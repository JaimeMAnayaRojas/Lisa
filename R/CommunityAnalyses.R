#######################
##  Analysis of effects of schistocephalus introductions into copepods habitats
## 
#####################


library(boral)
library(data.table)
library(ggplot2)
library(corrplot)
#library(rethinking)  

rm(list=ls(all=TRUE))

set.seed(125)

LOS = function(x){
  
  p =  100 * length(which(x > 0))/length(x)
  
  return(p)
  
}


# Put the data together ---------------------------------------------------

# Fit treatments
setwd("C:/Users/Lisa/Desktop/FoM AG Kurtz/Analysis/Lisa/")
Data = read.csv("Data/Rotifers.csv")

#Data$Exposed = factor(Data$Exposed)
#levels(Data$Exposed)
#levels(Data$Exposed) = c("C", "E+")

Data$infection = factor(Data$infection)
levels(Data$infection)

Data$Exposure = factor(Data$infection)
levels(Data$Exposure) = c("C", "E", "E", "E", "E")
levels(Data$Exposure)
Data$Infection = factor(Data$Exposed:Data$Infected) #?

head(Data)



# Rotifer composition model ----------------------------------------------
names(Data)[3:7]

Y = data.frame(Data[ , 3:8])
head(Y)


##
example_mcmc_control <- list(n.burnin = 25000, n.iteration = 50000, 
                             n.thin = 30, seed = set.seed(1234))

testpath0 <- file.path("jagsboralmodel0.txt")


# Random effects
# RE = data.frame(family = Data$Block)

# pure latent model 


fit0 = boral(y= Y,  family = "negative.binomial", 
             lv.control = list(num.lv = 2), save.model = TRUE, 
             mcmc.control = example_mcmc_control, model.name = testpath0)

##


head(Data)
X1 = model.matrix(~ Exposure * time.point, Data)[,-1]
testpath2 <- file.path(tempdir(), "jagsboralmodel2.txt")
fit1 = boral(y= Y, X =X1,  family = "negative.binomial", 
             lv.control = list(num.lv = 2), save.model = TRUE, 
             mcmc.control = example_mcmc_control, model.name = testpath2, 
             ssvs.index = 1, calc.ics = T, ssvs.index = 1 )



treatcors2 <- get.enviro.cor(fit1, est = "median")
rescors2 <- get.residual.cor(fit1, est = "median")
rescors0 <- get.residual.cor(fit0, est = "median")


# Often used in other areas of multivariate statistics, 
# the trace may be interpreted as
# the amount of covariation explained by the latent variables. 
(rescors0$trace/rescors2$trace-1) *100# 


# Correlation plots -------------------------------------------------------

corrplot(rescors0$sig.cor, type = "lower", diag = FALSE,
         title = " Residual Correlations",  mar = c(3,0.5,2,1),
         tl.srt = 45)

corrplot(rescors2$sig.cor, type = "lower", diag = FALSE,
         title = " Residual Correlations",  mar = c(3,0.5,2,1),
         tl.srt = 45)


corrplot(treatcors2$sig.cor, type = "lower", diag = FALSE,
         title = "Correlations due to covariates",  mar = c(3,0.5,2,1),
         tl.srt = 45)


## Calcuate the explain variance from the treatments

varP = calc.varpart(fit1, groupX = 1:(fit1$num.X+1)  )
varP

varP = data.frame(varP$varpart.X)
varP$Factor = c("B0",names(data.frame(X1)))
varP = varP[-1,]

df_plot = reshape2::melt(varP, value.name = "Factor")

names(df_plot)[2:3] = c("Species", "RI")
df_plot$Factor = factor(df_plot$Factor)

levels(df_plot$Factor)
levels(df_plot$Factor) <- c("Exposed","Exposed:Time", "Time")
levels(df_plot$Factor)
df_plot$Factor = factor(df_plot$Factor, levels = c( "Time", "Exposed","Exposed:Time"))
df_plot$RI = df_plot$RI *100
df_plot$Species = factor(df_plot$Species)



lisa_plot = ggplot(df_plot, aes(x = Species, y = RI, fill = Factor)) +
  geom_bar(stat = "identity", colour="black") + theme_bw() + theme(panel.grid.major = element_blank())

lisa_plot

# saves the plots
png("Plots/RI_rotifers.png")
lisa_plot+scale_fill_brewer(palette="Set2") + ylab("Relative importance (%)") + ylim(c(0,100))
graphics.off()
#



##################################################
### considering individual parasite families #####
##################################################


head(Data)
X1 = model.matrix(~ infection * time.point, Data)[,-1]
testpath3 <- file.path(tempdir(), "jagsboralmodel3.txt")
fit3 = boral(y= Y, X =X1,  family = "negative.binomial", 
             lv.control = list(num.lv = 2), save.model = TRUE, 
             mcmc.control = example_mcmc_control, model.name = testpath3, 
             ssvs.index = 1, calc.ics = T, ssvs.index = 1 )



treatcors3 <- get.enviro.cor(fit3, est = "median")
rescors3 <- get.residual.cor(fit3, est = "median")
rescors0 <- get.residual.cor(fit0, est = "median")


# Often used in other areas of multivariate statistics, 
# the trace may be interpreted as
# the amount of covariation explained by the latent variables. 
(rescors0$trace/rescors3$trace-1) *100# 


# Correlation plots -------------------------------------------------------

corrplot(rescors0$sig.cor, type = "lower", diag = FALSE,
         title = " Residual Correlations",  mar = c(3,0.5,2,1),
         tl.srt = 45)

corrplot(rescors3$sig.cor, type = "lower", diag = FALSE,
         title = " Residual Correlations",  mar = c(3,0.5,2,1),
         tl.srt = 45)


corrplot(treatcors3$sig.cor, type = "lower", diag = FALSE,
         title = "Correlations due to covariates",  mar = c(3,0.5,2,1),
         tl.srt = 45)


## Calcuate the explain variance from the treatments

varP = calc.varpart(fit3, groupX = 1:(fit3$num.X+1)  )
varP

varP = data.frame(varP$varpart.X)
varP$Factor = c("B0",names(data.frame(X1)))
varP = varP[-1,]

df_plot = reshape2::melt(varP, value.name = "Factor")

names(df_plot)[2:3] = c("Species", "RI")
df_plot$Factor = factor(df_plot$Factor)

levels(df_plot$Factor)
levels(df_plot$Factor) <- c("IBB.2","IBB.2:Time","IBB.4","IBB.4:Time","LK.1","LK.1:Time","LK.2","LK.2:Time", "Time")
levels(df_plot$Factor)
df_plot$Factor = factor(df_plot$Factor, levels = c( "IBB.2","IBB.2:Time","IBB.4","IBB.4:Time","LK.1","LK.1:Time","LK.2","LK.2:Time", "Time"))
df_plot$RI = df_plot$RI *100
df_plot$Species = factor(df_plot$Species)



lisa_plot = ggplot(df_plot, aes(x = Species, y = RI, fill = Factor)) +
  geom_bar(stat = "identity", colour="black") + theme_bw() + theme(panel.grid.major = element_blank())

lisa_plot

# saves the plots
png("Plots/RI_rotiferfam2.png")
lisa_plot+ ylab("Relative importance (%)") + ylim(c(0,100))
graphics.off()
#




######################
svg("RI_Overallcells.svg", width = 6, height = 5)
ggplot(df_plot, aes(x = Factor, y = RI, fill = Factor)) +
  geom_boxplot() + theme_bw() + theme(panel.grid.major = element_blank()) + 
  scale_fill_brewer(palette="Set2") + ylab("Overall Relative Importance (%)")
graphics.off()


png("Species_raw.png", width = 7, height = 7)

names(Data)
df = cbind(Data[,c("Exposure", "time.point")],Y)

df <- melt(df, id.vars = c("Exposure","time.point"), variable.name = "Species")

df$variable = factor(df$variable, levels = c("C1", "C2", "C5", "C6", "C7", "C10", "C13", "C14"))



df$set = factor(df$Species)
df$Time = factor(df$time.point)

levels(df$set)
levels(df$set) <- c(1,1,1,2,2,2)

df$value = log10(df$value)

ggplot(df, aes(x = Time, y = value, fill = Exposure)) +
  geom_boxplot(alpha=0.5) + theme_bw() + theme(panel.grid.major = element_blank()) +  
  scale_fill_brewer(palette="Set2") + facet_wrap( set ~ Species) +
  geom_jitter(size=0.6, alpha=0.5, width = 0.1)

graphics.off()

library(mvabund)
Ya = mvabund(Y)

Data$Time = factor(Data$time.point)

X = as.data.frame(model.matrix(~ Exposure * time.point, Data)[,-1])
head(X)
tas.nb <- manyglm(Ya ~ X$ExposureE + X$time.point + X$`ExposureE:time.point`, family = "negative.binomial")

anova(tas.nb, p.uni= "adjusted", test = "LR" )

