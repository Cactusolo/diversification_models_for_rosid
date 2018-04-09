library("ape")
library("geiger")
library("laser")

rm(list=ls())

#setwd("C:/Users/cactus/Dropbox (UFL)/Exploring the phylogeny and diversification of rosids through a supermatrix approach with limited sampling/data/LTT/LTT_4g/Laser_analyses")
#setwd("/Users/cactus/Dropbox (UFL)/Exploring the phylogeny and diversification of rosids through a supermatrix approach with limited sampling/data/LTT/LASER/")
#read tree
tree.4g <- read.tree("./data/rosids_4g_OG_ex.tre")
tree.5g <- read.tree("./data/rosids_5g_OG_ex.tre")

is.ultrametric(tree.4g)
is.ultrametric(tree.5g)
dir.create("results")

#Estimation of Speciation and Extinction Rates With Birth-Death Models
bd.4g <- birthdeath(tree.4g)
bd.5g <- birthdeath(tree.5g)

(bd.4g)
(bd.5g)
##################laser and model test###########################

Btime.4g <- getBtimes(file = "./data/rosids_4g_OG_ex.tre")
Btime.5g <- getBtimes(file = "./data/rosids_5g_OG_ex.tre")

result.rosid.4g <- fitdAICrc(Btime.4g, modelset = c("pureBirth", "bd", "DDX", "DDL", "yule2rate", "yule3rate"), ints = 10000)
write.csv(result.rosid.4g, "./results/fitdAICrc.4g.csv")
saveRDS(result.rosid.4g, "./results/fitdAICrc.4g.rds")

result.rosid.5g <- fitdAICrc(Btime.5g, modelset = c("pureBirth", "bd", "DDX", "DDL", "yule2rate", "yule3rate"), ints = 10000)
saveRDS(result.rosid.5g, "./results/fitdAICrc.5g.rds")
write.csv(result.rosid.5g, "./results/fitdAICrc.5g.csv")

SPVAR.4g <- fitSPVAR(Btime.4g, init=c(0.075639167,0,0))
SPVAR.4g <- append(SPVAR.4g, list(z="NA"), 6)
SPVAR.5g <- fitSPVAR(Btime.5g, init=c(0.126169014,0,0))
SPVAR.5g <- append(SPVAR.5g, list(z="NA"), 6)


EXVAR.4g <- fitEXVAR(Btime.4g, init=c(0.075639167,0,0))
EXVAR.4g <- append(EXVAR.4g, list(k=0), 4)
EXVAR.5g <- fitEXVAR(Btime.5g, init=c(0.126169014,0,0))
EXVAR.5g <- append(EXVAR.5g, list(k="NA"), 4)

model.4g <- rbind.data.frame(SPVAR.4g, EXVAR.4g)
model.5g <- rbind.data.frame(SPVAR.5g, EXVAR.5g)

#put all the models and corresponding AIC value in the table then caculate AkaikeWeight

Tal.model <- c(as.character(result.rosid.4g$model), as.character(model.4g$model))

Tal.4g.aci <- c(as.numeric(result.rosid.4g$AIC), as.numeric(model.4g$aic))
Tal.5g.aci <- c(as.numeric(result.rosid.5g$AIC), as.numeric(model.5g$aic))

Tal.model.4g.5g <- cbind.data.frame(Tal.model, Tal.4g.aci, Tal.5g.aci)
#for model where both speciation and extinction rates can vary through time.
# BOTHVAR.4g <- fitBOTHVAR(Btime.4g, init=c(0.3,0.5,0.1,0.5)) #test with start value and a butch other values which is not recorded

#Error in if (muFx(mu0, z, Tmax) > lambdaFx(lam0, k, Tmax)) stop("Error - mu exceeds lambda\n") : 
#missing value where TRUE/FALSE needed

#current data not suit for the BOTHVAR model?


#plot rate
# SPVAR.4g$z <- 10000
# plotRate(Btime.4g, SPVAR.4g)
# plotRate(Btime.4g, EXVAR.4g)

#caclulate AkaikeWeight

library("BMhyd")

#AIC <- read.csv("./result/AIC_4g_5g_fianl.csv", header=T, stringsAsFactors = FALSE)
AIC <- Tal.model.4g.5g
Delta.AIC.4g <- AIC$Tal.4g.aci - min(AIC$Tal.4g.aci)
AkaikeWeight.4g <- AkaikeWeight(Delta.AIC.4g)

Delta.AIC.5g <- AIC$Tal.5g.aci - min(AIC$Tal.5g.aci)
AkaikeWeight.5g <- AkaikeWeight(Delta.AIC.5g)

AW_4g_5g <- rbind(AkaikeWeight.4g, AkaikeWeight.5g)

#colnames(AW_4g_5g) <- c("pureBirth", "bd", "DDX", "DDL", "yule2rate", "yule3rate", "SPVAR", "EXVAR")
colnames(AW_4g_5g) <- Tal.model
write.csv(AW_4g_5g, "./results/AkaikeWeight4g.5g.csv")
