rm(list=ls())
library("ape")
# library("devtools")
# install_github("hmorlon/PANDA")
library("RPANDA")

#read tree
tree.4g <- read.tree("./data/rosids_4g_OG_ex.tre")
tree.5g <- read.tree("./data/rosids_5g_OG_ex.tre")

#loading temperature data
data(InfTemp)
dof<-smooth.spline(InfTemp[,1], InfTemp[,2])$df


# define functions for time, lmbda, and mu
# constant rate
lambda.cst <- function(t,y){y[1]} #t=time, y=lambda (λ)
mu.cst <- function(t,y){y[1]} #t=time, y=mu (u)
# mu fixed to 0 (no extinction)
mu.0 <- function(t,y){0}

#variable rate
lambda.var <- function(t,y){y[1]*exp(y[2]*t)}
mu.var <- function(t,y){y[1]*exp(y[2]*t)}

#variable linear rate
lambda.l <- function(t,y){y[1] + y[2]*t}
mu.l <- function(t,y){y[1] + y[2]*t}

#rate varying as an exponetial function of temperature (x) and/or time (t)
#lambda is varying as an exponetial function of temperature (x)
lambda.x <-function(t,x,y){y[1] * exp(y[2] * x)}
#lambda is varying as exponential, function of time - no temperature dependency
lambda.t <- function(t,x,y){y[1] * exp(y[2] * t)} # should the same as lambda.var
#rate varying as an exponential, function of temperature (x) & time (t)
lamb.x.t <-function(t,x,y){y[1] * exp(y[2] * t + y[3] * x)}
mu.0.t <-function(t,x,y){0} # mu fixed to 0 (no extinction)

#define a fucntion evaluate all the models for given tree and par

fit.multi.rpanda <- function(tree,par) {
  # caculate crown age
  tot_time <- max(node.age(tree)$ages)
  # caculate fraction
  fraction=length(tree$tip.label)/114477
  # 1) B constant
  bcst.d0 <- fit_bd(tree, tot_time, f.lamb=lambda.cst, f.mu=mu.0, lamb_par=par[[1]][1], mu_par=c(), f=fraction, cst.lamb=TRUE, fix.mu=TRUE, cond="crown", dt=1e-3)
  # 2) BD constant
  bcst.dcst <- fit_bd(tree, tot_time, f.lamb=lambda.cst, f.mu=mu.cst, lamb_par=par[[2]][1], mu_par=par[[2]][2], cst.lamb=TRUE, cst.mu=TRUE, cond="crown", f=fraction, dt=1e-3)
  # 3) B variable E ("B exponential")
  bvar.d0 <- fit_bd(tree, tot_time, f.lamb=lambda.var, f.mu=mu.0, lamb_par=par[[3]][c(1,2)], mu_par=c(), expo.lamb=TRUE, fix.mu=TRUE, cond="crown", f=fraction, dt=1e-3)
  # 4) B variable L ("B linear")
  bvar.l.d0 <- fit_bd(tree, tot_time, f.lamb=lambda.l, f.mu=mu.0, lamb_par=par[[4]][c(1,2)], mu_par=c(), fix.mu=TRUE, f=fraction, cond="crown", dt=1e-3)
  # 5) B variable E, D constant ("B exponential, D constant")
  bvar.dcst <- fit_bd(tree, tot_time, f.lamb=lambda.var, f.mu=mu.cst, lamb_par=par[[5]][c(1,2)], mu_par=par[[5]][3], expo.lamb=TRUE, cst.mu=TRUE,cond="crown", f=fraction, dt=1e-3)
  # 6) B variable L, D constant ("B linear, D constant")
  bvar.l.dcst <- fit_bd(tree, tot_time, f.lamb=lambda.l, f.mu=mu.cst, lamb_par=par[[6]][c(1,2)],mu_par=par[[6]][3], cst.mu=TRUE, cond="crown", f=fraction, dt=1e-3)
  # 7) B constant, D variable E ("B constant, D exponential")
  bcst.dvar <- fit_bd(tree, tot_time, f.lamb=lambda.cst, f.mu=mu.var, lamb_par=par[[7]][1], mu_par=par[[7]][c(2,3)], cst.lamb=TRUE, expo.mu=TRUE, cond="crown", f=fraction, dt=1e-3)
  # 8) B constant, D variable L ("B constant, D linear")
  bcst.dvar.l <- fit_bd(tree, tot_time, f.lamb=lambda.cst, f.mu=mu.l, lamb_par=par[[8]][1], mu_par=par[[8]][c(2,3)], cst.lamb=TRUE, cond="crown", f=fraction, dt=1e-3)
  # 9) B variable, D variable
  bvar.dvar <- fit_bd(tree, tot_time, f.lamb=lambda.var, f.mu=mu.var, lamb_par=par[[9]][c(1,2)], mu_par=par[[9]][c(3,4)], expo.lamb=TRUE, expo.mu=TRUE, cond="crown", f=fraction, dt=1e-3)
  # 10) B exponential, function of temperature
  bvar.x.d0 <- fit_env(tree, InfTemp, tot_time, f.lamb=lambda.x, f.mu=mu.0.t, lamb_par=par[[10]][c(1,2)], mu_par=c(), fix.mu=TRUE, df=dof, cond="crown", f=fraction, dt=1e-3)
  # 11) B exponential, function of temperature & time
  bvar.x.t.d0 <- fit_env(tree, InfTemp, tot_time, f.lamb=lambda.x.t, f.mu=mu.0.t, lamb_par=par[[11]][c(1,2,3)], mu_par=c(), fix.mu=TRUE, df=dof, cond="crown", f=fraction, dt=1e-3)
  # 12) "B exponential, function of time - no temperature dependency"
  bvar.t.d0 <- fit_env(tree, InfTemp, tot_time, f.lamb=lambda.t, f.mu=mu.0.t, lamb_par=par[[12]][c(1,2)], mu_par=c(), expo.lamb=TRUE, fix.mu=TRUE, df=dof, cond="crown", f=fraction, dt=1e-3)
  
  #return results as a list
  return(list(bcstd0="bcst.d0", bcst.dcst="bcstdcst", bvar.d0="bvar.d0", bvar.l.d0="bvar.l.d0", bvar.dcst="bvar.dcst", bvar.l.dcst="bvar.l.dcst", bcst.dvar="bcst.dvar", 
              bcst.dvar.l="bcst.dvar.l", bvar.dvar="bvar.dvar", bvar.x.d0="bvar.x.d0", bvar.x.t.d0="bvar.x.t.d0", bvar.t.d0="bvar.t.d0"))
}

#giving pars
rosid.par <- list(c(0.09), c(0.09,0.3), c(0.05, 0.01), c(0.05, 0.01), c(0.05, 0.01, 0.1), c(0.05, 0.01, 0.5), c(0.05, 0.005, 0.01), c(0.05, 0.005, 0.001), c(0.4,-0.05,0.1,0.05), c(0.10, 0.01), c(0.10, -0.01, 0.03), c(0.05, 0.01)) 

rosid.4g.res <- fit.multi.rpanda(tree.4g, rosid.par)
saveRDS(rosid.4g.res, file="rosid.4g.res.rds")

tree <- tree.4g
par <- rosid.par
