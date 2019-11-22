rm(list=ls())
require(ggdmc)
#setwd("E:/Dropbox/code/ggPMDC")
## Model 1 --------------------------------------------
## 27 elements with 20 levels
FR <- list(S = c("n","w","p"), cond=c("C","F", "H"), R=c("N", "W", "P"))
lev <- c("CnN","CwN", "CnW","CwW",
         "FnN","FwN","FpN", "FnW","FwW","FpW", "fa","FpP",
         "HnN","HwN","HpN", "HnW","HwW","HpW", "HpP",
         "FAKERATE")
map_mean_v <- ggdmc:::MakeEmptyMap(FR, lev)
map_mean_v[1:27] <- c(
  "CnN","CwN","FAKERATE", "FnN","FwN","FpN", "HnN","HwN","HpN",
  "CnW","CwW","FAKERATE", "FnW","FwW","FpW", "HnW","HwW","HpW",
  "FAKERATE","FAKERATE","FAKERATE", "fa","fa","FpP", "fa","fa","HpP")

model0 <- BuildModel(
  p.map     = list(A = "1", B = c("cond", "R"), t0 = "1", mean_v = c("MAPMV"),
                   sd_v = "1", st0 = "1", N = "cond"),
  match.map = list(M = list(n = "N", w = "W", p = "P"), MAPMV = map_mean_v),
  factors   = list(S = c("n","w","p"), cond = c("C","F", "H")),
  constants = c(N.C = 2, N.F = 3, N.H = 3, st0 = 0, B.C.P = Inf,
                mean_v.FAKERATE = 1, sd_v = 1),
  responses = c("N", "W", "P"),
  type      = "norm")
#
npar <- length(GetPNames(model0))
#
# ## Population distribution, rate effect on F
sim.p.vector <- c(A = .3, B.C.N = 1.3,  B.F.N = 1.3,  B.H.N = 1.3,
              B.C.W = 1.3,  B.F.W = 1.4,  B.H.W = 1.5,
              B.F.P = 1.1,  B.H.P = 1.3,

              t0=.1,

              mean_v.CnN = 2.8,  mean_v.CwN = -0.3, mean_v.CnW=-1,
              mean_v.CwW = 2.9,  mean_v.FnN = 2.8,  mean_v.FwN=-.3,

              mean_v.FpN = -1.6, mean_v.FnW = -1,   mean_v.FwW = 2.9,
              mean_v.FpW = .5 ,  mean_v.fa = -2.4,  mean_v.FpP = 2.5,

              mean_v.HnN = 2.8, mean_v.HwN = -.5,   mean_v.HpN = -.6,
              mean_v.HnW = -.7, mean_v.HwW = 3.0,   mean_v.HpW = 1.6,
              mean_v.HpP = 2.3)

prior.p.mean <- c(A = 1, B.C.N = 1,  B.F.N = 1,  B.H.N = 1,
              B.C.W = 1,  B.F.W = 1,  B.H.W = 1,
              B.F.P = 1,  B.H.P = 1,
              
              t0=.3,
              
              mean_v.CnN = 1,  mean_v.CwN = 0, mean_v.CnW=0,
              mean_v.CwW = 1,  mean_v.FnN = 1,  mean_v.FwN=0,
              
              mean_v.FpN = 0, mean_v.FnW = 0,   mean_v.FwW = 1,
              mean_v.FpW = 0 ,  mean_v.fa = 0,  mean_v.FpP = 1,
              
              mean_v.HnN = 1, mean_v.HwN = 0,   mean_v.HpN = 0,
              mean_v.HnW = 0, mean_v.HwW = 1,   mean_v.HpW = 0,
              mean_v.HpP = 1)


prior.p.scale <-c(A = 1, B.C.N = 1,  B.F.N = 1,  B.H.N = 1,
              B.C.W = 1,  B.F.W = 1,  B.H.W = 1,
              B.F.P = 1,  B.H.P = 1,

              t0=.3,

              mean_v.CnN = 2,  mean_v.CwN = 2, mean_v.CnW = 2,
              mean_v.CwW = 2,  mean_v.FnN = 2,  mean_v.FwN = 2,

              mean_v.FpN = 2, mean_v.FnW = 2,   mean_v.FwW = 2,
              mean_v.FpW = 2,  mean_v.fa = 2,  mean_v.FpP = 2,

              mean_v.HnN = 2, mean_v.HwN = 2,   mean_v.HpN = 2,
              mean_v.HnW = 2, mean_v.HwW = 2,   mean_v.HpW = 2,
              mean_v.HpP = 2)

p.prior <- BuildPrior(
  dists = rep("tnorm", 29),
  p1 = prior.p.mean,
  p2 = prior.p.scale,
  lower = c(rep(0, 9), .1, rep(NA, 19)),
  upper = c(rep(NA,9),  1, rep(NA, 19)))

print(p.prior)
# 
# dat0 <- simulate(model0, nsim = 1e3, ps = sim.p.vector)
# 
# save(dat0, file="recovery_data/data_example_1.RData")

load("recovery_data/data_example_1.RData")

dmi0 <- BuildDMI(dat0, model0)


#note ncore argument is useless for single subject fits at the moment
#just uses 1 core no matter what you input. 
# But will be great when need to dispatch one subject
#per core
fit0 <- StartNewsamples(dmi0, prior=p.prior,nmc=180, ncore=1, thin=5)
save(fit0, model0, dat0, dmi0, sim.p.vector, p.prior, file = "recovery_data/singleS_1.RData")
fit  <- run(fit0, thin=5, ncore=1, block=FALSE)
save(fit, fit0, model0, dat0, dmi0, sim.p.vector, p.prior, file = "recovery_data/ggsim_gg_singleS.RData")

##Cross fit with dmc sim

load("recovery_data/dmc_data_example_1.RData")

dmi0 <- BuildDMI(dmc_dat0, model0)
fit0 <- StartNewsamples(dmi0, prior=p.prior,nmc=180, ncore=1, thin=5)
save(fit0, model0, dmc_dat0, dmi0, sim.p.vector, p.prior, file = "recovery_data/singleS_1.RData")
fit  <- run(fit0, thin=5, ncore=1, block=FALSE)
save(fit, fit0, model0, dmc_dat0, dmi0, sim.p.vector, p.prior, file = "recovery_data/singleS_2.RData")


save(dmc_fit, file="recovery_data/dmcsim_gg_singleS.RData")