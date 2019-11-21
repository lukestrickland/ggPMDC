
library(ggdmc)

load("~/ggPMDC/tutorial_data/singleS_2.RData")
load("~/ggPMDC/tutorial_data/PM.samples.realparams.RData")


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

posts <- summary(fit)

cbind(posts$statistics[,1], sim.p.vector[rownames(posts$statistics[,1:2])],
      posts$statistics[,2])


source("dmc/dmc.R")
load_model("LBA", "lbaN_B.R")

posts <- summary.dmc(PM.samples.realparams)

cbind(posts$statistics[,1], sim.p.vector[rownames(posts$statistics[,1:2])],
      posts$statistics[,2])

