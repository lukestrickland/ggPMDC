load("tutorial_data/data_example_1.RData")

source("dmc/dmc.R")
load_model("LBA", "lbaN_B.R")

mapmeanv <- empty.map(list(S=c("n","w","p"),cond=c("C" ,"F", "H")
                           ,R=c("N","W","P")),
                      levels=c("CnN","CwN", "CnW","CwW",
                               "FnN","FwN","FpN","FnW","FwW","FpW","fa","FpP",
                               "HnN","HwN","HpN","HnW","HwW","HpW","HpP", 
                               "FAKERATE"))


mapmeanv[1:27] <- c(
  "CnN","CwN","FAKERATE",
  "FnN","FwN","FpN",
  "HnN","HwN","HpN",
  
  "CnW","CwW","FAKERATE",
  "FnW","FwW","FpW",
  "HnW","HwW","HpW",
  
  "FAKERATE","FAKERATE","FAKERATE",
  "fa","fa","FpP",
  "fa","fa","HpP"
  
)

#

#In this semi-factorial design, the is no PM accumulator under control 
#conditions. In order for model.dmc to know this, the user should do two things:
#1: input an 'N' argument corresponding to the number of accumulators for 
#each condition. Specify the number in constants. 
#2: set the thresholds to 'Inf' for non-existent accumulators. 

#In addition, set the dummy parameters "FAKERATE" to a constant.


model <- model.dmc(p.map = list(A="1",B=c("cond", "R"),t0="1",
                                mean_v=c("MAPMV"),sd_v="1",st0="1",N="cond"), 
                   match.map = list(M=list(n="N", w="W",p="P"),MAPMV=mapmeanv),
                   factors=list(S=c("n","w","p"),cond=c("C","F", "H")),
                   
                   constants = c(N.C=2,N.F=3, N.H =3, st0=0,B.C.P=Inf,
                                 mean_v.FAKERATE=1,sd_v=1), 
                   responses = c("N","W","P"), type="normN")

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

# A naive PM prior(similar to the study):
#uniform distribution for t0,
#truncated normals, truncated normal for A (0, 10),
#truncated normals for Bs (inf), normals for vs. 

p.prior <- prior.p.dmc(
  dists = rep("tnorm", 29),
  p1 = prior.p.mean,
  p2 = prior.p.scale,
  lower = c(rep(0, 9), .1, rep(NA, 19)),
  upper = c(rep(NA,9),  1, rep(NA, 19)))

#put the model to the test-can we recover real-ish PM parameter values 
# with the naive prior?
# 
# #generate data from "true" parameters that are different to the prior:
# sim.p.vector  <- c(t0=0.1,A=0.3,
#                    B.C.N = 1.3, B.C.W = 1.3, 
#                    B.H.N = 1.3,B.H.W = 1.5, B.H.P = 1.3, 
#                    B.F.N = 1.3,B.F.W = 1.4, B.F.P = 1.1,
#                    
#                    mean_v.CnN = 2.8, mean_v.CwN = -0.3,  
#                    mean_v.FnN = 2.8,mean_v.FwN = -0.3, mean_v.FpN = -1.6,
#                    mean_v.HnN = 2.8 , mean_v.HwN = -0.5, mean_v.HpN =-0.6, 
#                    
#                    mean_v.CnW= -1, mean_v.CwW = 2.9, 
#                    mean_v.FnW = -1, mean_v.FwW = 2.9, mean_v.FpW =0.5, 
#                    mean_v.HnW =-0.7,  mean_v.HwW = 3.0, mean_v.HpW =1.6,
#                    
#                    mean_v.fa= -2.4, mean_v.FpP =2.5 
#                    ,mean_v.HpP = 2.3)
# 
# 
# PM.data <- simulate.dmc(sim.p.vector, model, n=1e4)


###Get some samples
PM.data.model <- data.model.dmc(dat0,model)
PM.startpoints <- samples.dmc(nmc = 180,p.prior,PM.data.model, thin=20)
# Here we use RUN.dmc which is an automated method for getting converged samples.
# The function checks as number of criteria before accepting the final samples.
# However this automated method is not perfect, and the samples need to be
# inspected afterwards. 

system.time(
PM.samples.realparams <- RUN.dmc(PM.startpoints, cores= 8)
)

save(PM.samples.realparams, file="PM.samples.realparams.RData")
