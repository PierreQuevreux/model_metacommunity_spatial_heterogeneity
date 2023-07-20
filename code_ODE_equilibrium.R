# PACKAGES # ----
library(reshape2) # data formatting
library(tidyr) # data formatting
library(doParallel) # parallel computing
library(deSolve) #, ODE solver

source("code_functions.R") # functions describing the equations of the model
path_results="results/"

# PARAMETERS # ----
d_min=-5 # minimum log dispersal rate
d_max=5 # maximum log dispersal rate
d_interval=10^(seq(d_min,d_max,length.out=30))

g=1 # growth rate of the basal species
r=0 # respiration rate of predators
D=1 # self-regulation
e=0.65 # assimilation efficiency
m=c(0.0065,0.065,0.65,6.5,65) # predator to prey metabolic rate ratio
a=c(1/6.5,1/0.65,1/0.065) # attack rate relative to self regulation
sigma=1e-3 # perturbation variance
z=0.5 # scaling of perturbation variance with species equilibrium biomass

pert=list(c(1,1)) # c(species, patch)
disp=list(c(1,1)) # c(species1, species2,..., nSpecies)
asym=list(c(1.5,1.5),c(1,1)) # asymmetry of attack rate between patch #1 and #2

params_data_original<-expand.grid(simu_ID=0,
                                  g=g,
                                  r=r,
                                  D=D,
                                  e=e,
                                  m=m,
                                  a=a,
                                  sigma=sigma,
                                  z=z)
params_data_original$ma=params_data_original$m*params_data_original$a
params_data_original<-params_data_original[params_data_original$ma>0.05 & params_data_original$ma<=15,]
params_data_original$ea=params_data_original$e*params_data_original$a

########################## ----
# MAIN TEXT ############## ----
########################## ----
#### DISPERSAL OF PREDATORS #### ----
# parameters and simulations - gamma - omega # ----
nSpecies=2
nCommunity=2
gamma=3#seq(1,10,0.1) # asymmetry coefficient
#omega=seq(1,10,0.1) # asymmetry in growth rate
d=d_interval

params_data<-expand.grid(pert=list(list(c(1,1))), # perturbation of prey in patch 1
                         disp=list(c(0,1)), # dispersal of predators
                         gamma=gamma,
                         d=d,
                         model="pert_11")
params_data<-rbind(params_data,expand.grid(pert=list(list(c(1,2))), # perturbation of prey in patch 2
                                           disp=list(c(0,1)), # dispersal of predators
                                           gamma=gamma,
                                           d=d,
                                           model="pert_12"))
params_data$omega=params_data$gamma
params_data<-merge(params_data_original,params_data)
# sets asymmetry
for(i in 1:dim(params_data)[1]){
  params_data$asym[i]=set_asymmetry(nSpecies,params_data$gamma[i],params_data$omega[i])
}
# sets dispersal
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}

# select the values of ea and ma
params_data<-params_data[params_data$ea==1 & params_data$ma==1,]
params_data$simu_ID<-c(1:nrow(params_data))

# initial biomasses
B0<-equilibrium_analytical_TL2_disp_pred(params_data) # biomass at equilibrium

# integration parameters
t_max=1000
t_step=1

# simulations
start_time<-Sys.time()
print(start_time)
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)  
registerDoParallel(cl)
results<-foreach(i=1:nrow(params_data),.packages=c("deSolve")) %dopar% time_series(params_data[i,], B0[i,], t_max, t_step, nSpecies, nCommunity)
stopCluster(cl)
end_time<-Sys.time()
print(end_time)
end_time-start_time

TS<-create_TS(params_data, results, B0)
TS<-TS[,c("ea","ma","gamma","omega","model","time",colnames(B0))]

write.table(TS,paste0(path_results,"time_series.txt"),sep=",",row.names=F,col.names=T)
