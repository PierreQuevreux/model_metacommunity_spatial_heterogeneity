# PACKAGES # ----
library(ggplot2)
library(cowplot) # panel formatting
library(scales) # log scale
library(viridis) # colour blind friendly colour gradients
#library(magick)
library(reshape2) # data formatting
library(tidyr) # data formatting
#library(pracma)
library(doParallel) # parallel computing
library(deSolve) #, ODE solver
library(nleqslv) # root finding (function nleqslv)
library(numDeriv) # numerical Jacobian computing
library(rstudioapi) # to set the working directory

setwd(dirname(getActiveDocumentContext()$path)) # Set working directory to source file location
source("code_functions.R") # functions describing the equations of the model
source("code_figure_settings.R") # ggplot parameters and figure tuning
path_figure="Figures/"
path_data="Data/"

# PARAMETERS # ----
d_min=-5
d_max=5
d_step=0.1
d_interval=10^(seq(d_min,d_max,d_step))

g=1
r=0
D=1
e=0.65
m=c(0.0065,0.065,0.65,6.5,65)
a=c(1/6.5,1/0.65,1/0.065)
sigma=1e-3
z=0.5

pert=list(c(1,1)) # c(species, patch)
disp=list(c(1,1)) # c(species1, species2,..., nSpecies)
asym=list(c(1.5,1.5),c(1,1)) # asymmetry of attack rate in patch #1 and #2

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
params_data_original$ea<-as.factor(params_data_original$ea)
params_data_original$ma<-as.factor(params_data_original$ma)
levels(params_data_original$ma)<-c(ma01,ma1,ma10)
params_data_original$ma = factor(params_data_original$ma,levels(params_data_original$ma)[c(3,2,1)])
levels(params_data_original$ea)<-c(ea01,ea1,ea10)

########################## ----
# MAIN TEXT ############## ----
########################## ----
#### DISPERSAL OF PREDATORS #### ----
# parameters and simulations - gamma - omega # ----
nSpecies=2
nCommunity=2
gamma=seq(1,10,0.1) # asymmetry coefficient
#omega=seq(1,10,0.1) # asymmetry in growth rate
d=1e6 # high dispersal

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
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)

B<-equilibrium_analytical_TL2_disp_pred(params_data) # biomass at equilibrium
coexistence<-coexistence_TL2_disp_pred(params_data) # coexistence of all species

# index of the parameters leading to significant coexistence
index<-which(coexistence$coexistence_significative==1)
params_data<-params_data[index,]
B<-B[index,] # selects the biomasses allowing coexistence

# test
# i=3
# params_data_row<-params_data[i,]
# B<-B[i,]

# simulations
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1],.packages=NULL) %dopar% analysis(params_data[i,], B[i,], nSpecies, nCommunity)
stopCluster(cl)
list_data<-create_data(params_data, results, nSpecies, nCommunity)

## biomass # ----
data_B<-list_data$data_B
data_B<-data_B[data_B$model=="pert_11"
               & data_B$ea==paste(ea1) & data_B$ma==paste(ma1),
               which(names(data_B)%in%c("ea","ma","d","gamma","omega","model",list_data$B_names))]
data_B<-table_for_plot(data_B,6,"biomass")
data_B$biomass[data_B$biomass<1e-6]=NA

p1<-ggplot(data=data_B)+
      geom_line(aes(gamma,biomass,colour=species),size=1.5)+
      facet_wrap(~community,labeller=labeller(community=patch_labels))+
      corr_colour_TL_2+
      patch_line+
      theme +
      x_axis_gamma+
      xlab(label_gamma)+
      ylab("Biomass")

## biomass in isolated patches # ----
biomass<-cbind(params_data,as.data.frame(matrix(0,dim(params_data)[1],nSpecies*nCommunity)))
names(biomass)[(dim(params_data)[2]+1):dim(biomass)[2]]<-c(1:(nSpecies*nCommunity))

for(i in 1:dim(biomass)[1]){
  params<-as.list(c(g=params_data$g[i],
                    r=params_data$r[i],
                    D=params_data$D[i],
                    m=params_data$m[i],
                    a=params_data$a[i],
                    e=params_data$e[i]))
  params$asym=params_data$asym[[i]]
  params$disp=params_data$disp[[i]]
  biomass[i,c(1,2)+dim(params_data)[2]]<-equilibrium_isolated(params,nSpecies,1)
  biomass[i,c(3,4)+dim(params_data)[2]]<-equilibrium_isolated(params,nSpecies,2)
}
biomass<-biomass[biomass$ea==paste(ea1) & biomass$ma==paste(ma1)
                 & biomass$model=="pert_11",]

# biomasses with dispersal
data_B<-list_data$data_B
data_B<-data_B[data_B$model=="pert_11"
               & data_B$ea==paste(ea1) & data_B$ma==paste(ma1),
               which(names(data_B)%in%c("ea","ma","d","gamma","omega","model",list_data$B_names))]
data_B[,dim(data_B)[2]-(3:0)]<-data_B[,dim(data_B)[2]-(3:0)]/biomass[,dim(biomass)[2]-(3:0)]
data_B<-table_for_plot(data_B,6,"biomass")
data_B$biomass[data_B$biomass<1e-6]=NA

p2<-ggplot(data=data_B)+
      geom_line(aes(gamma,biomass,colour=species),size=1.5)+
      geom_hline(yintercept=1,linetype="dashed")+
      facet_wrap(~community,labeller=labeller(community=patch_labels))+
      corr_colour_TL_2+
      patch_line+
      theme +
      x_axis_gamma+
      xlab(label_gamma)+
      ylab("Scaled biomass")

# final figure # ----
graph<-ggdraw(xlim = c(0, 1.05), ylim = c(0, 2)) +
  draw_plot(p1, 0.05, 1, 1, 1)+
  draw_plot(p2, 0.05, 0, 1, 1)+
  draw_plot_label(c("A","B"), c(0,0), c(2,1), size = 30)
ggsave(paste(path_figure,"figure_biomass.pdf",sep=""), graph, width = 8, height = 8, device=cairo_pdf)

## biomass barplot # ----
biomass<-cbind(params_data,as.data.frame(matrix(0,dim(params_data)[1],nSpecies*nCommunity)))
names(biomass)[(dim(params_data)[2]+1):dim(biomass)[2]]<-c(1:(nSpecies*nCommunity))

for(i in 1:dim(biomass)[1]){
  params<-as.list(c(g=params_data$g[i],
                    r=params_data$r[i],
                    D=params_data$D[i],
                    m=params_data$m[i],
                    a=params_data$a[i],
                    e=params_data$e[i]))
  params$asym=params_data$asym[[i]]
  params$disp=params_data$disp[[i]]
  biomass[i,c(1,2)+dim(params_data)[2]]<-equilibrium_isolated(params,nSpecies,1)
  biomass[i,c(3,4)+dim(params_data)[2]]<-equilibrium_isolated(params,nSpecies,2)
}
names(biomass)[dim(biomass)[2]-(3:0)]<-list_data$B_names
biomass<-biomass[biomass$ea==paste(ea1) & biomass$ma==paste(ma1)
                 & biomass$model=="pert_11",
                 which(names(biomass)%in%c("ea","ma","d","gamma","omega","model",list_data$B_names))]
biomass<-table_for_plot(biomass,6,"biomass")
biomass$biomass<--biomass$biomass

data_B<-list_data$data_B
data_B<-data_B[data_B$model=="pert_11"
               & data_B$ea==paste(ea1) & data_B$ma==paste(ma1),
               which(names(data_B)%in%c("ea","ma","d","gamma","omega","model",list_data$B_names))]
data_B<-table_for_plot(data_B,6,"biomass")

data_B<-rbind(data_B,biomass)
data_B$species<-as.numeric(data_B$species)

p1<-ggplot(data=data_B[data_B$gamma==3,])+
      geom_rect(aes(xmin=species-0.5,xmax=species+0.5,ymin=0.0001,ymax=biomass,fill=as.factor(species)))+
      facet_wrap(~community)+
      coord_flip()+
      scale_fill_manual(values=c("dodgerblue3","chocolate1"))+
      theme_void()

ggsave(paste(path_figure,"figure_spillover.pdf",sep=""), p1, width = 6, height = 2, device = cairo_pdf)

## correlation inter-patch # ----
data_C<-list_data$data_C
databis<-data_C[data_C$ea==paste(ea1) & data_C$ma==paste(ma1),
                which(names(data_C)%in%c("ea","ma","d","gamma","model","C_11_12",'C_21_22'))]
h_value=databis$C_11_12[databis$gamma==1][1]
databis<-table_for_plot(databis,5,"correlation")

p1<-ggplot(data=databis)+
      geom_line(aes(gamma,correlation,colour=species,linetype=model),size=1.5)+
      geom_hline(yintercept=h_value,linetype="dashed")+
      #facet_grid(ma~ea, labeller=label_parsed)+
      perturbation_prey_line+
      corr_colour_TL_2+
      theme+theme(legend.position="bottom",legend.direction="vertical")+
      x_axis_gamma+
      ylim(-1,1)+
      xlab(label_gamma)+
      ylab(label_correlation)
ggsave(paste(path_figure,"figure_correlation.pdf",sep=""),p1, width = 6, height = 7, device = cairo_pdf)

## biomass CV at population scale # ----
data_CV_pop<-list_data$data_CV_pop
data_CV_pop<-data_CV_pop[data_CV_pop$ea==paste(ea1) & data_CV_pop$ma==paste(ma1),
                         which(names(data_CV_pop)%in%c("ea","ma","d","gamma","omega","model",list_data$CV_pop_names))]
data_CV_pop<-table_for_plot(data_CV_pop,6,"CV")
levels(data_CV_pop$model)<-c(pert_11,pert_12)

p2<-ggplot(data=data_CV_pop)+
      geom_line(aes(gamma,CV,colour=species,linetype=community),size=1.5)+
      facet_wrap(~model)+
      corr_colour_TL_2+
      patch_line+
      theme +
      x_axis_gamma+
      y_axis_log10+
      xlab(label_gamma)+
      ylab(label_CV_pop)

## biomass CV at metapopulation scale # ----
data_CV_species<-list_data$data_CV_species
data_CV_species<-data_CV_species[data_CV_species$ea==paste(ea1) & data_CV_species$ma==paste(ma1),
                                 which(names(data_CV_species)%in%c("ea","ma","d","gamma","omega","model",list_data$CV_species_names))]
data_CV_species<-table_for_plot(data_CV_species,6,"CV")
data_CV_species$species=data_CV_species$community
data_CV_species$community=NULL
levels(data_CV_species$model)<-c(pert_11,pert_12)

p3<-ggplot(data=data_CV_species)+
      geom_line(aes(gamma,CV,colour=species),size=1.5)+
      facet_wrap(~model)+
      corr_colour_TL_2+
      theme +
      x_axis_gamma+
      y_axis_log10+
      xlab(label_gamma)+
      ylab(label_CV_metapop)

## biomass CV at metacommunity scale # ----
data_CV<-list_data$data_CV
levels(data_CV$model)<-c(pert_11,pert_12)

p4<-ggplot(data=data_CV[data_CV$ea==paste(ea1) & data_CV$ma==paste(ma1),])+
      geom_line(aes(gamma,CV_tot),size=1.5)+
      facet_wrap(~model)+
      corr_colour_TL_2+
      theme +
      x_axis_gamma+
      y_axis_log10+
      xlab(label_gamma)+
      ylab(label_CV_metacom)

# final figure # ----
graph<-ggdraw(xlim = c(0, 2.08), ylim = c(0, 2)) +
  draw_plot(p2, 1.08, 1, 1, 1)+
  draw_plot(p3, 0.08, 0, 1, 1)+
  draw_plot(p4, 1.08, 0, 0.85, 1)+
  draw_plot_label(c("A","B","C","D"), c(0,1,0,1), c(2,2,1,1), size = 30)
ggsave(paste(path_figure,"figure_CV.pdf",sep=""), graph, width = 16, height = 10, device = cairo_pdf)

## terms of the jacobian matrix # ----
data_J<-list_data$data_J
data_J<-data_J[data_J$model=="pert_11" & data_J$ea==paste(ea1) & data_J$ma==paste(ma1),
                         which(names(data_J)%in%c("ea","ma","d","gamma","omega","model",c("J_21_11","J_11_21","J_22_12","J_12_22")))]
data_J[data_J$gamma==3,]
data_J<-melt(data_J,
             id.vars = names(data_J)[1:6],
             variable.name = "interaction",
             value.name = "jacobian")

ylim=c(min(data_J$jacobian),max(data_J$jacobian))
colour=viridis(4)

p1<-ggplot(data=data_J[data_J$interaction=="J_21_11",])+
  geom_line(aes(gamma,jacobian),colour=colour[1],size=1.5)+
  theme +
  x_axis_gamma+
  ylim(ylim)+
  xlab(expression(gamma))+
  ylab(label_jacobian)+
  ggtitle("J_21_11")

p2<-ggplot(data=data_J[data_J$interaction=="J_11_21",])+
  geom_line(aes(gamma,jacobian),colour=colour[3],size=1.5)+
  theme +
  x_axis_gamma+
  ylim(ylim)+
  xlab(expression(gamma))+
  ylab(label_jacobian)+
  ggtitle("J_11_21")

p3<-ggplot(data=data_J[data_J$interaction=="J_22_12",])+
  geom_line(aes(gamma,jacobian),colour=colour[2],size=1.5)+
  theme +
  x_axis_gamma+
  ylim(ylim)+
  xlab(expression(gamma))+
  ylab(label_jacobian)+
  ggtitle("J_22_12")

p4<-ggplot(data=data_J[data_J$interaction=="J_12_22",])+
  geom_line(aes(gamma,jacobian),colour=colour[4],size=1.5)+
  theme +
  x_axis_gamma+
  ylim(ylim)+
  xlab(expression(gamma))+
  ylab(label_jacobian)+
  ggtitle("J_12_22")

graph<-ggdraw(xlim = c(0, 2), ylim = c(0, 2)) +
  draw_plot(p1, 0, 1, 1, 1)+
  draw_plot(p2, 1, 1, 1, 1)+
  draw_plot(p3, 0, 0, 1, 1)+
  draw_plot(p4, 1, 0, 1, 1)
ggsave(paste(path_figure,"figure_jacobian.pdf",sep=""), graph, width = 6, height = 5, device = cairo_pdf)

## time series # ----
nSpecies=2
nCommunity=2
gamma=3 # asymmetry coefficient
d=1e6

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
params_data<-merge(params_data_original,params_data)
params_data$omega=params_data$gamma
params_data<-params_data[params_data$ea==paste(ea1) & params_data$ma==paste(ma1),]
# sets asymmetry
for(i in 1:dim(params_data)[1]){
  params_data$asym[i]=list(list(rep(params_data$gamma[i],nSpecies), # patch with modified attack rate
                                rep(1,nSpecies))) # reference patch
}
# sets dispersal
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}

# biomasses computed analytically (for d>>1, B_21=B_22)
B0<-equilibrium_analytical_TL2_disp_pred(params_data)

# integration time and pulse perturbations depending on parameters
t_max=10
t_step=0.1
pert_factor=1.2

# simulations
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1],.packages=c("deSolve")) %dopar% time_series_perturbation(params_data[i,], B0[i,], t_max, t_step, pert_factor,nSpecies, nCommunity)
stopCluster(cl)

TS<-create_TS(params_data, results, B0)
TS<-TS[,which(names(TS)%in%c("time","model","gamma",colnames(B0)))]
TS<-time_series_for_plot(TS, 2)

levels(TS$model)<-c("prey perturbed in patch #1 (fast patch)","prey perturbed in patch #2 (slow patch)")
#TS$gamma<-as.factor(TS$gamma)
#levels(TS$gamma)<-c(expression(italic("\u03B3")*"=3"))

p1<-ggplot(data=TS)+
      geom_line(aes(time,biomass,colour=species,linetype=community),size=1.5)+
      facet_wrap(~model)+
      corr_colour_TL_2+
      patch_line+
      theme+
      scale_x_continuous(breaks=seq(0,t_max,2))+
      scale_y_continuous(breaks=seq(0,pert_factor,0.1))+
      xlab("Time")+
      ylab("Scaled biomass")
ggsave(paste(path_figure,"figure_TS.pdf",sep=""), p1, width = 11, height = 5, device = cairo_pdf)

########################## ----
# APPENDIX ############## ----
########################## ----
#### DESCRIPTION OF PARAMETERS #### ----
nSpecies=4
nCommunity=1
params_data<-params_data_original
biomass<-cbind(params_data,as.data.frame(matrix(0,dim(params_data)[1],nSpecies)))
names(biomass)[(dim(params_data)[2]+1):dim(biomass)[2]]<-c(1:nSpecies)

for(i in 1:dim(biomass)[1]){
  params<-as.list(c(g=params_data$g[i],
                    r=params_data$r[i],
                    D=params_data$D[i],
                    m=params_data$m[i],
                    a=params_data$a[i],
                    e=params_data$e[i]))
  biomass[i,(dim(params_data)[2]+1):dim(biomass)[2]]<-equilibrium_symmetric(params,nSpecies,nCommunity)
}

biomass<-melt(biomass[,which(names(biomass)%in%c("ea","ma",names(biomass)[(dim(params_data)[2]+1):dim(biomass)[2]]))],
              id.vars = c("ea","ma"),
              variable.name = "species",
              value.name = "biomass")
biomass$x<-as.numeric(biomass$species)

p1<-ggplot(data=biomass)+
      geom_rect(aes(xmin=x-0.5,xmax=x+0.5,ymin=1e-4,ymax=biomass,fill=species))+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_4_fill+
      y_axis_log10_short+
      coord_flip()+
      theme+
      xlab('Trophic level')+
      ylab("Biomass")

params_data$x=0.5
params_data$y=0.5

p2<-ggplot(data=params_data)+
      geom_raster(aes(x,y,fill=m))+
      facet_grid(ma~ea, labeller=label_parsed)+
      scale_fill_viridis(name=expression(paste("Metabolic\nrate ratio",italic(m))),
                         breaks = unique(params_data$m),
                         trans = "log10")+
      theme_raster+theme(axis.text=element_blank(),
                         axis.ticks=element_blank(),
                         legend.key=element_blank(),
                         legend.text.align = 1)+
      scale_x_continuous(trans=log10_trans())+
      scale_y_continuous(trans=log10_trans())+
      xlab("")+
      ylab("")

graph<-ggdraw(xlim = c(0, 2), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(p2, 1, 0, 1, 1)+
  draw_plot_label(c("A","B"), c(0,0.98), c(1,1), size = 30)
ggsave(paste(path_figure,"supp_parameters.pdf",sep=""), graph, width = 18, height = 8, device = cairo_pdf)

# effect of attack rate a on biomass # ----
nSpecies=2
nCommunity=1
params_data<-params_data_original[params_data_original$ea==paste(ea01) & params_data_original$ma==paste(ma01),]
a_interval<-data.frame(a_interval=seq(0.1,20,0.01))
params_data<-merge(params_data,a_interval)
params_data$a<-params_data$a_interval
params_data$a_interval<-NULL
rm(a_interval)

biomass<-cbind(params_data,as.data.frame(matrix(0,dim(params_data)[1],nSpecies)))
names(biomass)[(dim(params_data)[2]+1):dim(biomass)[2]]<-c(1:nSpecies)

for(i in 1:dim(biomass)[1]){
  params<-as.list(c(g=params_data$g[i],
                    r=params_data$r[i],
                    D=params_data$D[i],
                    m=params_data$m[i],
                    a=params_data$a[i],
                    e=params_data$e[i]))
  biomass[i,(dim(params_data)[2]+1):dim(biomass)[2]]<-equilibrium_symmetric(params,nSpecies,nCommunity)
}
biomass<-melt(biomass[,which(names(biomass)%in%c("a",names(biomass)[(dim(params_data)[2]+1):dim(biomass)[2]]))],
              id.vars = c("a"),
              variable.name = "species",
              value.name = "biomass")

p1<-ggplot(data=biomass)+
  geom_line(aes(a,biomass,colour=species),size=1.5)+
  corr_colour_TL_2+
  theme+
  x_axis_log10+
  xlab(expression("Interaction strength "*italic(a)))+
  ylab("Biomass")

ggsave(paste(path_figure,"supp_biomass_attack_rate.pdf",sep=""),p1, width = 6, height = 4, device = cairo_pdf)

#### DISPERSAL OF PREDATORS #### ----
# coexistence # ----
gamma=seq(1,10,0.1) # asymmetry coefficient
omega=seq(1,10,0.1) # asymmetry in growth rate

params_data<-expand.grid(gamma=gamma,
                         omega=omega)
params_data<-merge(params_data_original,params_data)

data<-params_data
data$lambda<-data$e*data$a^2*data$m
# criterion on B1
data$crit_patch_1=data$omega*(2+data$lambda)/data$lambda
# criteria on B1'
data$crit_delta=sqrt(8/data$lambda)
data$crit_gamma_1=(data$lambda*data$omega+sqrt(data$lambda*(data$lambda*data$omega^2-8)))/2/data$lambda
data$crit_gamma_2=(data$lambda*data$omega-sqrt(data$lambda*(data$lambda*data$omega^2-8)))/2/data$lambda
data$coexistence=0

# biomasses
data$B_11=data$g/data$D*(2*data$omega+data$lambda*data$omega*(data$gamma^2+1)-data$lambda*data$gamma*(1+data$omega*data$gamma))/(2+data$lambda*(data$gamma^2+1))
data$B_12=data$g/data$D*(2+data$lambda*(data$gamma^2+1)-data$lambda*(1+data$omega*data$gamma))/(2+data$lambda*(data$gamma^2+1))
data$B_2=2*data$e*data$a*data$g/data$D*(1+data$omega*data$gamma)/(2+data$lambda*(data$gamma^2+1))/2
data$B_11[data$B_11<1e-3]=NA
data$B_12[data$B_12<1e-3]=NA

for(i in 1:dim(data)[1]){
  if(data$gamma[i]<data$crit_patch_1[i] & (data$omega[i]<data$crit_delta[i]
                                           | data$omega[i]>data$crit_delta[i] &
                                           (data$gamma[i]>data$crit_gamma_1[i] | data$gamma[i]<data$crit_gamma_2[i]))){
    data$coexistence[i]=1
  }
}
data$coexistence<-as.factor(data$coexistence)
levels(data$coexistence)<-c("no coexistence","coexistence")

p1<-ggplot(data=data)+
      geom_raster(aes(omega,gamma,fill=B_11))+
      facet_grid(ma~ea, labeller=label_parsed)+
      scale_fill_viridis(name=expression("B"[1]^(1)),
                         trans = "log10")+
      theme_raster+
      xlab(label_omega)+
      ylab(label_gamma)

p2<-ggplot(data=data)+
      geom_raster(aes(omega,gamma,fill=B_12))+
      facet_grid(ma~ea, labeller=label_parsed)+
      scale_fill_viridis(name=expression("B"[1]^(2)),
                         trans = "log10")+
      theme_raster+
      xlab(label_omega)+
      ylab(label_gamma)

p3<-ggplot(data=data)+
      geom_raster(aes(omega,gamma,fill=B_2))+
      facet_grid(ma~ea, labeller=label_parsed)+
      scale_fill_viridis(name=expression("B"[2]),
                         trans = "log10")+
      theme_raster+
      xlab(label_omega)+
      ylab(label_gamma)

p4<-ggplot(data=data)+
      geom_raster(aes(omega,gamma,fill=coexistence))+
      facet_grid(ma~ea, labeller=label_parsed)+
      scale_fill_manual(values = c("lightgrey","black"))+
      theme_raster+theme(legend.title = element_blank())+
      xlab(label_omega)+
      ylab(label_gamma)

graph<-ggdraw(xlim = c(0, 2.15), ylim = c(0, 2)) +
  draw_plot(p1, 0, 1, 1, 1)+
  draw_plot(p2, 1, 1, 1, 1)+
  draw_plot(p3, 0, 0, 1, 1)+
  draw_plot(p4, 1, 0, 1.15, 1)+
  draw_plot_label(c("A","B","C","D"), c(0,1,0,1), c(2,2,1,1), size = 30)
ggsave(paste(path_figure,"supp_coexistence_disp2.pdf",sep=""),graph, width = 20, height = 16, device = cairo_pdf)

# parameters and simulations - gamma - omega - perturbation of prey # ----
nSpecies=2
nCommunity=2
gamma=seq(1,10,0.1) # asymmetry coefficient
#omega=seq(1,10,0.1) # asymmetry in growth rate
d=1e6 # high dispersal

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
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)

B<-equilibrium_analytical_TL2_disp_pred(params_data) # biomass at equilibrium
coexistence<-coexistence_TL2_disp_pred(params_data) # coexistence of all species

# index of the parameters leading to significant coexistence
index<-which(coexistence$coexistence_significative==1)
params_data<-params_data[index,]
B<-B[index,] # selects the biomasses allowing coexistence

# simulations
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1],.packages=NULL) %dopar% analysis(params_data[i,], B[i,], nSpecies, nCommunity)
stopCluster(cl)
list_data<-create_data(params_data, results, nSpecies, nCommunity)

## correlation intra-patch # ----
data_C<-list_data$data_C
databis<-data_C[which(params_data$ea==paste(ea1) & params_data$ma==paste(ma1)),
                which(names(data_C)%in%c("ea","ma","d","gamma","model","C_11_21",'C_12_22'))]
databis<-melt(databis,
              id.vars = c("ea","ma","d","gamma","model"),
              variable.name = "community",
              value.name = "correlation")
levels(databis$community)<-c(1,2)
levels(databis$model)<-c(pert_11,pert_12)

p1<-ggplot(data=databis)+
      geom_line(aes(gamma,correlation,linetype=community),size=1.5)+
      geom_vline(xintercept=7,linetype="dashed",colour="red")+
      geom_label(data=data.frame(x=rep(7,2),y=rep(0.8,2),label=c("B","C"),model=levels(databis$model)),aes(x,y,label=label),size=8)+
      facet_wrap(~model)+
      patch_line+
      perturbation_prey_colour+
      theme+
      x_axis_gamma+
      ylim(-1,1)+
      xlab(label_gamma)+
      ylab(label_correlation_intra)

databis<-data_C[which(data_C$model%in%c("pert_11","pert_12")
                      & data_C$ea==paste(ea1) & data_C$ma==paste(ma1)
                      & data_C$gamma==7),
                which(names(data_C)%in%c("model","C_11_21","C_12_22","C_11_12",'C_21_22'))]
databis<-melt(databis,
              id.vars = c("model"),
              variable.name = "couple",
              value.name = "correlation")

p2<-ggplot(data=databis)+
      geom_point(aes(couple,correlation,colour=correlation),size=8)+
      facet_wrap(~model)+
      scale_colour_gradient2(low = "red",
                             mid = "white",
                             high = "blue",
                             midpoint = 0,
                             limits = c(-1,1))+
      ylim(-1,1)

graph<-ggdraw(xlim = c(0, 2), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(p2, 1, 0, 1, 1)+
  draw_plot_label(c("A","B"), c(0,1), c(1,1), size = 30)
ggsave(paste(path_figure,"supp_correlation_intra.pdf",sep=""), graph, width = 17, height = 6, device = cairo_pdf)


## correlation inter-patch # ----
data_C<-list_data$data_C
databis<-data_C[which(names(data_C)%in%c("ea","ma","d","gamma","model","C_11_12",'C_21_22'))]
databis<-table_for_plot(databis,5,"correlation")

p1<-ggplot(data=databis)+
      geom_line(aes(gamma,correlation,colour=species,linetype=model),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed)+
      perturbation_prey_line+
      corr_colour_TL_2+
      theme+
      x_axis_gamma+
      ylim(-1,1)+
      xlab(label_gamma)+
      ylab(label_correlation)
ggsave(paste(path_figure,"supp_correlation_disp2.pdf",sep=""),p1, width = 12, height = 8, device = cairo_pdf)

## biomass CV # ----
data_CV_species<-list_data$data_CV_species
data_CV_species<-data_CV_species[,which(names(data_CV_species)%in%c("ea","ma","d","gamma","omega","model",list_data$CV_species_names))]
data_CV_species<-table_for_plot(data_CV_species,6,"CV")
data_CV_species$species=data_CV_species$community
data_CV_species$community=NULL

p1<-ggplot(data=data_CV_species)+
      geom_line(aes(gamma,CV,colour=species,linetype=model),size=1.5)+
      corr_colour_TL_2+
      perturbation_prey_line+
      theme
legend<-get_legend(p1)

p1<-ggplot(data=data_CV_species)+
      geom_line(aes(gamma,CV,colour=species,linetype=model),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed, scales="free")+
      corr_colour_TL_2+
      perturbation_prey_line+
      theme+theme(legend.position="none")+
      x_axis_gamma+
      y_axis_log10+
      xlab(label_gamma)+
      ylab(label_CV_metapop)

data_CV<-list_data$data_CV

p2<-ggplot(data=data_CV)+
      geom_line(aes(gamma,CV_tot,linetype=model),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed, scales="free")+
      corr_colour_TL_2+
      perturbation_prey_line+
      theme+theme(legend.position="none")+
      x_axis_gamma+
      y_axis_log10+
      xlab(label_gamma)+
      ylab(label_CV_metacom)

graph<-ggdraw(xlim = c(0, 2.3), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(p2, 1, 0, 1, 1)+
  draw_plot(legend, 2.1, 0.25, 0.1, 0.5)+
  draw_plot_label(c("A","B"), c(0,1), c(1,1), size = 30)
ggsave(paste(path_figure,"supp_CV_disp2_pert1.pdf",sep=""), graph, width = 16, height = 7, device = cairo_pdf)

## biomass # ----
data_B<-list_data$data_B
data_B<-data_B[data_B$model=="pert_11",
               which(names(data_B)%in%c("ea","ma","d","gamma","omega","model",list_data$B_names))]
data_B<-table_for_plot(data_B,6,"biomass")
data_B$biomass[data_B$biomass<1e-6]=NA

p1<-ggplot(data=data_B)+
  geom_line(aes(gamma,biomass,colour=species,linetype=community),size=1.5)+
  corr_colour_TL_2+
  patch_line+
  theme
legend<-get_legend(p1)

p1<-ggplot(data=data_B)+
      geom_line(aes(gamma,biomass,colour=species,linetype=community),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed, scales="free")+
      corr_colour_TL_2+
      patch_line+
      theme+theme(legend.position="none") +
      x_axis_gamma+
      xlab(label_gamma)+
      ylab("Biomass")

## biomass in isolated patches # ----
biomass<-cbind(params_data,as.data.frame(matrix(0,dim(params_data)[1],nSpecies*nCommunity)))
names(biomass)[(dim(params_data)[2]+1):dim(biomass)[2]]<-c(1:(nSpecies*nCommunity))

for(i in 1:dim(biomass)[1]){
  params<-as.list(c(g=params_data$g[i],
                    r=params_data$r[i],
                    D=params_data$D[i],
                    m=params_data$m[i],
                    a=params_data$a[i],
                    e=params_data$e[i]))
  params$asym=params_data$asym[[i]]
  params$disp=params_data$disp[[i]]
  biomass[i,c(1,2)+dim(params_data)[2]]<-equilibrium_isolated(params,nSpecies,1)
  biomass[i,c(3,4)+dim(params_data)[2]]<-equilibrium_isolated(params,nSpecies,2)
}
biomass<-biomass[biomass$model=="pert_11",]

# biomasses with dispersal
data_B<-list_data$data_B
data_B<-data_B[data_B$model=="pert_11",
               which(names(data_B)%in%c("ea","ma","d","gamma","omega","model",list_data$B_names))]
data_B[,dim(data_B)[2]-(3:0)]<-data_B[,dim(data_B)[2]-(3:0)]/biomass[,dim(biomass)[2]-(3:0)]
data_B<-table_for_plot(data_B,6,"biomass")
data_B$biomass[data_B$biomass<1e-6]=NA

p2<-ggplot(data=data_B)+
      geom_line(aes(gamma,biomass,colour=species,linetype=community),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed, scales="free")+
      corr_colour_TL_2+
      patch_line+
      theme+theme(legend.position="none") +
      x_axis_gamma+
      xlab(label_gamma)+
      ylab("Scaled biomass")

# final figure # ----
graph<-ggdraw(xlim = c(0, 2.15), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(p2, 1, 0, 1, 1)+
  draw_plot(legend, 2.05, 0.25, 0.05, 0.5)+
  draw_plot_label(c("A","B"), c(0,0.99), c(1,1), size = 30)
ggsave(paste(path_figure,"supp_biomass_disp2.pdf",sep=""), graph, width = 18, height = 8, device = cairo_pdf)

## asymptotic resilience # ----
data_resilience<-list_data$data_resilience # lead eigenvalue

p1<-ggplot(data=data_resilience[data_resilience$model=="pert_11",])+
      geom_line(aes(gamma,resilience),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed, scales="free")+
      fill_resilience+
      theme+
      x_axis_gamma+
      xlab(label_gamma)+
      ylab(label_resilience)

## lead eigenvector # ----
data_E<-list_data$data_E # lead eigenvector
data_E<-data_E[data_E$model=="pert_11",
               which(names(data_E)%in%c("model","gamma","omega","ea","ma",list_data$E_names))]
data_E<-table_for_plot(data_E,5,"contribution")

p2<-ggplot(data=data_E)+
  geom_line(aes(gamma,contribution,colour=species,linetype=community),size=1.5)+
  corr_colour_TL_2+
  patch_line+
  theme
legend<-get_legend(p2)

p2<-ggplot(data=data_E)+
      geom_line(aes(gamma,contribution,colour=species,linetype=community),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed, scales="free")+
      corr_colour_TL_2+
      patch_line+
      theme+theme(legend.position="none")+
      x_axis_gamma+
      ylim(0,1)+
      xlab(label_gamma)+
      ylab(label_contribution_2lines)

# final figure # ----
graph<-ggdraw(xlim = c(0, 2.15), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(p2, 1, 0, 1, 1)+
  draw_plot(legend, 2.05, 0.25, 0.05, 0.5)+
  draw_plot_label(c("A","B"), c(0,1), c(1,1), size = 30)
ggsave(paste(path_figure,"supp_eigen_disp2.pdf",sep=""), graph, width = 18, height = 9, device = cairo_pdf)

## time series gamma=10 # ----
nSpecies=2
nCommunity=2
gamma=c(1,10) # asymmetry coefficient
d=1e6

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
params_data<-merge(params_data_original,params_data)
params_data$omega=params_data$gamma
params_data<-params_data[params_data$ea==paste(ea1) & params_data$ma==paste(ma1),]
# sets asymmetry
for(i in 1:dim(params_data)[1]){
  params_data$asym[i]=list(list(rep(params_data$gamma[i],nSpecies), # patch with modified attack rate
                                rep(1,nSpecies))) # reference patch
}
# sets dispersal
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}

# biomasses computed analytically (for d>>1, B_21=B_22)
B0<-equilibrium_analytical_TL2_disp_pred(params_data)

# integration time and pulse perturbations depending on parameters
t_max=10
t_step=0.1
pert_factor=1.2

# simulations
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1],.packages=c("deSolve")) %dopar% time_series_perturbation(params_data[i,], B0[i,], t_max, t_step, pert_factor,nSpecies, nCommunity)
stopCluster(cl)

TS<-create_TS(params_data, results, B0)
TS<-TS[,which(names(TS)%in%c("time","model","gamma",colnames(B0)))]
TS<-time_series_for_plot(TS, 2)

levels(TS$model)<-c(expression("prey perturbed"*" in patch #1"),expression("prey perturbed"*" in patch #2"))
TS$gamma<-as.factor(TS$gamma)
levels(TS$gamma)<-c(expression(italic("\u03B3")*"=1"),expression(italic("\u03B3")*"=10"))

p1<-ggplot(data=TS)+
      geom_line(aes(time,biomass,colour=species,linetype=community),size=1.5)+
      facet_grid(gamma~model,scale="free", labeller=label_parsed)+
      corr_colour_TL_2+
      patch_line+
      theme+
      scale_x_continuous(breaks=seq(0,t_max,2))+
      scale_y_continuous(breaks=seq(0,pert_factor,0.1))+
      xlab("Time")+
      ylab("Scaled biomass")

ggsave(paste(path_figure,"supp_TS_disp2_gamma10.pdf",sep=""), p1, width = 8, height = 5, device = cairo_pdf)

## time series # ----
nSpecies=2
nCommunity=2
gamma=3 # asymmetry coefficient
d=1e6

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
params_data<-merge(params_data_original,params_data)
params_data$omega=params_data$gamma
# sets asymmetry
for(i in 1:dim(params_data)[1]){
  params_data$asym[i]=list(list(rep(params_data$gamma[i],nSpecies), # patch with modified attack rate
                                rep(1,nSpecies))) # reference patch
}
# sets dispersal
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}

# biomasses computed analytically (for d>>1, B_21=B_22)
B0<-equilibrium_analytical_TL2_disp_pred(params_data)

# integration time and pulse perturbations depending on parameters
t_max=10
t_step=0.1
pert_factor=1.2

# simulations
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1],.packages=c("deSolve")) %dopar% time_series_perturbation(params_data[i,], B0[i,], t_max, t_step, pert_factor,nSpecies, nCommunity)
stopCluster(cl)

TS<-create_TS(params_data, results, B0)
TS<-TS[,which(names(TS)%in%c("time","model","gamma","ea","ma",colnames(B0)))]
TS<-time_series_for_plot(TS, 4)

levels(TS$model)<-c(pert_11,pert_12)

p1<-ggplot(data=TS[TS$model==pert_11,])+
  geom_line(aes(time,biomass,colour=species,linetype=community),size=1.5)+
  corr_colour_TL_2+
  patch_line+
  theme
legend<-get_legend(p1)

p1<-ggplot(data=TS[TS$model==pert_11,])+
      geom_line(aes(time,biomass,colour=species,linetype=community),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed, scales="free")+
      corr_colour_TL_2+
      patch_line+
      theme+theme(legend.position="none")+
      scale_x_continuous(breaks=seq(0,t_max,2))+
      scale_y_continuous(breaks=seq(0,pert_factor,0.1))+
      xlab("Time")+
      ylab("Scaled biomass")+
      ggtitle("prey perturbed in patch #1")

p2<-ggplot(data=TS[TS$model==pert_12,])+
      geom_line(aes(time,biomass,colour=species,linetype=community),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed, scales="free")+
      corr_colour_TL_2+
      patch_line+
      theme+theme(legend.position="none")+
      scale_x_continuous(breaks=seq(0,t_max,2))+
      scale_y_continuous(breaks=seq(0,pert_factor,0.1))+
      xlab("Time")+
      ylab("Scaled biomass")+
      ggtitle("prey perturbed in patch #2")

graph<-ggdraw(xlim = c(0, 2.15), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(p2, 1, 0, 1, 1)+
  draw_plot(legend, 2.05, 0.25, 0.05, 0.5)+
  draw_plot_label(c("A","B"), c(0,1), c(1,1), size = 30)
ggsave(paste(path_figure,"supp_TS_disp2.pdf",sep=""), graph, width = 14, height = 7, device = cairo_pdf)

# parameters and simulations - gamma - omega - effect of omega # ----
nSpecies=2
nCommunity=2
gamma=seq(1,10,0.1) # asymmetry coefficient
omega=seq(1,10,3) # asymmetry in growth rate
d=1e6 # high dispersal

params_data<-expand.grid(pert=list(list(c(1,1))), # perturbation of prey in patch 1
                         disp=list(c(0,1)), # dispersal of predators
                         gamma=gamma,
                         omega=omega,
                         d=d,
                         model="pert_11")
params_data<-rbind(params_data,expand.grid(pert=list(list(c(1,2))), # perturbation of prey in patch 2
                                           disp=list(c(0,1)), # dispersal of predators
                                           gamma=gamma,
                                           omega=omega,
                                           d=d,
                                           model="pert_12"))
params_data<-merge(params_data_original,params_data)
params_data<-params_data[params_data$ea==paste(ea01) & params_data$ma==paste(ma01),]
# sets asymmetry
for(i in 1:dim(params_data)[1]){
  params_data$asym[i]=set_asymmetry(nSpecies,params_data$gamma[i],params_data$omega[i])
}
# sets dispersal
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)

B<-equilibrium_analytical_TL2_disp_pred(params_data) # biomass at equilibrium
coexistence<-coexistence_TL2_disp_pred(params_data) # coexistence of all species

# index of the parameters leading to significant coexistence
index<-which(coexistence$coexistence_significative==1)
params_data<-params_data[index,]
B<-B[index,] # selects the biomasses allowing coexistence

# simulations
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1],.packages=NULL) %dopar% analysis(params_data[i,], B[i,], nSpecies, nCommunity)
stopCluster(cl)
list_data<-create_data(params_data, results, nSpecies, nCommunity)

## biomass # ----
data_B<-list_data$data_B
data_B<-data_B[data_B$model=="pert_11",
               which(names(data_B)%in%c("ea","ma","d","gamma","omega","model",list_data$B_names))]
data_B<-table_for_plot(data_B,6,"biomass")
data_B$biomass[data_B$biomass<1e-6]=NA
data_B$omega=as.factor(data_B$omega)

p1<-ggplot(data=data_B)+
  geom_line(aes(gamma,biomass,colour=species,linetype=community,alpha=omega),size=1.5)+
  corr_colour_TL_2+
  patch_line+
  omega_alpha+
  theme
legend<-get_legend(p1)

p1<-ggplot(data=data_B)+
  geom_line(aes(gamma,biomass,colour=species,linetype=community,alpha=omega),size=1.5)+
  corr_colour_TL_2+
  patch_line+
  omega_alpha+
  theme+theme(legend.position="none")+
  x_axis_gamma+
  xlab(label_gamma)+
  ylab("Biomass")

## biomass in isolated patches # ----
biomass<-cbind(params_data,as.data.frame(matrix(0,dim(params_data)[1],nSpecies*nCommunity)))
names(biomass)[(dim(params_data)[2]+1):dim(biomass)[2]]<-c(1:(nSpecies*nCommunity))

for(i in 1:dim(biomass)[1]){
  params<-as.list(c(g=params_data$g[i],
                    r=params_data$r[i],
                    D=params_data$D[i],
                    m=params_data$m[i],
                    a=params_data$a[i],
                    e=params_data$e[i]))
  params$asym=params_data$asym[[i]]
  params$disp=params_data$disp[[i]]
  biomass[i,c(1,2)+dim(params_data)[2]]<-equilibrium_isolated(params,nSpecies,1)
  biomass[i,c(3,4)+dim(params_data)[2]]<-equilibrium_isolated(params,nSpecies,2)
}
biomass<-biomass[biomass$model=="pert_11",]

# biomasses with dispersal
data_B<-list_data$data_B
data_B<-data_B[data_B$model=="pert_11",
               which(names(data_B)%in%c("ea","ma","d","gamma","omega","model",list_data$B_names))]
data_B[,dim(data_B)[2]-(3:0)]<-data_B[,dim(data_B)[2]-(3:0)]/biomass[,dim(biomass)[2]-(3:0)]
data_B<-table_for_plot(data_B,6,"biomass")
data_B$biomass[data_B$biomass<1e-6]=NA
data_B$omega=as.factor(data_B$omega)

p2<-ggplot(data=data_B)+
  geom_line(aes(gamma,biomass,colour=species,linetype=community,alpha=omega),size=1.5)+
  corr_colour_TL_2+
  patch_line+
  omega_alpha+
  theme+theme(legend.position="none")+
  x_axis_gamma+
  xlab(label_gamma)+
  ylab("Scaled biomass")

# final figure # ----
graph<-ggdraw(xlim = c(0, 2.5), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(p2, 1, 0, 1, 1)+
  draw_plot(legend, 2.2, 0.25, 0.05, 0.5)+
  draw_plot_label(c("A","B"), c(0,0.99), c(1,1), size = 30)
ggsave(paste(path_figure,"supp_biomass_disp2_omega.pdf",sep=""), graph, width = 14, height = 6, device = cairo_pdf)

## correlation and CV # ----
data_C<-list_data$data_C
databis<-data_C[which(names(data_C)%in%c("ea","ma","d","gamma","omega","model","C_11_12",'C_21_22'))]
databis<-table_for_plot(databis,6,"correlation")
databis$omega=as.factor(databis$omega)

p1<-ggplot(data=databis)+
  geom_line(aes(gamma,correlation,colour=species,linetype=model,alpha=omega),size=1.5)+
  perturbation_prey_line+
  corr_colour_TL_2+
  omega_alpha+
  theme
legend<-get_legend(p1)

p1<-ggplot(data=databis)+
  geom_line(aes(gamma,correlation,colour=species,linetype=model,alpha=omega),size=1.5)+
  perturbation_prey_line+
  corr_colour_TL_2+
  omega_alpha+
  theme+theme(legend.position="none")+
  x_axis_gamma+
  ylim(-1,1)+
  xlab(label_gamma)+
  ylab(label_correlation)

data_CV_species<-list_data$data_CV_species
data_CV_species<-data_CV_species[,which(names(data_CV_species)%in%c("ea","ma","d","gamma","omega","model",list_data$CV_species_names))]
data_CV_species<-table_for_plot(data_CV_species,6,"CV")
data_CV_species$species=data_CV_species$community
data_CV_species$community=NULL
data_CV_species$omega=as.factor(data_CV_species$omega)

p2<-ggplot(data=data_CV_species)+
  geom_line(aes(gamma,CV,colour=species,linetype=model,alpha=omega),size=1.5)+
  corr_colour_TL_2+
  perturbation_prey_line+
  omega_alpha+
  theme+theme(legend.position="none")+
  x_axis_gamma+
  y_axis_log10+
  xlab(label_gamma)+
  ylab(label_CV_metapop)

graph<-ggdraw(xlim = c(0, 2.5), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(p2, 1, 0, 1, 1)+
  draw_plot(legend, 2.2, 0.25, 0.05, 0.5)+
  draw_plot_label(c("A","B"), c(0,0.99), c(1,1), size = 30)
ggsave(paste(path_figure,"supp_correlation_CV_disp2_omega.pdf",sep=""), graph, width = 14, height = 6.5, device = cairo_pdf)

# parameters and simulations - gamma - omega - perturbation of predators # ----
nSpecies=2
nCommunity=2
gamma=seq(1,10,0.1) # asymmetry coefficient
#omega=seq(1,10,0.1) # asymmetry in growth rate
d=1e6 # high dispersal

params_data<-expand.grid(pert=list(list(c(2,1))), # perturbation of predators in patch 1
                         disp=list(c(0,1)), # dispersal of predators
                         gamma=gamma,
                         d=d,
                         model="pert_21")
params_data<-rbind(params_data,expand.grid(pert=list(list(c(2,2))), # perturbation of predators in patch 2
                                           disp=list(c(0,1)), # dispersal of predators
                                           gamma=gamma,
                                           d=d,
                                           model="pert_22"))
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
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)

B<-equilibrium_analytical_TL2_disp_pred(params_data) # biomass at equilibrium
coexistence<-coexistence_TL2_disp_pred(params_data) # coexistence of all species

# index of the parameters leading to significant coexistence
index<-which(coexistence$coexistence_significative==1)
params_data<-params_data[index,]
B<-B[index,] # selects the biomasses allowing coexistence

# simulations
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1],.packages=NULL) %dopar% analysis(params_data[i,], B[i,], nSpecies, nCommunity)
stopCluster(cl)
list_data<-create_data(params_data, results, nSpecies, nCommunity)

## correlation inter-patch #### ----
data_C<-list_data$data_C
databis<-data_C[,which(names(data_C)%in%c("ea","ma","d","gamma","model","C_11_12",'C_21_22'))]
databis<-table_for_plot(databis,5,"correlation")

p1<-ggplot(data=databis)+
      geom_line(aes(gamma,correlation,colour=species,linetype=model),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed)+
      perturbation_pred_line+
      corr_colour_TL_2+
      theme+
      x_axis_gamma+
      #ylim(-1,1)+
      xlab(label_gamma)+
      ylab(label_correlation)

ggsave(paste(path_figure,"supp_correlation_inter_disp2_pert2.pdf",sep=""),p1, width = 12, height = 8, device = cairo_pdf)

## biomass CV # ----
data_CV_species<-list_data$data_CV_species
data_CV_species<-data_CV_species[,which(names(data_CV_species)%in%c("ea","ma","d","gamma","omega","model",list_data$CV_species_names))]
data_CV_species<-table_for_plot(data_CV_species,6,"CV")
data_CV_species$species=data_CV_species$community
data_CV_species$community=NULL

p1<-ggplot(data=data_CV_species)+
  geom_line(aes(gamma,CV,colour=species,linetype=model),size=1.5)+
  corr_colour_TL_2+
  perturbation_prey_line+
  theme
legend<-get_legend(p1)

p1<-ggplot(data=data_CV_species)+
      geom_line(aes(gamma,CV,colour=species,linetype=model),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed, scales="free")+
      corr_colour_TL_2+
      perturbation_prey_line+
      theme+theme(legend.position="none")+
      x_axis_gamma+
      y_axis_log10+
      xlab(label_gamma)+
      ylab(label_CV_metapop)

data_CV<-list_data$data_CV

p2<-ggplot(data=data_CV)+
      geom_line(aes(gamma,CV_tot,linetype=model),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed, scales="free")+
      corr_colour_TL_2+
      perturbation_prey_line+
      theme+theme(legend.position="none")+
      x_axis_gamma+
      y_axis_log10+
      xlab(label_gamma)+
      ylab(label_CV_metacom)

graph<-ggdraw(xlim = c(0, 2.3), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(p2, 1, 0, 1, 1)+
  draw_plot(legend, 2.1, 0.25, 0.1, 0.5)+
  draw_plot_label(c("A","B"), c(0,1), c(1,1), size = 30)
ggsave(paste(path_figure,"supp_CV_disp2_pert2.pdf",sep=""), graph, width = 16, height = 7, device = cairo_pdf)

# parameters and simulations - TL3 - gamma - omega - perturbation of prey # ----
nSpecies=3
nCommunity=2
gamma=seq(1,10,0.1) # asymmetry coefficient
omega=1 # asymmetry in growth rate
d=1e6 # high dispersal

params_data<-expand.grid(pert=list(list(c(1,1))), # perturbation of prey in patch 1
                         disp=list(c(0,0,1)), # dispersal of top predators
                         gamma=gamma,
                         d=d,
                         model="pert_11")
params_data<-rbind(params_data,expand.grid(pert=list(list(c(1,2))), # perturbation of prey in patch 2
                                           disp=list(c(0,0,1)), # dispersal of predators
                                           gamma=gamma,
                                           d=d,
                                           model="pert_12"))
params_data$omega=omega#params_data$gamma
params_data<-merge(params_data_original,params_data)
# sets asymmetry
for(i in 1:dim(params_data)[1]){
  params_data$asym[i]=set_asymmetry(nSpecies,params_data$gamma[i],params_data$omega[i])
}
# sets dispersal
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)

B<-matrix(0,dim(params_data)[1],nSpecies*nCommunity)
for(i in 1:dim(params_data)[1]){
  params<-get_params(params_data,i)
  B[i,]<-equilibrium_disp_pred(params,nSpecies,nCommunity)
}
# index of the parameters leading to significant coexistence
index<-which(B<0,arr.ind=T)[,1]
if(length(index)>0){
  B<-B[-index,] # selects the biomasses allowing coexistence
  params_data<-params_data[-index,]
}

# simulations
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1],.packages=NULL) %dopar% analysis(params_data[i,], B[i,], nSpecies, nCommunity)
stopCluster(cl)
list_data<-create_data(params_data, results, nSpecies, nCommunity)

## correlation inter-patch # ----
data_C<-list_data$data_C
databis<-data_C[,which(names(data_C)%in%c("ea","ma","d","gamma","model","C_11_12",'C_21_22','C_31_32'))]
databis<-table_for_plot(databis,5,"correlation")

p1<-ggplot(data=databis)+
      geom_line(aes(gamma,correlation,colour=species,linetype=model),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed)+
      perturbation_prey_line+
      corr_colour_TL_3+
      theme+
      x_axis_gamma+
      ylim(-1,1.01)+
      xlab(label_gamma)+
      ylab(label_correlation)

ggsave(paste(path_figure,"supp_correlation_inter_disp2_TL3.pdf",sep=""),p1, width = 12, height = 8, device = cairo_pdf)

## biomass CV # ----
data_CV_species<-list_data$data_CV_species
data_CV_species<-data_CV_species[,which(names(data_CV_species)%in%c("ea","ma","d","gamma","omega","model",list_data$CV_species_names))]
data_CV_species<-table_for_plot(data_CV_species,6,"CV")
data_CV_species$species=data_CV_species$community
data_CV_species$community=NULL

p1<-ggplot(data=data_CV_species)+
  geom_line(aes(gamma,CV,colour=species,linetype=model),size=1.5)+
  corr_colour_TL_3+
  perturbation_prey_line+
  theme
legend<-get_legend(p1)

p1<-ggplot(data=data_CV_species)+
  geom_line(aes(gamma,CV,colour=species,linetype=model),size=1.5)+
  facet_grid(ma~ea, labeller=label_parsed, scales="free")+
  corr_colour_TL_3+
  perturbation_prey_line+
  theme+theme(legend.position="none")+
  x_axis_gamma+
  y_axis_log10+
  xlab(label_gamma)+
  ylab(label_CV_metapop)

data_CV<-list_data$data_CV

p2<-ggplot(data=data_CV)+
  geom_line(aes(gamma,CV_tot,linetype=model),size=1.5)+
  facet_grid(ma~ea, labeller=label_parsed, scales="free")+
  corr_colour_TL_3+
  perturbation_prey_line+
  theme+theme(legend.position="none")+
  x_axis_gamma+
  y_axis_log10+
  xlab(label_gamma)+
  ylab(label_CV_metacom)

graph<-ggdraw(xlim = c(0, 2.3), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(p2, 1, 0, 1, 1)+
  draw_plot(legend, 2.1, 0.25, 0.1, 0.5)+
  draw_plot_label(c("A","B"), c(0,1), c(1,1), size = 30)
ggsave(paste(path_figure,"supp_CV_disp2_pert1_TL3.pdf",sep=""), graph, width = 16, height = 7, device = cairo_pdf)

#### DISPERSAL OF PREY #### ----
# coexistence # ----
gamma=seq(1,10,0.1) # asymmetry coefficient
omega=seq(1,10,0.1) # asymmetry in growth rate

params_data<-expand.grid(gamma=gamma,
                         omega=omega)
params_data<-merge(params_data_original,params_data)

data<-params_data
data$lambda<-data$e*data$a^2*data$m

data$B_1<-data$g/data$D*(1+data$omega)/(1+data$lambda/2*(data$gamma^2+1))
data$B_21<-data$g/data$D*data$e*data$a*data$gamma*(1+data$omega)/(2+data$lambda*(data$gamma^2+1))
data$B_22<-data$g/data$D*data$e*data$a*(1+data$omega)/(2+data$lambda*(data$gamma^2+1))

p1<-ggplot(data=data)+
  geom_raster(aes(omega,gamma,fill=B_1))+
  facet_grid(ma~ea, labeller=label_parsed)+
  scale_fill_viridis(name=expression("B"[1]),
                     trans = "log10")+
  theme_raster+
  xlab(label_omega)+
  ylab(label_gamma)

p2<-ggplot(data=data)+
  geom_raster(aes(omega,gamma,fill=B_21))+
  facet_grid(ma~ea, labeller=label_parsed)+
  scale_fill_viridis(name=expression("B"[2]^(1)),
                     trans = "log10")+
  theme_raster+
  xlab(label_omega)+
  ylab(label_gamma)

p3<-ggplot(data=data)+
  geom_raster(aes(omega,gamma,fill=B_22))+
  facet_grid(ma~ea, labeller=label_parsed)+
  scale_fill_viridis(name=expression("B"[2]^(2)),
                     trans = "log10")+
  theme_raster+
  xlab(label_omega)+
  ylab(label_gamma)

p4<-ggplot(data=data)+
  geom_raster(aes(omega,gamma,fill=B_21/B_22))+
  facet_grid(ma~ea, labeller=label_parsed)+
  scale_fill_viridis(name=expression("B"[2]^(1)*"/B"[2]^(2)))+
  theme_raster+
  xlab(label_omega)+
  ylab(label_gamma)

graph<-ggdraw(xlim = c(0, 2.15), ylim = c(0, 2)) +
  draw_plot(p1, 0, 1, 1, 1)+
  draw_plot(p2, 1, 1, 1, 1)+
  draw_plot(p3, 0, 0, 1, 1)+
  draw_plot(p4, 1, 0, 1, 1)+
  draw_plot_label(c("A","B","C","D"), c(0,1,0,1), c(2,2,1,1), size = 30)
ggsave(paste(path_figure,"supp_coexistence_disp1.pdf",sep=""),graph, width = 20, height = 16, device = cairo_pdf)

# parameters and simulations - gamma - omega - perturbation of predators # ----
nSpecies=2
nCommunity=2
gamma=seq(1,10,0.1) # asymmetry coefficient
#omega=seq(1,10,0.1) # asymmetry in growth rate
d=1e6 # high dispersal

params_data<-expand.grid(pert=list(list(c(2,1))), # perturbation of predators in patch 1
                         disp=list(c(1,0)), # dispersal of prey
                         gamma=gamma,
                         d=d,
                         model="pert_21")
params_data<-rbind(params_data,expand.grid(pert=list(list(c(2,2))), # perturbation of predators in patch 2
                                           disp=list(c(1,0)), # dispersal of prey
                                           gamma=gamma,
                                           d=d,
                                           model="pert_22"))
params_data$omega=1#params_data$gamma
params_data<-merge(params_data_original,params_data)
# sets asymmetry
for(i in 1:dim(params_data)[1]){
  params_data$asym[i]=set_asymmetry(nSpecies,params_data$gamma[i],params_data$omega[i])
}
# sets dispersal
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)

B<-matrix(0,dim(params_data)[1],nSpecies*nCommunity)
for(i in 1:dim(params_data)[1]){
  params<-get_params(params_data,i)
  B[i,]<-equilibrium_disp_prey(params,nSpecies,nCommunity)
}

# simulations
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1],.packages=NULL) %dopar% analysis(params_data[i,], B[i,], nSpecies, nCommunity)
stopCluster(cl)
list_data<-create_data(params_data, results, nSpecies, nCommunity)

## correlation inter-patch # ----
data_C<-list_data$data_C
databis<-data_C[which(names(data_C)%in%c("ea","ma","d","gamma","model","C_11_12",'C_21_22'))]
databis<-table_for_plot(databis,5,"correlation")

p1<-ggplot(data=databis)+
      geom_line(aes(gamma,correlation,colour=species,linetype=model),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed)+
      perturbation_pred_line+
      corr_colour_TL_2+
      theme+
      x_axis_gamma+
      ylim(-1,1)+
      xlab(label_gamma)+
      ylab(label_correlation)
ggsave(paste(path_figure,"supp_correlation_disp1.pdf",sep=""),p1, width = 12, height = 8, device = cairo_pdf)

## biomass CV # ----
data_CV_species<-list_data$data_CV_species
data_CV_species<-data_CV_species[,which(names(data_CV_species)%in%c("ea","ma","d","gamma","omega","model",list_data$CV_species_names))]
data_CV_species<-table_for_plot(data_CV_species,6,"CV")
data_CV_species$species=data_CV_species$community
data_CV_species$community=NULL

p1<-ggplot(data=data_CV_species)+
  geom_line(aes(gamma,CV,colour=species,linetype=model),size=1.5)+
  corr_colour_TL_2+
  perturbation_pred_line+
  theme
legend<-get_legend(p1)

p1<-ggplot(data=data_CV_species)+
      geom_line(aes(gamma,CV,colour=species,linetype=model),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed, scales="free")+
      corr_colour_TL_2+
      perturbation_pred_line+
      theme+theme(legend.position="none")+
      x_axis_gamma+
      y_axis_log10+
      xlab(label_gamma)+
      ylab(label_CV_metapop)

data_CV<-list_data$data_CV

p2<-ggplot(data=data_CV)+
      geom_line(aes(gamma,CV_tot,linetype=model),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed, scales="free")+
      corr_colour_TL_2+
      perturbation_pred_line+
      theme+theme(legend.position="none")+
      x_axis_gamma+
      y_axis_log10+
      xlab(label_gamma)+
      ylab(label_CV_metacom)

graph<-ggdraw(xlim = c(0, 2.3), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(p2, 1, 0, 1, 1)+
  draw_plot(legend, 2.1, 0.25, 0.1, 0.5)+
  draw_plot_label(c("A","B"), c(0,1), c(1,1), size = 30)
ggsave(paste(path_figure,"supp_CV_disp1_pert2.pdf",sep=""), graph, width = 16, height = 7, device = cairo_pdf)

## biomass # ----
data_B<-list_data$data_B
data_B<-data_B[data_B$model=="pert_21",
               which(names(data_B)%in%c("ea","ma","d","gamma","omega","model",list_data$B_names))]
data_B<-table_for_plot(data_B,6,"biomass")
data_B$biomass[data_B$biomass<1e-6]=NA

p1<-ggplot(data=data_B)+
  geom_line(aes(gamma,biomass,colour=species,linetype=community),size=1.5)+
  corr_colour_TL_2+
  patch_line+
  theme
legend<-get_legend(p1)

p1<-ggplot(data=data_B)+
      geom_line(aes(gamma,biomass,colour=species,linetype=community),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed, scales="free")+
      corr_colour_TL_2+
      patch_line+
      theme+theme(legend.position="none") +
      x_axis_gamma+
      xlab(label_gamma)+
      ylab("Biomass")

## biomass in isolated patches # ----
biomass<-cbind(params_data,as.data.frame(matrix(0,dim(params_data)[1],nSpecies*nCommunity)))
names(biomass)[(dim(params_data)[2]+1):dim(biomass)[2]]<-c(1:(nSpecies*nCommunity))

for(i in 1:dim(biomass)[1]){
  params<-as.list(c(g=params_data$g[i],
                    r=params_data$r[i],
                    D=params_data$D[i],
                    m=params_data$m[i],
                    a=params_data$a[i],
                    e=params_data$e[i]))
  params$asym=params_data$asym[[i]]
  params$disp=params_data$disp[[i]]
  biomass[i,c(1,2)+dim(params_data)[2]]<-equilibrium_isolated(params,nSpecies,1)
  biomass[i,c(3,4)+dim(params_data)[2]]<-equilibrium_isolated(params,nSpecies,2)
}
biomass<-biomass[biomass$model=="pert_21",]

# biomasses with dispersal
data_B<-list_data$data_B
data_B<-data_B[data_B$model=="pert_21",
               which(names(data_B)%in%c("ea","ma","d","gamma","omega","model",list_data$B_names))]
data_B[,dim(data_B)[2]-(3:0)]<-data_B[,dim(data_B)[2]-(3:0)]/biomass[,dim(biomass)[2]-(3:0)]
data_B<-table_for_plot(data_B,6,"biomass")
data_B$biomass[data_B$biomass<1e-6]=NA

p2<-ggplot(data=data_B)+
      geom_line(aes(gamma,biomass,colour=species,linetype=community),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed, scales="free")+
      corr_colour_TL_2+
      patch_line+
      theme+theme(legend.position="none") +
      x_axis_gamma+
      xlab(label_gamma)+
      ylab("Scaled biomass")

# final figure # ----
graph<-ggdraw(xlim = c(0, 2.15), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(p2, 1, 0, 1, 1)+
  draw_plot(legend, 2.05, 0.25, 0.05, 0.5)+
  draw_plot_label(c("A","B"), c(0,0.99), c(1,1), size = 30)
ggsave(paste(path_figure,"supp_biomass_disp1.pdf",sep=""), graph, width = 18, height = 8, device = cairo_pdf)

## asymptotic resilience # ----
data_resilience<-list_data$data_resilience # lead eigenvalue

p1<-ggplot(data=data_resilience[data_resilience$model=="pert_21",])+
      geom_line(aes(gamma,resilience),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed, scales="free")+
      fill_resilience+
      theme+
      x_axis_gamma+
      xlab(label_gamma)+
      ylab(label_resilience)

## lead eigenvector # ----
data_E<-list_data$data_E # lead eigenvector
data_E<-data_E[data_E$model=="pert_21",
               which(names(data_E)%in%c("model","gamma","omega","ea","ma",list_data$E_names))]
data_E<-table_for_plot(data_E,5,"contribution")

p2<-ggplot(data=data_E)+
  geom_line(aes(gamma,contribution,colour=species,linetype=community),size=1.5)+
  corr_colour_TL_2+
  patch_line+
  theme
legend<-get_legend(p2)

p2<-ggplot(data=data_E)+
      geom_line(aes(gamma,contribution,colour=species,linetype=community),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed, scales="free")+
      corr_colour_TL_2+
      patch_line+
      theme+theme(legend.position="none")+
      x_axis_gamma+
      ylim(0,1)+
      xlab(label_gamma)+
      ylab(label_contribution_2lines)

# final figure # ----
graph<-ggdraw(xlim = c(0, 2.15), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(p2, 1, 0, 1, 1)+
  draw_plot(legend, 2.05, 0.25, 0.05, 0.5)+
  draw_plot_label(c("A","B"), c(0,1), c(1,1), size = 30)
ggsave(paste(path_figure,"supp_eigen_disp1.pdf",sep=""), graph, width = 18, height = 9, device = cairo_pdf)

## time series # ----
nSpecies=2
nCommunity=2
gamma=3 # asymmetry coefficient
d=1e6

params_data<-expand.grid(pert=list(list(c(2,1))), # perturbation of predators in patch 1
                         disp=list(c(1,0)), # dispersal of prey
                         gamma=gamma,
                         d=d,
                         model="pert_21")
params_data<-rbind(params_data,expand.grid(pert=list(list(c(2,2))), # perturbation of predators in patch 2
                                           disp=list(c(1,0)), # dispersal of prey
                                           gamma=gamma,
                                           d=d,
                                           model="pert_22"))
params_data<-merge(params_data_original,params_data)
params_data$omega=1#params_data$gamma
# sets asymmetry
for(i in 1:dim(params_data)[1]){
  params_data$asym[i]=list(list(rep(params_data$gamma[i],nSpecies), # patch with modified attack rate
                                rep(1,nSpecies))) # reference patch
}
# sets dispersal
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}

B0<-matrix(0,dim(params_data)[1],nSpecies*nCommunity)
colnames(B0)<-c("B_11","B_21","B_12","B_22")
for(i in 1:dim(params_data)[1]){
  params<-get_params(params_data,i)
  B0[i,]<-equilibrium_disp_prey(params,nSpecies,nCommunity)
}

# integration time and pulse perturbations depending on parameters
t_max=10
t_step=0.1
pert_factor=1.2

# simulations
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1],.packages=c("deSolve")) %dopar% time_series_perturbation(params_data[i,], B0[i,], t_max, t_step, pert_factor,nSpecies, nCommunity)
stopCluster(cl)

TS<-create_TS(params_data, results, B0)
TS<-TS[,which(names(TS)%in%c("time","model","gamma","ea","ma",colnames(B0)))]
TS<-time_series_for_plot(TS, 4)

levels(TS$model)<-c(pert_21,pert_22)

p1<-ggplot(data=TS[TS$model==pert_21,])+
  geom_line(aes(time,biomass,colour=species,linetype=community),size=1.5)+
  corr_colour_TL_2+
  patch_line+
  theme
legend<-get_legend(p1)

p1<-ggplot(data=TS[TS$model==pert_21,])+
      geom_line(aes(time,biomass,colour=species,linetype=community),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed, scales="free")+
      corr_colour_TL_2+
      patch_line+
      theme+theme(legend.position="none")+
      scale_x_continuous(breaks=seq(0,t_max,2))+
      #scale_y_continuous(breaks=seq(0,pert_factor,0.1))+
      xlab("Time")+
      ylab("Scaled biomass")+
      ggtitle("predators perturbed in patch #1")

p2<-ggplot(data=TS[TS$model==pert_22,])+
      geom_line(aes(time,biomass,colour=species,linetype=community),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed, scales="free")+
      corr_colour_TL_2+
      patch_line+
      theme+theme(legend.position="none")+
      scale_x_continuous(breaks=seq(0,t_max,2))+
      #scale_y_continuous(breaks=seq(0,pert_factor,0.1))+
      xlab("Time")+
      ylab("Scaled biomass")+
      ggtitle("predators perturbed in patch #2")

graph<-ggdraw(xlim = c(0, 2.15), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(p2, 1, 0, 1, 1)+
  draw_plot(legend, 2.05, 0.25, 0.05, 0.5)+
  draw_plot_label(c("A","B"), c(0,1), c(1,1), size = 30)
ggsave(paste(path_figure,"supp_TS_disp1.pdf",sep=""), graph, width = 14, height = 7, device = cairo_pdf)

#### ENVIRONMENTAL PERTURBATION #### ----
# parameters and simulations - gamma - omega # ----
nSpecies=2
nCommunity=2
gamma=seq(1,10,0.1) # asymmetry coefficient
#omega=seq(1,10,0.1) # asymmetry in growth rate
d=1e6 # high dispersal

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
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)
params_data$z=1 # environmental perturbations

B<-equilibrium_analytical_TL2_disp_pred(params_data) # biomass at equilibrium
coexistence<-coexistence_TL2_disp_pred(params_data) # coexistence of all species

# index of the parameters leading to significant coexistence
index<-which(coexistence$coexistence_significative==1)
params_data<-params_data[index,]
B<-B[index,] # selects the biomasses allowing coexistence

# simulations
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1],.packages=NULL) %dopar% analysis(params_data[i,], B[i,], nSpecies, nCommunity)
stopCluster(cl)
list_data<-create_data(params_data, results, nSpecies, nCommunity)

## correlation inter-patch # ----
data_C<-list_data$data_C
databis<-data_C[data_C$ea==paste(ea1) & data_C$ma==paste(ma1),
                which(names(data_C)%in%c("ea","ma","d","gamma","model","C_11_12",'C_21_22'))]
h_value=databis$C_11_12[databis$gamma==1][1]
databis<-table_for_plot(databis,5,"correlation")

p1<-ggplot(data=databis)+
  geom_line(aes(gamma,correlation,colour=species,linetype=model),size=1.5)+
  geom_hline(yintercept=h_value,linetype="dashed")+
  #facet_grid(ma~ea, labeller=label_parsed)+
  perturbation_prey_line+
  corr_colour_TL_2+
  theme+
  x_axis_gamma+
  ylim(-1,1)+
  xlab(label_gamma)+
  ylab(label_correlation)

## biomass CV at population scale # ----
data_CV_pop<-list_data$data_CV_pop
data_CV_pop<-data_CV_pop[data_CV_pop$ea==paste(ea1) & data_CV_pop$ma==paste(ma1),
                         which(names(data_CV_pop)%in%c("ea","ma","d","gamma","omega","model",list_data$CV_pop_names))]
data_CV_pop<-table_for_plot(data_CV_pop,6,"CV")
levels(data_CV_pop$model)<-c(pert_11,pert_12)

p2<-ggplot(data=data_CV_pop)+
  geom_line(aes(gamma,CV,colour=species,linetype=community),size=1.5)+
  facet_wrap(~model)+
  corr_colour_TL_2+
  patch_line+
  theme +
  x_axis_gamma+
  y_axis_log10+
  xlab(label_gamma)+
  ylab(label_CV_pop)

## biomass CV at metapopulation scale # ----
data_CV_species<-list_data$data_CV_species
data_CV_species<-data_CV_species[data_CV_species$ea==paste(ea1) & data_CV_species$ma==paste(ma1),
                                 which(names(data_CV_species)%in%c("ea","ma","d","gamma","omega","model",list_data$CV_species_names))]
data_CV_species<-table_for_plot(data_CV_species,6,"CV")
data_CV_species$species=data_CV_species$community
data_CV_species$community=NULL
levels(data_CV_species$model)<-c(pert_11,pert_12)

p3<-ggplot(data=data_CV_species)+
  geom_line(aes(gamma,CV,colour=species),size=1.5)+
  facet_wrap(~model)+
  corr_colour_TL_2+
  theme +
  x_axis_gamma+
  y_axis_log10+
  xlab(label_gamma)+
  ylab(label_CV_metapop)

## biomass CV at metacommunity scale # ----
data_CV<-list_data$data_CV
levels(data_CV$model)<-c(pert_11,pert_12)

p4<-ggplot(data=data_CV[data_CV$ea==paste(ea1) & data_CV$ma==paste(ma1),])+
  geom_line(aes(gamma,CV_tot),size=1.5)+
  facet_wrap(~model)+
  corr_colour_TL_2+
  theme +
  x_axis_gamma+
  y_axis_log10+
  xlab(label_gamma)+
  ylab(label_CV_metacom)

# final figure # ----
graph<-ggdraw(xlim = c(0, 2.08), ylim = c(0, 2)) +
  draw_plot(p1, 0.08, 1, 1, 1)+
  draw_plot(p2, 1.08, 1, 1, 1)+
  draw_plot(p3, 0.08, 0, 1, 1)+
  draw_plot(p4, 1.08, 0, 1, 1)+
  draw_plot_label(c("A","B","C","D"), c(0,1,0,1), c(2,2,1,1), size = 30)
ggsave(paste(path_figure,"supp_env_pert_11.pdf",sep=""), graph, width = 16, height = 10, device = cairo_pdf)

#### CORRELATED ENVIRONMENTAL PERTURBATION #### ----
# parameters and simulations - gamma - omega # ----
nSpecies=2
nCommunity=2
gamma=seq(1,10,0.1) # asymmetry coefficient
#omega=seq(1,10,0.1) # asymmetry in growth rate
d=1e6 # high dispersal

params_data<-expand.grid(pert=list(list(c(1,1),c(1,2))), # perturbation of prey in patch 1 and 2
                         disp=list(c(0,1)), # dispersal of predators
                         gamma=gamma,
                         d=d,
                         corr=c(0,0.5,1), # correlation of perturbations
                         model="pert_11_12")
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
# sets perturbation covariance
VE<-matrix(0,nSpecies*nCommunity,nSpecies*nCommunity)
for (i in 1:nSpecies){
  for (j in 1:nCommunity){
    VE[i+nSpecies*(j-1),seq(i,i+nSpecies*(nCommunity-1),nSpecies)]=sigma^2 # correlated perturbations
  }
}
params_data$VE<-list(VE)
for(i in 1:dim(params_data)[1]){
  params_data$VE[[i]]<-params_data$VE[[i]]*params_data$corr[i] # sets the 
  diag(params_data$VE[[i]])<-rep(sigma^2,nSpecies*nCommunity)
}
params_data$z=1 # environmental perturbations

B<-equilibrium_analytical_TL2_disp_pred(params_data) # biomass at equilibrium
coexistence<-coexistence_TL2_disp_pred(params_data) # coexistence of all species

# index of the parameters leading to significant coexistence
index<-which(coexistence$coexistence_significative==1)
params_data<-params_data[index,]
B<-B[index,] # selects the biomasses allowing coexistence

# simulations
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1],.packages=NULL) %dopar% analysis(params_data[i,], B[i,], nSpecies, nCommunity)
stopCluster(cl)
list_data<-create_data(params_data, results, nSpecies, nCommunity)

## variance ratio # ----
data_V<-list_data$data_V
databis<-data_V[data_V$ea==paste(ea1) & data_V$ma==paste(ma1),
                which(names(data_V)%in%c("ea","ma","d","gamma","corr","V_11_11",'V_12_12'))]
databis<-table_for_plot(databis,5,"variance")
databis$variance<-databis$variance/sigma
databis$corr<-as.factor(databis$corr)
databis$community<-as.factor(databis$community)
levels(databis$community)<-c("patch #1 (fast)","patch #2 (slow)")

p1<-ggplot(data=databis)+
  geom_line(aes(gamma,variance,alpha=corr),colour="dodgerblue3",size=1.5)+
  facet_wrap(~community)+
  corr_alpha+
  theme+
  x_axis_gamma+
  y_axis_log10+
  xlab(label_gamma)+
  ylab("Ratio of prey biomass variance\nto perturbation variance")
ggsave(paste(path_figure,"supp_env_pert_variance_ratio.pdf",sep=""),p1, width = 8, height = 4, device = cairo_pdf)

## correlation inter-patch # ----
data_C<-list_data$data_C
databis<-data_C[data_C$ea==paste(ea1) & data_C$ma==paste(ma1),
                which(names(data_C)%in%c("ea","ma","d","gamma","corr","C_11_12",'C_21_22'))]
h_value=databis$C_11_12[databis$gamma==1][1]
databis<-table_for_plot(databis,5,"correlation")
databis$corr<-as.factor(databis$corr)

p1<-ggplot(data=databis)+
  geom_line(aes(gamma,correlation,colour=species,alpha=corr),size=1.5)+
  corr_alpha+
  corr_colour_TL_2+
  theme
legend<-get_legend(p1)

p1<-ggplot(data=databis)+
      geom_line(aes(gamma,correlation,colour=species,alpha=corr),size=1.5)+
      #geom_hline(yintercept=h_value,linetype="dashed")+
      #facet_grid(ma~ea, labeller=label_parsed)+
      corr_alpha+
      corr_colour_TL_2+
      theme+theme(legend.position="none")+
      x_axis_gamma+
      ylim(-1,1)+
      xlab(label_gamma)+
      ylab(label_correlation)

## biomass CV at population scale # ----
data_CV_pop<-list_data$data_CV_pop
data_CV_pop<-data_CV_pop[data_CV_pop$ea==paste(ea1) & data_CV_pop$ma==paste(ma1),
                         which(names(data_CV_pop)%in%c("ea","ma","d","gamma","omega","corr",list_data$CV_pop_names))]
data_CV_pop<-table_for_plot(data_CV_pop,6,"CV")
data_CV_pop$corr<-as.factor(data_CV_pop$corr)
data_CV_pop$community<-as.factor(data_CV_pop$community)
levels(data_CV_pop$community)<-c("patch #1 (fast)","patch #2 (slow)")

p2<-ggplot(data=data_CV_pop)+
      geom_line(aes(gamma,CV,colour=species,alpha=corr),size=1.5)+
      facet_wrap(~community)+
      corr_colour_TL_2+
      corr_alpha+
      theme+theme(legend.position="none")+
      x_axis_gamma+
      y_axis_log10+
      xlab(label_gamma)+
      ylab(label_CV_pop)

## biomass CV at metapopulation scale # ----
data_CV_species<-list_data$data_CV_species
data_CV_species<-data_CV_species[data_CV_species$ea==paste(ea1) & data_CV_species$ma==paste(ma1),
                                 which(names(data_CV_species)%in%c("ea","ma","d","gamma","omega","corr",list_data$CV_species_names))]
data_CV_species<-table_for_plot(data_CV_species,6,"CV")
data_CV_species$species=data_CV_species$community
data_CV_species$community=NULL
data_CV_species$corr<-as.factor(data_CV_species$corr)

p3<-ggplot(data=data_CV_species)+
      geom_line(aes(gamma,CV,colour=species,alpha=corr),size=1.5)+
      corr_colour_TL_2+
      corr_alpha+
      theme+theme(legend.position="none")+
      x_axis_gamma+
      y_axis_log10+
      xlab(label_gamma)+
      ylab(label_CV_metapop)

## biomass CV at metacommunity scale # ----
data_CV<-list_data$data_CV
data_CV$corr<-as.factor(data_CV$corr)

p4<-ggplot(data=data_CV[data_CV$ea==paste(ea1) & data_CV$ma==paste(ma1),])+
      geom_line(aes(gamma,CV_tot,alpha=corr),size=1.5)+
      corr_alpha+
      theme+theme(legend.position="none")+
      x_axis_gamma+
      y_axis_log10+
      xlab(label_gamma)+
      ylab(label_CV_metacom)

# final figure # ----
graph<-ggdraw(xlim = c(0, 2.5), ylim = c(0, 2)) +
  draw_plot(p1, 0.08, 1, 1, 1)+
  draw_plot(p2, 1.16, 1, 1, 1)+
  draw_plot(p3, 0.08, 0, 1, 1)+
  draw_plot(p4, 1.16, 0, 1, 1)+
  draw_plot(legend, 2.3, 0.75, 0.05, 0.5)+
  draw_plot_label(c("A","B","C","D"), c(0,1.08,0,1.08), c(2,2,1,1), size = 30)
ggsave(paste(path_figure,"supp_env_pert.pdf",sep=""), graph, width = 14, height = 10, device = cairo_pdf)

## time series # ----
nSpecies=2
nCommunity=2
gamma=3 # asymmetry coefficient
d=1e6

params_data<-expand.grid(pert=list(list(c(1,1),c(1,2))), # perturbation of prey in patch 1
                         disp=list(c(0,1)), # dispersal of predators
                         gamma=c(1,3),
                         d=d,
                         model="pert_11_12")
params_data<-merge(params_data_original,params_data)
params_data$omega=params_data$gamma
params_data<-params_data[params_data$ea==paste(ea1) & params_data$ma==paste(ma1),]
# sets asymmetry
for(i in 1:dim(params_data)[1]){
  params_data$asym[i]=list(list(rep(params_data$gamma[i],nSpecies), # patch with modified attack rate
                                rep(1,nSpecies))) # reference patch
}
# sets dispersal
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}

# biomasses computed analytically (for d>>1, B_21=B_22)
B0<-equilibrium_analytical_TL2_disp_pred(params_data)

# integration time and pulse perturbations depending on parameters
t_max=10
t_step=0.1
pert_factor=1.2

# simulations
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1],.packages=c("deSolve")) %dopar% time_series_perturbation(params_data[i,], B0[i,], t_max, t_step, pert_factor,nSpecies, nCommunity)
stopCluster(cl)

TS<-create_TS(params_data, results, B0)
TS<-TS[,which(names(TS)%in%c("time","model","gamma",colnames(B0)))]
TS<-time_series_for_plot(TS, 2)

p1<-ggplot(data=TS)+
  geom_line(aes(time,biomass,colour=species,linetype=community),size=1.5)+
  facet_wrap(~gamma)+
  corr_colour_TL_2+
  patch_line+
  theme+
  scale_x_continuous(breaks=seq(0,t_max,2))+
  scale_y_continuous(breaks=seq(0,pert_factor,0.1))+
  xlab("Time")+
  ylab("Scaled biomass")
ggsave(paste(path_figure,"supp_env_corr_TS.pdf",sep=""), p1, width = 11, height = 5, device = cairo_pdf)
