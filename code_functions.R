### Jacobian matrix ####
# Jacobian single patch J(k)
jacobian_single<-function(B,params,nSpecies,CommunityID){
  with(params,{
    asym<-params$asym[[CommunityID]]
    J<-matrix(0,nrow=nSpecies,ncol=nSpecies)
    J[1,1]=D*(g*asym[1]/D-2*B[1])
    if(nSpecies>1){
      J[1,1]=J[1,1]-D*m*a*asym[2]*B[2] # basal species
      J[nSpecies,nSpecies]=D*m^(nSpecies-1)*(-r/D-2*B[nSpecies]+e*a*asym[nSpecies]*B[nSpecies-1]) # top species
      for(i in 1:(nSpecies-1)){
        J[i,i+1]=-D*m^(i-1)*m*a*asym[i+1]*B[i] # effect of predators
      }
      for(i in 2:nSpecies){
        J[i,i-1]=D*m^(i-1)*e*a*asym[i]*B[i] # effect of prey
      }
    }
    if(nSpecies>2){
      for(i in 2:(nSpecies-1)){
        J[i,i]=D*m^(i-1)*(-r/D-2*B[i]+e*a*asym[i]*B[i-1]-m*a*asym[i+1]*B[i+1]) # intermediate species
      }
    }
    return(J)
  })
}

# Jacobian of the intrapatch dynamics J'
jacobian_intra_patch<-function(B,params,nSpecies,nCommunity){
  Z<-matrix(0,nrow=nSpecies,ncol=nSpecies)
  M<-NULL
  L<-NULL
  J<-NULL
  for(i in 1:nCommunity){
    L<-NULL
    J<-jacobian_single(B[nSpecies*(i-1)+c(1:nSpecies)],params,nSpecies,i)
    for(j in 1:nCommunity){
      if(j==i){
        L<-cbind(L,J)
      }
      else{
        L<-cbind(L,Z)
      }
    }
    M<-rbind(M,L)
  }
  return(M)
}

# Jacobian of the dispersal dynamics P' (linear dispersal components)
jacobian_disp<-function(B,params,nSpecies,nCommunity){
  with(params,{
    dim=nSpecies*nCommunity
    P<-matrix(0,nrow=dim,ncol=dim)
    # effects of species i on itself
    for(i in 1:nSpecies){
      for(j in 1:nCommunity){
        for(k in 1:nCommunity){
          row=(j-1)*nSpecies+i
          col=(k-1)*nSpecies+i
          if(j==k){
            P[row,col]=-disp[i]
          }
          else{
            P[row,col]=disp[i]/(nCommunity-1)
          }
        }
      }
    }
    # Finalising the matrix
    for(i in 1:nSpecies){
      for(j in 1:nCommunity){
        P[(j-1)*nSpecies+i,]=D*m^(i-1)*P[(j-1)*nSpecies+i,]
      }
    }
    return(P)
  })
}

# Final Jacobian for the multiRoot algorithm
jacobian_multiroot<-function(B,params){
  J<-jacobian_intra_patch(B,params,params$nSpecies,params$nCommunity)
  P<-jacobian_disp(B,params,params$nSpecies,params$nCommunity)
  return(J+P)
}

### Biomass at equilibrium ####
# analytical calculation of biomass - TL2 predator dispersal
equilibrium_analytical_TL2_disp_pred<-function(params_data){
  params_data$lambda<-params_data$e*params_data$a^2*params_data$m
  B0<-matrix(0,dim(params_data)[1],4)
  colnames(B0)<-c("B_11","B_21","B_12","B_22")
  # biomasses computed analytically (for d>>1, B_21=B_22)
  B0[,1]=params_data$g/params_data$D*
    (2*params_data$omega + params_data$lambda*params_data$omega - params_data$lambda*params_data$gamma)/
    (2+params_data$lambda*(params_data$gamma^2+1))
  B0[,2]=2*params_data$e*params_data$a*params_data$g/params_data$D*
    (1+params_data$omega*params_data$gamma)/
    (2+params_data$lambda*(params_data$gamma^2+1))/2
  B0[,3]=params_data$g/params_data$D*
    (2+params_data$lambda*params_data$gamma*(params_data$gamma-params_data$omega))/
       (2+params_data$lambda*(params_data$gamma^2+1))
  B0[,4]=B0[,2]
  B0[,1][B0[,1]<1e-6]=NA
  B0[,3][B0[,3]<1e-6]=NA
  return(B0)
}

# species coexistence - TL2 predator dispersal
coexistence_TL2_disp_pred<-function(params_data){
  data<-params_data
  data$lambda<-data$e*data$a^2*data$m
  # criterion on prey in patch #1
  data$crit_patch_1=data$omega*(2+data$lambda)/data$lambda
  # criteria on prey in patch #2
  data$crit_delta=sqrt(8/data$lambda)
  data$crit_gamma_1=(data$lambda*data$omega + sqrt(data$lambda*(data$lambda*data$omega^2-8)))/2/data$lambda
  data$crit_gamma_2=(data$lambda*data$omega - sqrt(data$lambda*(data$lambda*data$omega^2-8)))/2/data$lambda
  # initiate the column with no coexistence
  data$coexistence=0
  
  for(i in 1:nrow(data)){
    if(data$gamma[i]<data$crit_patch_1[i] & # prey biomass in patch #1 positive
       (data$omega[i]<data$crit_delta[i] # prey biomass in patch #2 is always positive
        | data$omega[i]>data$crit_delta[i] & # two roots for gamma (negative biomass between the roots)
        (data$gamma[i]>data$crit_gamma_1[i] | data$gamma[i]<data$crit_gamma_2[i]))){
      data$coexistence[i]=1
    }
  }
  data$coexistence<-as.factor(data$coexistence)
  levels(data$coexistence)[levels(data$coexistence)=="1"]="coexistence"
  levels(data$coexistence)[levels(data$coexistence)=="0"]="no coexistence"
  
  # biomass
  B<-as.data.frame(equilibrium_analytical_TL2_disp_pred(params_data))
  
  # index of the parameters leading to significant coexistence
  data$coexistence_significative=0
  data$coexistence_significative[which(data$coexistence=="coexistence" & is.na(B$B_11)==F & is.na(B$B_12)==F)]=1
  return(data)
}

# biomass
equilibrium_ode<-function(params_data_row, B0, t_max, t_step, nSpecies, nCommunity){
  params<-set_params(params_data_row)
  time<-seq(0,t_max,t_step)
  TS<-time_series(params_data_row, B0, t_max, t_step, nSpecies, nCommunity)
  B<-as.numeric(TS[nrow(TS),2:ncol(TS)])
  return(B)
}

# Biomasses at equilibrium in the symmetric case
equilibrium_symmetric<-function(params,nSpecies,nCommunity){ # compute the biomasses at equilibrium
  with(params,{
    A<-diag(rep(-1,nSpecies*nCommunity))
    for(j in 1:nCommunity){
      for(i in 2:nSpecies){
        A[(j-1)*nSpecies+i,(j-1)*nSpecies+i-1]=e*a
      }
      for(i in 1:(nSpecies-1)){
        A[(j-1)*nSpecies+i,(j-1)*nSpecies+i+1]=-m*a
      }
    }
    B<-matrix(r/D,
              nrow = nSpecies*nCommunity,
              ncol = 1,
              byrow = TRUE)
    for(j in 1:nCommunity){
      B[(j-1)*nSpecies+1,1]=-g/D
    }
    C<-matrix(0,
              nrow = nSpecies*nCommunity,
              ncol = 1,
              byrow = TRUE)
    C<-solve(A) %*% B
    return(as.numeric(C))
  })
}

# Biomasses at equilibrium in an isolated community (in the case of asymmetry)
equilibrium_isolated<-function(params,nSpecies,communityID){ # compute the biomasses at equilibrium
  with(params,{
    A<-diag(rep(-1,nSpecies))
    for(j in 1:nCommunity){
      for(i in 2:nSpecies){
        A[i,i-1]=e*a*asym[[communityID]][i]
      }
      for(i in 1:(nSpecies-1)){
        A[i,i+1]=-m*a*asym[[communityID]][i+1]
      }
    }
    B<-matrix(r/D,
              nrow = nSpecies,
              ncol = 1,
              byrow = TRUE)
    for(j in 1:nCommunity){
      B[1,1]=-g*asym[[communityID]][1]/D
    }
    C<-matrix(0,
              nrow = nSpecies,
              ncol = 1,
              byrow = TRUE)
    C<-solve(A) %*% B
    return(as.numeric(C))
  })
}

# Biomasses at equilibrium when top predators disperse a lot
equilibrium_disp_pred<-function(params,nSpecies,nCommunity){ # compute the biomasses at equilibrium
  with(params,{
    # e=params$e
    # a=params$a
    # m=params$m
    # asym=params$asym
    # r=params$r
    # g=params$g
    # D=params$D
    nSpecies=nSpecies-1 # top predators are added at the end
    A<-diag(rep(-1,nSpecies*nCommunity+1)) # +1 is the line of top predators (last line of the matrix)
    for(j in 1:nCommunity){
      if(nSpecies>=2){
        for(i in 2:nSpecies){
          A[(j-1)*nSpecies+i,(j-1)*nSpecies+i-1]=e*a*asym[[j]][i]
        }
      }
      for(i in 1:nSpecies){
        A[(j-1)*nSpecies+i,(j-1)*nSpecies+i+1]=-m*a*asym[[j]][i+1]
        if(i==nSpecies){
          A[(j-1)*nSpecies+i,(j-1)*nSpecies+i+1]=0
          A[(j-1)*nSpecies+i,nSpecies*nCommunity+1]=-m*a*asym[[j]][i+1]/nCommunity # because we consider the total biomass of top predators
        }
      }
    }
    for(j in 1:nCommunity){
      A[dim(A)[1],nSpecies+(j-1)*nSpecies]=e*a*asym[[j]][nSpecies+1]
    }
    B<-matrix(r/D,
              nrow = nSpecies*nCommunity+1,
              ncol = 1,
              byrow = TRUE)
    for(j in 1:nCommunity){
      B[(j-1)*nSpecies+1,1]=-g*asym[[j]][1]/D
    }
    B[nrow(B),1]=2*r/D
    C<-matrix(0,
              nrow = nSpecies*nCommunity+1,
              ncol = 1,
              byrow = TRUE)
    C<-solve(A) %*% B
    C<-as.numeric(C)
    B<-rep(C[length(C)]/nCommunity,(nSpecies+1)*nCommunity)
    for(j in 1:nCommunity){
      B[c(1:nSpecies)+(nSpecies+1)*(j-1)]=C[c(1:nSpecies)+nSpecies*(j-1)]
    }
    return(B)
  })
}

# Biomasses at equilibrium when basal species disperse a lot
equilibrium_disp_prey<-function(params,nSpecies,nCommunity){ # compute the biomasses at equilibrium
  with(params,{
    # e=params$e
    # a=params$a
    # m=params$m
    # asym=params$asym
    # r=params$r
    # g=params$g
    # D=params$D
    nSpecies=nSpecies-1 # top predators are added at the end
    A<-diag(rep(-1,nSpecies*nCommunity+1)) # +1 is the line of basal species (first line of the matrix)
    for(j in 1:nCommunity){
      if(nSpecies>=2){
        for(i in 1:(nSpecies-1)){
          A[(j-1)*nSpecies+i+1,(j-1)*nSpecies+i+2]=-m*a*asym[[j]][i+2] # +1 because of basal species (and +1 because it is the asymmetry coeff of predators i+1)
        }
      }
      if(nSpecies>=1){
        for(i in 1:nSpecies){
          A[(j-1)*nSpecies+i+1,(j-1)*nSpecies+i]=e*a*asym[[j]][i+1]
          if(i==1){
            A[(j-1)*nSpecies+i+1,(j-1)*nSpecies+i]=0
            A[(j-1)*nSpecies+i+1,1]=e*a*asym[[j]][i+1]/nCommunity # because we consider the total biomass of basal species
          }
        }
      }
    }
    for(j in 1:nCommunity){
      A[1,1+(j-1)*nSpecies+1]=-m*a*asym[[j]][2]
    }
    B<-matrix(r/D,
              nrow = nSpecies*nCommunity+1,
              ncol = 1,
              byrow = TRUE)
    B[1,1]=0
    for(j in 1:nCommunity){
      B[1,1]=B[1,1]-g*asym[[j]][1]/D
    }
    C<-matrix(0,
              nrow = nSpecies*nCommunity+1,
              ncol = 1,
              byrow = TRUE)
    C<-solve(A) %*% B
    C<-as.numeric(C)
    B<-rep(C[1]/nCommunity,(nSpecies+1)*nCommunity)
    for(j in 1:nCommunity){
      B[c(2:(nSpecies+1))+(nSpecies+1)*(j-1)]=C[c(1:nSpecies)+nSpecies*(j-1)+1]
    }
    return(B)
  })
}

# Biomasses at equilibrium in a tri-trophic food chain in which the intermediate species disperse a lot
equilibrium_TL3_disp_2<-function(params){ # compute the biomasses at equilibrium
  with(params,{
    e=params$e
    a=params$a
    m=params$m
    asym=params$asym
    r=params$r
    g=params$g
    D=params$D
    A<-diag(rep(-1,5))
    A[1,5]=-asym[[1]][2]*m*a/2
    A[2,5]=asym[[1]][3]*e*a/2
    A[3,5]=-m*a/2
    A[4,5]=e*a/2
    A[5,1]=asym[[1]][2]*e*a/2
    A[5,2]=-asym[[1]][3]*m*a/2
    A[5,3]=e*a/2
    A[5,4]=-m*a/2
    B<-matrix(r/D,
              nrow = 5,
              ncol = 1,
              byrow = TRUE)
    B[1,1]=-asym[[1]][1]*g/D
    B[3,1]=-g/D
    B[5,1]=2*B[5,1]
    C<-matrix(0,
              nrow = 5,
              ncol = 1,
              byrow = TRUE)
    C<-solve(A) %*% B
    B<-c(C[1,1],C[5,1]/2,C[2,1],C[3,1],C[5,1]/2,C[4,1])
    return(B)
  })
}

### Stability analysis ####
# T matrix
T_matrix<-function(pert, B, z, nSpecies, nCommunity){
  T<-matrix(0,nSpecies*nCommunity,nSpecies*nCommunity)
  coord=NULL
  for(i in 1:length(pert)){
    coord=(pert[[i]][2]-1)*nSpecies+pert[[i]][1] # pert->list of vector containing the trophic level and the patch of the perturbed species (e.g. (1,2) species 1 in patch 2)
    T[coord,coord]=1
  }
  T<-T*diag(B)^z
  return(T)
}

# Lyapunov equation
lyapunov<-function(J, T, VE, nSpecies, nCommunity){
  TVT<-T%*%VE%*%t(T)
  TVT<-matrix(array(TVT),ncol=1)
  kron<-kronecker(J,diag(rep(1,nSpecies*nCommunity))) + kronecker(diag(rep(1,nSpecies*nCommunity)),J)
  return(-solve(kron)%*%TVT)
}

# importance of dispersal relative to local demographic processes
dispersal_importance<-function(B, params, nSpecies, nCommunity){
  with(params,{
    # effects of trophic interactions
    trophic<-rep(0,nSpecies*nCommunity)
    ID=0
    for(j in 1:nCommunity){
      ID=(j-1)*nSpecies
      trophic[ID+1]=B[ID+1]*(g/D + B[ID+1])
      if(nSpecies>1){
        for(i in 2:nSpecies){
          trophic[ID+i]=B[ID+i]*(r/D + B[ID+i] + e*a*asym[[j]][i]*B[ID+i-1])
        }
        for(i in 1:(nSpecies-1)){
          trophic[ID+i]=trophic[ID+i] + B[ID+i]*m*a*asym[[j]][i+1]*B[ID+i+1]
        }
      }
    }
    dispersal<-rep(0,nSpecies*nCommunity)
    d_matrix<-matrix(0,nSpecies,nCommunity) # matrix containing the dispersal terms
    B_matrix<-matrix(B,nSpecies,nCommunity) # matrix containg the biomasses
    # final equations
    for(i in 1:nSpecies){
      d_matrix[i,]=disp[i]*B_matrix[i,]
    }
    for(i in 1:nSpecies){
      for(j in 1:nCommunity){
        dispersal[(j-1)*nSpecies+i]=d_matrix[i,j]+(sum(d_matrix[i,])-d_matrix[i,j])/(nCommunity-1)
      }
    }
    return(dispersal/(trophic+dispersal))
  })
}

# CV of regional populations of species
CV_species<-function(V,B,nSpecies,nCommunity){
  V<-matrix(V,nSpecies*nCommunity)
  B_species<-rep(0,nSpecies)
  CV_species<-rep(0,nSpecies)
  for(i in 1:nSpecies){
    index=seq(i,i+nSpecies*(nCommunity-1),nSpecies) # index of species i in vectors
    B_species[i]=sum(B[index]) # total biomass of species i
    index<-expand.grid(index,index) # index of variance and covariances of populations of species i in matrix V
    CV_species[i]=sum(V[index$Var1,index$Var2]) # variance of the total biomass of species i
  }
  CV_species<-sqrt(CV_species)/B_species # CV
  return(CV_species)
}

# Parallelised analytic calculation 
analysis<-function(params_data_row, B, nSpecies, nCommunity){
  params<-set_params(params_data_row)
  J<-jacobian_intra_patch(B,params,nSpecies,nCommunity)
  P<-jacobian_disp(B,params,nSpecies,nCommunity)
  T<-T_matrix(params$pert,B,params$z,nSpecies,nCommunity)
  V<-lyapunov(J+P,T,params$VE,nSpecies,nCommunity)
  C<-cov2cor(matrix(as.numeric(V),nSpecies*nCommunity,nSpecies*nCommunity))
  M1<-dispersal_importance(B,params,nSpecies,nCommunity)
  CV_pop<-sqrt(diag(matrix(V,nSpecies*nCommunity)))/B # biomass CV of each population
  CV_spe<-CV_species(V,B,nSpecies,nCommunity) # biomass CV of the regional population of each species
  CV<-sqrt(c(sum(diag(matrix(V,nSpecies*nCommunity))),sum(V)))/sum(B) # average CV of populations and CV of total biomass
  Press<-0#-solve(diag(B))%*%solve(J+P)%*%diag(B)
  eigen<-eigen(J+P) # eigenvalues and eigenvectors
  resilience<-max(Re(eigen$values)) # lead eigen value
  lead<-which(Re(eigen$values)==resilience) # index of the lead eigenvalue
  E<-abs(Re(eigen$vectors[,lead])) # lead eigenvector
  E<-E/sum(E) # relative contribution of each element
  jacobian<-J+P
  return(list(B=B,
              V=as.numeric(V),
              C=as.numeric(C),
              M1=M1,
              CV_pop=CV_pop,
              CV_species=CV_spe,
              CV=CV,
              Press=as.numeric(Press),
              resilience=-resilience,
              E=E,
              jacobian=as.numeric(jacobian)))
}

### ODE ####
# ODE of the system for the numerical simulation of time series
ODE<-function(B, params){
  with(params,{
    dB<-rep(0,nSpecies*nCommunity)
    # intra-patch dynamics
    ID=0
    for(j in 1:nCommunity){
      ID=(j-1)*nSpecies
      dB[ID+1]=dB[ID+1]+D*B[ID+1]*(g*asym[[j]][1]/D - B[ID+1])
      if(nSpecies>1){
        for(i in 2:nSpecies){
          dB[ID+i]=dB[ID+i]+D*m^(i-1)*B[ID+i]*(-r/D - B[ID+i] + e*a*asym[[j]][i]*B[ID+i-1])
        }
        for(i in 1:(nSpecies-1)){
          dB[ID+i]=dB[ID+i]+D*m^(i-1)*B[ID+i]*(- m*a*asym[[j]][i+1]*B[ID+i+1])
        }
      }
    }
    # dispersal dynamics
    d_matrix<-matrix(0,nSpecies,nCommunity) # matrix containing the dispersal terms
    B_matrix<-matrix(B,nSpecies,nCommunity) # matrix containing the biomasses
    # final equations
    for(i in 1:nSpecies){
      d_matrix[i,]=D*m^(i-1)*disp[i]*B_matrix[i,]
    }
    for(i in 1:nSpecies){
      for(j in 1:nCommunity){
        ID=(j-1)*nSpecies
        dB[ID+i]=dB[ID+i]-(1+1/(nCommunity-1))*d_matrix[i,j]+sum(d_matrix[i,])/(nCommunity-1)
      }
    }
    return(dB)
  })
}

# function for the ODE solver
ODE_function<-function(t, B, params){
  dB<-ODE(B, params)
  dB[B<=0]=0 # extinction
  return(list(dB))
}

# root function evaluating if a species has a negative biomass
root_function<-function(t, B, params) {
  B_root=B-1e-8 # detection of negative values
  return(B_root)
}

# event function setting negative biomass equal to zero
event_function<-function(t, B, params) {
  B[which(B<1e-8)]=0 # set negative biomass equal to zero
  return(B)
}

# time series
time_series<-function(params_data_row, B0, t_max, t_step, nSpecies, nCommunity){
  params<-set_params(params_data_row)
  time<-seq(0, t_max, t_step)
  #TS<-as.data.frame(rk(B0, time, ODE_function, params, method="rk45dp7", hini=0.001))
  TS<-as.data.frame(lsoda(B0, time, ODE_function, params,
                          hini=0.001, # initial step size
                          rootfun=root_function))#, events = list(func=event_function, root = TRUE))) # solver with adaptive time step
  if(length(time)>1000){
    TS<-TS[seq(1,length(time),floor(length(time)/1000)),] # keeps only 1000 points
  }
  return(TS)
}

# time series after a pulse perturbation
time_series_perturbation<-function(params_data_row, B0, t_max, t_step, pert_factor, nSpecies, nCommunity){
  params<-set_params(params_data_row)
  time<-seq(0, t_max, t_step)
  # perturbation of predators in patch #1
  B<-B0
  coord=NULL
  pert<-params$pert# selects the list of perturbed species
  for(i in 1:length(pert)){
    coord=(pert[[i]][2]-1)*nSpecies+pert[[i]][1] # pert->list of vector containing the trophic level and the patch of the perturbed species (e.g. (1,2) species 1 in patch 2)
    B[coord]=B[coord]*pert_factor # generate the pulse perturbation
  }
  TS<-as.data.frame(lsoda(B, time, ODE_function, params,
                          hini=0.0001, # initial step size
                          rootfun=root_function))#, events = list(func=event_function, root = TRUE))) # solver with adaptive time step
  if(length(time)>1000){
    TS<-TS[seq(1,length(time),floor(length(time)/1000)),] # keeps only 1000 points
  }
  return(TS)
}

### Formatting parameters and results ####
# Create the params list
set_params<-function(params_data_row){ # one line of params_data as argument
  params<-as.list(c(g=params_data_row$g,
                    r=params_data_row$r,
                    D=params_data_row$D,
                    m=params_data_row$m,
                    a=params_data_row$a,
                    e=params_data_row$e,
                    asym=params_data_row$asym,
                    disp=params_data_row$disp,
                    VE=params_data_row$VE,
                    pert=params_data_row$pert,
                    z=params_data_row$z))
  return(params)
}

# Create the output dataframe from the results list
create_data<-function(params_data,results,nSpecies,nCommunity){
  # biomass
  data_B<-as.data.frame(matrix(0,dim(params_data)[1],dim(params_data)[2]+nSpecies*nCommunity))
  names(data_B)[1:dim(params_data)[2]]=names(params_data)
  data_B[,1:dim(params_data)[2]]=params_data
  B_names<-expand.grid(c("B"),c(1:nSpecies),c(1:nCommunity))
  B_names$Var1<-paste(B_names$Var1,B_names$Var2,sep = "_")
  B_names<-paste(B_names$Var1,B_names$Var3,sep = "")
  names(data_B)[(dim(params_data)[2]+1):dim(data_B)[2]]=B_names
  # variance
  data_V<-as.data.frame(matrix(0,dim(params_data)[1],dim(params_data)[2]+(nSpecies*nCommunity)^2))
  names(data_V)[1:dim(params_data)[2]]=names(params_data)
  data_V[,1:dim(params_data)[2]]=params_data
  V_names<-expand.grid(c(1:nSpecies),c(1:nCommunity))
  V_names<-paste(V_names$Var1,V_names$Var2,sep = "")
  V_names<-expand.grid(c("V"),V_names,V_names)
  V_names<-paste(V_names$Var1,V_names$Var2,V_names$Var3,sep = "_")
  names(data_V)[(dim(params_data)[2]+1):dim(data_V)[2]]=V_names
  # correlation
  data_C<-as.data.frame(matrix(0,dim(params_data)[1],dim(params_data)[2]+(nSpecies*nCommunity)^2))
  names(data_C)[1:dim(params_data)[2]]=names(params_data)
  data_C[,1:dim(params_data)[2]]=params_data
  C_names<-expand.grid(c(1:nSpecies),c(1:nCommunity))
  C_names<-paste(C_names$Var1,C_names$Var2,sep = "")
  C_names<-expand.grid(c("C"),C_names,C_names)
  C_names<-paste(C_names$Var1,C_names$Var2,C_names$Var3,sep = "_")
  names(data_C)[(dim(params_data)[2]+1):dim(data_C)[2]]=C_names
  # relative importance of dispersal M1
  data_M1<-as.data.frame(matrix(0,dim(params_data)[1],dim(params_data)[2]+nSpecies*nCommunity))
  names(data_M1)[1:dim(params_data)[2]]=names(params_data)
  data_M1[,1:dim(params_data)[2]]=params_data
  M1_names<-expand.grid(c("M"),c(1:nSpecies),c(1:nCommunity))
  M1_names$Var1<-paste(M1_names$Var1,M1_names$Var2,sep = "_")
  M1_names<-paste(M1_names$Var1,M1_names$Var3,sep = "")
  names(data_M1)[(dim(params_data)[2]+1):dim(data_M1)[2]]=M1_names
  # biomass CV of each population
  data_CV_pop<-as.data.frame(matrix(0,dim(params_data)[1],dim(params_data)[2]+nSpecies*nCommunity))
  names(data_CV_pop)[1:dim(params_data)[2]]=names(params_data)
  data_CV_pop[,1:dim(params_data)[2]]=params_data
  CV_pop_names<-expand.grid(c("CV"),c(1:nSpecies),c(1:nCommunity))
  CV_pop_names$Var1<-paste(CV_pop_names$Var1,CV_pop_names$Var2,sep = "_")
  CV_pop_names<-paste(CV_pop_names$Var1,CV_pop_names$Var3,sep = "")
  names(data_CV_pop)[(dim(params_data)[2]+1):dim(data_CV_pop)[2]]=CV_pop_names
  # biomass CV of each species at regional scale
  data_CV_species<-as.data.frame(matrix(0,dim(params_data)[1],dim(params_data)[2]+nSpecies))
  names(data_CV_species)[1:dim(params_data)[2]]=names(params_data)
  data_CV_species[,1:dim(params_data)[2]]=params_data
  CV_species_names<-paste("CV_",c(1:nSpecies))
  names(data_CV_species)[(dim(params_data)[2]+1):dim(data_CV_species)[2]]=CV_species_names
  # aggregated CV
  data_CV<-params_data
  data_CV$CV_pop=0 # average biomass CV
  data_CV$CV_tot=0 # total biomass CV
  # response to a press perturbation
  data_P<-as.data.frame(matrix(0,dim(params_data)[1],dim(params_data)[2]+(nSpecies*nCommunity)^2))
  names(data_P)[1:dim(params_data)[2]]=names(params_data)
  data_P[,1:dim(params_data)[2]]=params_data
  P_names<-expand.grid(c(1:nSpecies),c(1:nCommunity))
  P_names<-paste(P_names$Var1,P_names$Var2,sep = "")
  P_names<-expand.grid(c("P"),P_names,P_names)
  P_names<-paste(P_names$Var1,P_names$Var2,P_names$Var3,sep = "_")
  names(data_P)[(dim(params_data)[2]+1):dim(data_P)[2]]=P_names
  # asymptotic resilience
  data_resilience<-params_data
  data_resilience$resilience=0 # real part of the lead eigne value
  # contribution to the lead eigenvector
  data_E<-as.data.frame(matrix(0,dim(params_data)[1],dim(params_data)[2]+nSpecies*nCommunity))
  names(data_E)[1:dim(params_data)[2]]=names(params_data)
  data_E[,1:dim(params_data)[2]]=params_data
  E_names<-expand.grid(c("E"),c(1:nSpecies),c(1:nCommunity))
  E_names$Var1<-paste(E_names$Var1,E_names$Var2,sep = "_")
  E_names<-paste(E_names$Var1,E_names$Var3,sep = "")
  names(data_E)[(dim(params_data)[2]+1):dim(data_E)[2]]=E_names
  # jacobian matrix
  data_J<-as.data.frame(matrix(0,dim(params_data)[1],dim(params_data)[2]+(nSpecies*nCommunity)^2))
  names(data_J)[1:dim(params_data)[2]]=names(params_data)
  data_J[,1:dim(params_data)[2]]=params_data
  J_names<-expand.grid(c(1:nSpecies),c(1:nCommunity))
  J_names<-paste(J_names$Var1,J_names$Var2,sep = "")
  J_names<-expand.grid(c("J"),J_names,J_names)
  J_names<-paste(J_names$Var1,J_names$Var2,J_names$Var3,sep = "_")
  names(data_J)[(dim(params_data)[2]+1):dim(data_J)[2]]=J_names
  # fill in tha data frames
  for (i in 1:dim(params_data)[1]){
    data_B[i,which(names(data_B)%in%B_names)]<-unlist(results[[i]]$B)
    data_V[i,which(names(data_V)%in%V_names)]<-unlist(results[[i]]$V)
    data_C[i,which(names(data_C)%in%C_names)]<-unlist(results[[i]]$C)
    data_M1[i,which(names(data_M1)%in%M1_names)]<-unlist(results[[i]]$M1)
    data_CV_pop[i,which(names(data_CV_pop)%in%CV_pop_names)]<-unlist(results[[i]]$CV_pop)
    data_CV_species[i,which(names(data_CV_species)%in%CV_species_names)]<-unlist(results[[i]]$CV_species)
    data_CV[i,which(names(data_CV)%in%c("CV_pop","CV_tot"))]=unlist(results[[i]]$CV)
    data_P[i,which(names(data_P)%in%P_names)]<-unlist(results[[i]]$Press)
    data_resilience[i,which(names(data_resilience)%in%c("resilience"))]=unlist(results[[i]]$resilience)
    data_E[i,which(names(data_E)%in%E_names)]<-unlist(results[[i]]$E)
    data_J[i,which(names(data_J)%in%J_names)]<-unlist(results[[i]]$jacobian)
  }
  return(list(data_B=data_B,
              data_V=data_V,
              data_C=data_C,
              data_M1=data_M1,
              data_CV_pop=data_CV_pop,
              data_CV_species=data_CV_species,
              data_CV=data_CV,
              data_P=data_P,
              data_resilience=data_resilience,
              data_E=data_E,
              data_J=data_J,
              B_names=B_names,
              V_names=V_names,
              C_names=C_names,
              M1_names=M1_names,
              CV_pop_names=CV_pop_names,
              CV_species_names=CV_species_names,
              P_names=P_names,
              E_names=E_names,
              J_names=J_names))
}

# Create the output dataframe from the results list containing the time series
create_TS<-function(params_data, results, B0){
  TS<-NULL
  for(i in 1:dim(params_data)[1]){
    databis<-results[[i]]
    if(is.null(B0)==F){
      for(j in 1:dim(databis)[1]){
        databis[j,2:dim(databis)[2]] = databis[j,2:dim(databis)[2]]/B0[i,] # rescaling with B0
      }
    }
    databis<-merge(params_data[i,],databis)
    TS<-rbind(TS,databis)
  }
  return(TS)
}

# Make a table to plot a correlation matrix
table_for_matrix<-function(table,nparams){
  table<-melt(table,
              id.vars = names(table)[1:nparams],
              variable.name = "species",
              value.name = "correlation")
  table<-table %>% separate(species,c(NA,"species_1","species_2"),sep="_")
  table$species_1<-as.factor(table$species_1)
  table$species_2<-as.factor(table$species_2)
  table<-table %>% separate(species_1,c("species_1","community_1"),sep=1)
  table<-table %>% separate(species_2,c("species_2","community_2"),sep=1)
  return(table)
}

# Make a table ready to use for ggplot
table_for_plot<-function(table,nparams,value.name){
  table<-melt(table,
              id.vars = names(table)[1:nparams],
              variable.name = "species",
              value.name = value.name)
  table<-table %>% separate(species,into=c(NA,"species"),sep="_")
  table<-table %>% separate(species,into=c("species","community"),sep=c(1))
  return(table)
}

# set asymmetry for two patches
set_asymmetry<-function(nSpecies,gamma,omega){
  asym_1<-c(omega,rep(gamma,nSpecies-1)) # patch with modified growth rate and attack rate
  asym_2<-rep(1,nSpecies) # reference patch
  return(list(list(asym_1,asym_2)))
}

# return an aggregated dataframe ready to plot
time_series_for_plot<-function(TS, n_params){
  TS<-melt(TS,
           id.vars = names(TS)[1:(1+n_params)],
           variable.name = "species",
           value.name = "biomass")
  TS<-TS %>% separate(species,into=c(NA, NA, "species","community"),sep=c(1,2,3))
  return(TS)
}

# get the params list
get_params<-function(params_data,i){
  params<-as.list(c(g=params_data$g[i],
                    r=params_data$r[i],
                    D=params_data$D[i],
                    m=params_data$m[i],
                    a=params_data$a[i],
                    e=params_data$e[i],
                    gamma=params_data$gamma[i],
                    omega=params_data$omega[i]))
  params$asym=params_data$asym[[i]]
  params$disp=params_data$disp[[i]]
  return(params)
}
