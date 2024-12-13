#function for simulating different parasite-mediated anorexia scenario
#Output is two-layered list
#Each element of the first layer is a replicate
#Each element of the second layer is a time point
simulation = function(rep,perc=1,feedback=T,tm,constFT=NA,SDresid,
                      indM,ParFEC,simindpar,r,K,MIT,FP,
                      suri,sure,beta,shed,N0,A0,EP0,FEC0low=NA,FEC0up=NA,core){
  #rep: # replicates; perc: percent herbivore individuals showing changing feeding rate (FR)
  #feedback: feedback effect of feeding change on parasite ingestion? tm: maximum time step (30days)
  #constFT: parasite ingestion rate when no feedback effect
  #SDresid: standard deviation of feeding rate deviant from a GLM prediction
  #r: plant growth rate; K: plant carrying capacity; MIT: maximum intake rate
  #FP: fecal production; suri: parasite survival within host; 
  #sure: parasite survival in the environment
  #beta: transmission probability given 100% of encounter time; N0: # herbivore individuals
  #A0: initial plant biomass; EP0: initial environmental parasite abundance
  #FEC0low and FEC0up: lower and upper bounds of initial fecal egg count if exist
  #core: # cores used for simulation
  
  #Par: parameter values for FEC, Intp, vegSD, vegSD-slope and GLM
  library(parallel);library(doParallel);library(tidyverse)
  
  #function for feeding rate prediction
  FT.pred = function(M,FECvalue,SDresid=NA){
    pred = data.frame(FEC=FECvalue) %>% 
      data.frame(mean=predict(M,.,se=T,type="response")$fit,
                 se=predict(M,.,se=T,type="response")$se.fit) %>% 
      select(mean,se)
    if(!is.na(SDresid)){
      pred$se=SDresid
    }
    out = rnorm(1,pred$mean,pred$se)
    #print(out)
    if(out<0){out=0}
    if(out>1){out=1}
    return(out)
  }
  
  #function for generating individual intercept and coefficient for the anorexic effect
  rnglm = function(min=1,max=9){
    repeat{a = ceiling(runif(n=1,min=(min-1),max=max))
    if(a>=min&a<=max){break}}
    return(a)
  }
  
  cat(paste0("Start: ",Sys.time(),"  "))
  registerDoParallel(cores=core)
  

  #Beginning of simulation
  sim1 = foreach(n=1:rep,.packages=c("tidyverse")) %dopar% {
    sim2 = list()
    #create starting herbivore individuals
    set.seed(round(runif(1,1,1e5)))
    N = N0
    ind = vector(mode="list", N+2) #+2 for environmental parasite and plant
    
    #plant biomass
    ind[[1]] = c(NA,NA,NA,NA,NA,NA,NA,NA,NA,A0)
    
    #environmental parasite
    ind[[2]] = c(NA,NA,NA,NA,NA,NA,NA,NA,EP0,NA)
    
    Model = list()
    #determine individual characteristics
    #print(2)
    if(perc!=0){
      for(i in 3:(round(N*perc)+2)){
        nsimindpar=rnglm(min=1,max=nrow(simindpar))
        Intp = simindpar[nsimindpar,"Intp"]
        Slope = simindpar[nsimindpar,"Slope"] 
        if(is.na(FEC0low)){
          FEC = rnbinom(1,size=ParFEC[1],mu=ParFEC[2])
        }
        if(!is.na(FEC0low)){
          lowFEC = qnbinom(FEC0low,size=ParFEC[1],mu=ParFEC[2])
          upFEC = qnbinom(FEC0up,size=ParFEC[1],mu=ParFEC[2])
          repeat{
            FEC = rnbinom(1,size=ParFEC[1],mu=ParFEC[2])
            if(FEC>=lowFEC&FEC<=upFEC){break}
          }
        }
        nGLM = rnglm(min=1,max=9)
        Model[[i]] = indM[[nGLM]]
        Model[[i]]$coefficients=c(Intp,Slope)
        FT = FT.pred(M=Model[[i]],FECvalue=FEC,SDresid=SDresid)
        Alive = 1
        Age = 0
        ind[[i]]= c(Alive,Age,nGLM,Intp,Slope,FEC,1,FT,NA,NA)
        #c(1.alive,2.age,3.model,4.intercept,5.coefficient,6.fecal egg count,
        #  7.feeding rate change?,8.foraging time,9.environmental parasite,10.plant) 
      }
    }
    if(perc!=1){
      for(i in (round(N*perc)+3):(N+2)){
        nsimindpar=rnglm(min=1,max=nrow(simindpar))
        Intp = simindpar[nsimindpar,"Intp"]
        Slope = NA
        if(is.na(FEC0low)){
          FEC = rnbinom(1,size=ParFEC[1],mu=ParFEC[2])
        }
        if(!is.na(FEC0low)){
          lowFEC = qnbinom(FEC0low,size=ParFEC[1],mu=ParFEC[2])
          upFEC = qnbinom(FEC0up,size=ParFEC[1],mu=ParFEC[2])
          repeat{
            FEC = rnbinom(1,size=ParFEC[1],mu=ParFEC[2])
            if(FEC>=lowFEC&FEC<=upFEC){break}
          }
        }
        nGLM = rnglm(min=1,max=9)
        Model[[i]] = indM[[nGLM]]
        Model[[i]]$coefficients[1]=Intp
        FT = FT.pred(M=Model[[i]],FECvalue=0,SDresid=SDresid)
        Alive = 1
        Age = 0
        ind[[i]]= c(Alive,Age,nGLM,Intp,Slope,FEC,0,FT,NA,NA)
        #c(1.alive,2.age,3.model,4.intercept,5.coefficient,6.fecal egg count,
        #  7.feeding rate change?,8.foraging time,9.environmental parasite,10.plant) 
      }
    }
    
    #make empty vectors to record population statistics
    #print(3)
    time = seq(tm*30+1)
    
    Out = vector(mode="list", (tm*30)+1)
    Out[[1]] = data.frame(ID=rep(NA,N+2),Alive=rep(NA,N+2),Age=rep(NA,N+2),GLM=rep(NA,N+2),
                          Intp=rep(NA,N+2),Slope=rep(NA,N+2),FEC=rep(NA,N+2),
                          Change=rep(NA,N+2),FT=rep(NA,N+2),EP=rep(NA,N+2),A=rep(NA,N+2))
    
    for(i in 1:(N+2)){
      if(i>2){
        Out[[1]]$ID[i]  = i-2  
      }
      Out[[1]]$Alive[i]  = ind[[i]][1]
      Out[[1]]$Age[i]  = ind[[i]][2]
      Out[[1]]$GLM[i]  = ind[[i]][3]
      Out[[1]]$Intp[i]  = ind[[i]][4]
      Out[[1]]$Slope[i]  = ind[[i]][5]
      Out[[1]]$FEC[i] = ind[[i]][6]
      Out[[1]]$Change[i] = ind[[i]][7]
      Out[[1]]$FT[i] = ind[[i]][8]
      Out[[1]]$EP[i] = ind[[i]][9]
      Out[[1]]$A[i]  = ind[[i]][10] 
    }
    
    for(i in 1:(tm*30)){ #loop through timesteps
      #ingested parasite
      ingEP=rep(NA,(N+2))
      for(k in 3:(N+2)){
        if(feedback==T){
          ingEP[k]=beta*ind[[2]][9]*ind[[k]][8]/FP 
        }else{
          ingEP[k]=beta*ind[[2]][9]*constFT/FP
        }
      }
      sumNi = Reduce('+',lapply(ind, function(x) replace(x,is.na(x),0)))[1]
      sumFECi = Reduce('+',lapply(ind, function(x) replace(x,is.na(x),0)))[6]
      sumFTi = Reduce('+',lapply(ind, function(x) replace(x,is.na(x),0)))[8]
      
      #plant
      ind[[1]][10] = max(0,(ind[[1]][10] + r*(1-ind[[1]][10]/K)*ind[[1]][10] - sumFTi*MIT/0.8)) 

      #environmental parasite
      ind[[2]][9] = ind[[2]][9]*sure + sumFECi*115 - sum(ingEP,na.rm=T)*115
      
      #herbivores      
      for(j in 3:(N+2)){ #loop for each individual
        if(ind[[j]][1]==0){ NULL }
        if(ind[[j]][1]==1){
          #Age
          ind[[j]][2] = ind[[j]][2]+1
          #FEC
          ind[[j]][6] = sum(ind[[j]][6]*suri,ingEP[j])
          #FT
          if(ind[[j]][7]==1){
            ind[[j]][8] = FT.pred(M=Model[[j]],FECvalue=ind[[j]][6],SDresid=SDresid)
            
          }
          if(ind[[j]][7]==0){
            ind[[j]][8] = FT.pred(M=Model[[j]],FECvalue=0,SDresid=SDresid)
          }
        }
      }
      
      Out[[i+1]] = data.frame(ID=rep(NA,N+2),Alive=rep(NA,N+2),Age=rep(NA,N+2),GLM=rep(NA,N+2),
                              Intp=rep(NA,N+2),Slope=rep(NA,N+2),FEC=rep(NA,N+2),
                              Change=rep(NA,N+2),FT=rep(NA,N+2),EP=rep(NA,N+2),A=rep(NA,N+2))
      
      #Population stats
      for(k in 1:(N+2)){
        if(k > 2){
          Out[[i+1]]$ID[k]  = k-2  
        }
        Out[[i+1]]$Alive[k]  = ind[[k]][1]
        Out[[i+1]]$Age[k]  = ind[[k]][2]
        Out[[i+1]]$GLM[k]  = ind[[k]][3]
        Out[[i+1]]$Intp[k]  = ind[[k]][4]
        Out[[i+1]]$Slope[k] = ind[[k]][5]
        Out[[i+1]]$FEC[k] = ind[[k]][6]
        Out[[i+1]]$Change[k] = ind[[k]][7]
        Out[[i+1]]$FT[k] = ind[[k]][8]
        Out[[i+1]]$EP[k] = ind[[k]][9]
        Out[[i+1]]$A[k]  = ind[[k]][10] 
      }
    }
    sim2 = Out
    return(sim2)
  }
  sim = sim1
  cat(paste0("End: ",Sys.time(),"\n"))
  beepr::beep(2)
  return(sim)
}  

load("input.RData")
#Four objects
#indM: a list of the nine GLMs from the nine gazelle empirical data
#simindpar: a data frame including 100,000 combinations of GLM intercepts 
#and coefficients randomly generated by package SimMultiCorrData
#ParFEC: parameters for a negative binomial for initial fecal egg count
#SDresid: average standard deviation of feeding rate deviant from a GLM prediction

r = 1.08      #plant growth rate
K = 200       #plant carrying capacity
MIT = 0.86    #herbivore intake per time unit
FP = 230

suri = 0.979  #parasite survival inside hosts
sure = 0.917  #parasite survival in the environment
beta = 1e-4   #probability of infection

N0 = 100      #number of herbivore individuals
A0 = K*1.0    #plant starting biomass
EP0 = 1e8   #starting environmental parasite
tm = 12

example = simulation(rep=20,perc=1,feedback=T,tm=tm,N0=N0,A0=A0,EP0=EP0,
                     SDresid=SDresid,indM=indM,ParFEC=ParFEC,simindpar=simindpar,
                     r=r,K=K,MIT=MIT,FP=FP,suri=suri,sure=sure,beta=beta,core=20)
length(example) #20 elements because of 20 replicates
length(example[[1]]) #361 elements per replicate because of timestep 0-360
example[[1]][[1]]
#102 rows including plants, environmental parasites and 100 herbivores
#columns: 1.alive or not (for counting),2.age,3.model number from the nine original GLMs,
#4.intercept,5.coefficient,6.fecal egg count,7.feeding rate change?,8.foraging time,
#9.environmental parasite,10.plant