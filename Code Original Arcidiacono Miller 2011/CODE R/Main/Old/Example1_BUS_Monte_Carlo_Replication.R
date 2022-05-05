##########################################################################
#from "intconddata"
#This program recreates Table 1 columns 2,3,5, and 6 from
#"Conditional Choice Probability Estimation of Dynamic Discrete Choice Models 
#with Unobserved Heterogeneity" by Arcidiacono and Miller (2011)

#The original code is by the above authors

#Coversion to R was done by Wayne Taylor

#Version 1/5/2016

#Note: The intent of this code was to remain consistent with the original code in terms of structure and naming conventions
#Therefore, the code has not been optimtized for speed

################################################################################

################################################################################
# This code is from Wayne Taylor
# Link to the Github : https://github.com/waynejtaylor/Single-Agent-Dynamic-Choice/blob/master/AM2011Table1cols2356.R
# The intent is to replicate the results of Arcidiacono and Miller (2011)
# This is a conversion to R of the authors code that is available on Matlab
# Link to authors' code : https://www.econometricsociety.org/publications/econometrica/2011/11/01/conditional-choice-probability-estimation-dynamic-discrete
# Wayne's code replicate Table 1 columns 2,3,5 and 6

# Maxime Brun, Loïc Cantin and Thibauld Ingrand
# we removed FIML results so that the code only replicates columns 3 and 6 
# We add comments to the code

###############################################################################

library(Rcpp) # Library used to read code in C++ from R

#path = 'D:/2021-22 Academic Year/IPP/Advanced Econometrics - Simulation Methods/'
path = 'C:/Users/33689/Documents/GitHub/Exam-Advanced-Econometrics-Simulation-Methods/Code Original Arcidiacono Miller 2011/CODE R/Files/'

source(paste(path,'xgrid.R',sep="")) # Function used to compute the transition matrix
source(paste(path,'wlogitd.R',sep="")) # Function used for reduced form logit
source(paste(path,'wlogit.R',sep="")) # Function maximized by the structural parameters
sourceCpp(paste(path,'genbus4.cpp',sep="")) # Function used to generate the data
sourceCpp(paste(path,'fvdataBOTH.cpp',sep="")) # Function used to compute value function part of the likelihood
source(paste(path,'intcond.R',sep="")) # Used to compute initial conditions
source(paste(path,'intcondP.R',sep="")) # Used to compute q_ns


set.seed(1)

# Specify the real parameters of the Data Generating Process (DGP)
alpha=c(2,-.15,1,.9,.4) #Intercept (theta_0), mileage (theta_1), heterogeneity (theta_2), discount factor (beta), Pi

# Tolerance of the EM Algorithm
tol=.001 

# Code Specification
MCiter=1       # Number of Monte Carlo iterations
hetero = TRUE   # Is heterogeneity unobserved? FALSE = column 3 ; TRUE = column 6
T=200          # Time periods
if(hetero) T=T/10 # If heterogeneity is unobserved, reduce the time period
N=1000        # Observations per time period

# Ouput vectors General
Bccp=NULL #CCP parameter storage
Tccp=NULL #CCP timing

# Additional output vectors (when we allow for unobserved heterogeneity)
if(hetero){
  Lccp=NULL
  LFl=NULL
  Iccp=NULL
  Binit=NULL
}

# Create transition matrices
zval=seq(.25,1.25,.01) # support of x2
zbin=length(zval) # cardinal of x2 support
xval=seq(0,25,.125) # support of x1
xbin=length(xval) # cardinal of x1 support
xtran=matrix(0,zbin*xbin,xbin) # to be filled with transition probabilities
xtranc=array(0,c(xbin,xbin,zbin)) # to be filled with transition cumulative distribution function (CDF) 
for(z in 1:zbin){
  temp=xgrid(zval[z],xval)  # Call xgrid
  xtran[(1+(z-1)*xbin):(z*xbin),] = temp$xtran # store the transition probability
  xtranc[,,z] = temp$xtranc # store the transition CDF
}

# In the RCCP format
xtrancRcpp = matrix(xtranc,xbin,xbin*zbin)

# State space cardinal
tbin=xbin*zbin # 201*101

# Create vector of z and x values for each state
zvalr=kronecker(zval,rep(1,xbin)) # repeat xbin times zval
xvalr=kronecker(rep(1,zbin),xval)/10 # repeat zbin times xval

# Regressors for the reduced form logits which covers the state space
RX1=cbind(rep(1,zbin*xbin),xvalr,zvalr,xvalr*zvalr,xvalr*xvalr,zvalr*zvalr) # W1t in Supplemental (p. 12/17)

## Monte Carlo simulations

# Starting values for CCP
alphaf= c(alpha[1:3],log(alpha[4])-log(1-alpha[4]))
alphac= alpha[1:4] # theta_0,theta_1,theta_2,beta

MC=1
while(MC <= MCiter){
  
  # Generating the data 
  genbusout=genbusRcpp(alpha,N,T,xtran,xtrancRcpp,xbin,zbin,xval,zval) # Simulate the data
  
  Y=genbusout$Y # Decisions (1 is keep, 0 is replace) : d2_nt => N,T matrix
  X=genbusout$X # Mileage : x1_nt => N,T matrix 
  Z=genbusout$Z # Permanent route characteristic : x2_n => N,1 matrix 
  State=genbusout$State # Brand : s_n => N,1 matrix
  FVT=genbusout$FVT # Value function : V_nt => N,t matrix
  
  Xstate=genbusout$Xstate # Position of X in xval
  Zstate=genbusout$Zstate # Position of Z in zval
  
  
  # Vectorize the simulated data in the format N*T (when s is observed)
  y2=as.vector(Y)
  x2=as.vector(X[,1:T])/10 # Why divide by 10??
  z2=kronecker(rep(1,T),Z)
  s2=kronecker(rep(1,T),State)
  t2=kronecker(1:T,rep(1,N))/10 # Why divide by 10?
  
  # Vectorize when s in unobserved : format N*T*2
  if(hetero){
    y2=c(y2,y2)
    x2=c(x2,x2)
    z2=c(z2,z2)
    s2=c(rep(0,N*T),rep(1,N*T)) # Brand is restated because unobserved
    t2=c(t2,t2)
    stemp=c(rep(0,N),rep(1,N))
  }
  
  # Estimating with data CCPs
  
  tic = proc.time()[3] #start the timer
   
  # Setting up data for reduced form logit : interation between W1t and W2t (see Supplemental p. 13/17)
  xx=cbind(rep(1,N*T),x2,z2,x2*z2,x2*x2,z2*z2,s2,s2*x2,s2*z2,s2*x2*z2,s2*x2*x2,s2*z2*z2) # 12 regressors
  xx=cbind(xx,matrix(rep(t2,12),ncol=12)*xx,matrix(rep(t2*t2,12),ncol=12)*xx) # 36 regressors
  
  # Estimating reduced form logit
  # Obtain b1 : a set of parameters that will be used to estimate p1.
  if(!hetero){ # if the brand is observed
    b1=glm(((y2==0)*1)~xx-1,family='binomial')$coef
  } else { # if the brand is unobserved
    PType=.5*rep(1,2*N*T)
    oPType=rep(0,2*N*T)
    Pi2=c(.5,.5)  # Prior probability
    
    b1 = rep(0,ncol(xx))
    b1 = optim(b1,wlogitd,Y=(y2==0),X=xx,P=PType,method="BFGS")$par # coefficient estimates
    # For a binomial GLM prior weights are used to give the number of trials when the response is the proportion of successes
    # So we cannot send them into the "weights" argument
  }
  
  # Calculating fv terms
  # fvt1 : value function part of likelihood -> the term multiplied by beta
  if(!hetero){ # s observed
    fvt1 = fvdataRcpp(b1,RX1,tbin,xbin,Zstate,Xstate,xtran,N,T,State)  
  } else { # s unobserved
    fvt1 = fvdataRcpp(b1,RX1,tbin,xbin,Zstate,Xstate,xtran,N,T,rep(1,N),hetero)  
  }
  
  # Estimating the structural parameters
  xccp = cbind(rep(1,N*T),x2*10,s2) # Regressors for theta_0,theta_1,theta_2
  
  if(!hetero){ # s observed
    bccp = glm(y2~cbind(xccp,fvt1)-1,family='binomial')$coef # Use a GLM
  } else { # s unobserved --> EM ALGORIHM
    
    # When s in unobserved, we need the EM algorithm
    # Starting the EM algorithm
    
    j=0
    
    bccp = alphac # Initial parameter values
    
    intcondX=cbind(rep(1,N), X[1:N,1],Z[1:N,1]) # Regressors at time 1
    binit=rep(0,3) # Initial value of b1
    
    cond=0 
    lp=NULL
    while(cond==0){
      
      # EXPECTATION
      
      # Updating PType
      ## first getting the type-specific likelihoods
      
      oPType=PType
      
      # Replaces the authors' Matlab function "likeCPP"
      U1 = cbind(xccp,fvt1)%*%bccp # is v2-v1
      Like = (y2*exp(U1)+(1-y2))/(1+exp(U1)) # N*T*2 vector
      
      Like2=array(Like,c(N,T,2)) # Likelihood contributions by bus, time and brand : N*T*S
      base=apply(Like2,c(1,3),prod) # Likelihood contribution by bus and brand : N*S
      
      # UPDATE the m-th pi(x|x)
      # Now getting the initial condition parameters (see supplemental B.1.4)
      
      intcond_optim=optim(binit,intcond,like=base,X=intcondX,method="BFGS") 
      binit = intcond_optim$par # Obtain delta(m+1) => to be used to update the m-th pi(s|x)
      lp=c(lp,intcond_optim$value) # Value of the likelihood
      
      # UPDATE the m-th q_ns
      # And the PType's
      
      PType=intcondP(binit,base,intcondX) # update q_ns : N bus * 2 states
      PType=kronecker(rep(1,T),PType) # q_nst : N bus * 2 states * T periods
      PType=as.vector(PType) # q_nst as a vector
      
      # UPDATE the m_th p1(x,s)
      # Estimating reduced form logit (See supplemental B.1.3)
      b1 = optim(b1,wlogitd,Y=(y2==0),X=xx,P=PType,method="BFGS")$par # y2==0 <=> d1t=1. See supplemental B.1.3.
      # b1 will be used to obtain p1 in the likelihood
      
      
      # MAXIMIZATION
      # Calculating fv terms --> The term multiplied by beta that involes p1
      fvt1=fvdataRcpp(b1,RX1,tbin,xbin,Zstate,Xstate,xtran,N,T,rep(1,N),hetero) #  Value function part of the likelihood
      
      # Structural parameters
      bccp = optim(bccp,wlogit,Y=y2,X=cbind(xccp,fvt1),P=PType,method="BFGS")$par # Last step of the algorithm
      
      #CHECKING CONVERGENCE
      if(j>26){
        junk=abs((lp[j]-lp[j-25])/lp[j])<tol
        junk2=abs((lp[j-1]-lp[j-26])/lp[j-1])<tol
        
        cond=junk2*junk
        
        if(j>1000) cond=1
      }
      
      j=j+1
      cat("j: ",j,fill=TRUE)
    }
  }
  
  toc = proc.time()[3]-tic
  
  Tccp=c(Tccp,toc) # Store the computation time
  Bccp=rbind(Bccp,bccp) # Store the coefficient estimates
  
  if(hetero){
    Iccp=c(Iccp,j)
    Binit=rbind(Binit,binit)  
  }
  
  cat("MC ",MC, " completed",fill=TRUE)
  MC = MC+1
}
Bccp
