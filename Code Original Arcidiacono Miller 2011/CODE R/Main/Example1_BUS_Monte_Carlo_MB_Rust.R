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

##########################################################################


library(Rcpp)
library(readxl)
library(zoo) # To retreat data
library(dplyr)
library(reshape2)
library(truncnorm)


# Set your path up to Exam-Advanced-Econometrics-Simulation-Methods
#path = 'D:/2021-22 Academic Year/IPP/Advanced Econometrics - Simulation Methods/'
path = 'C:/Users/33689/Documents/GitHub/Exam-Advanced-Econometrics-Simulation-Methods/'

source(paste(path,'Code Original Arcidiacono Miller 2011/CODE R/Files/xgrid_Rust.R',sep=""))
source(paste(path,'Code Original Arcidiacono Miller 2011/CODE R/Files/wlogitd.R',sep=""))
source(paste(path,'Code Original Arcidiacono Miller 2011/CODE R/Files/wlogit.R',sep=""))
source(paste(path,'Code Original Arcidiacono Miller 2011/CODE R/Files/likebusML4.R',sep=""))
#sourceCpp(paste(path,'genbus4.cpp',sep="")) not necessary for estimation
sourceCpp(paste(path,'Code Original Arcidiacono Miller 2011/CODE R/Files/fvdataBOTH_Rust.cpp',sep=""))
source(paste(path,'Code Original Arcidiacono Miller 2011/CODE R/Files/intcond.R',sep=""))
source(paste(path,'Code Original Arcidiacono Miller 2011/CODE R/Files/intcondP.R',sep=""))

# Load data
df_rust <- read_excel(paste(path,'Rust Data/rust-data.xlsx',sep=""))

# Retreat NA values
df_rust$Bus_ID <- na.locf(df_rust$Bus_ID) # Fill NA values with Bus ID
df_rust$mileage[is.na(df_rust$mileage)] <- 0 # In column mileage, replace NA with 0
df_rust$usage[is.na(df_rust$usage)] <- 0 # Same for column usage


# Round mileage
df_rust$mileage <- floor(df_rust$mileage/1000) # express mileage in thousands

# Modify the decision values to be consistent with the original code
df_rust$decision <- df_rust$decision+1 # Replace decision values 0 -> 1
df_rust$decision[df_rust$decision==2] <- 0 # Replace decision values 2 -> 0

# For each bus, compute the total number of observations available and if there is a replacement
df_rust <- df_rust %>%
  group_by(Bus_ID) %>% # for each bus 
  mutate(n_periods = n()) %>% # store the number of observations in n_periods
  mutate(remplacement = sum(decision)<n()) %>%
  ungroup()

# Keep bus with replacement
df_rust <- df_rust[df_rust$remplacement==TRUE,]

# Retreat type
df_rust$s=(df_rust$type=="5308A-75")+(df_rust$type=="5308A-74")+(df_rust$type=="5308A-72")


# Define the period and duration for the estimation
t_start=5 # first period to be considered
T=100 # number of periods we take into account
t_min=t_start+T+5 # the bus must be operated for at least for t_min period

# Filter the data frame to keep only the buses that match our period selection
df_rust <- df_rust[df_rust$n_periods>=t_min,] # The bus must be operated for at least t_min periods
df_rust <- df_rust[df_rust$period>=t_start,] # Keep only the periods above t_start
df_rust <- df_rust[df_rust$period<=t_start+T-1,] # Keep only the T periods after t_start

set.seed(1)

#Initial parameter values
# #Intercept (theta_0), mileage (theta_1), heterogeneity (theta_2), discount factor (beta), Pi
# Fonctionne : alpha=c(30,-0.001,4,0.9,.4)
#alpha=c(7.21661496,-0.01632606,1,0.01104702,.4)

alpha=c(2,-0.01,4,0.9,.5)


# tol=.0000001
tol=.000001

FIML = FALSE   #estimate FIML too? (it takes much longer than CCP)
hetero = TRUE #Is heterogeneity observed? FALSE = cols 1 and 2 TRUE = cols 5 and 6

# Bccp=NULL #CCP parameter storage
# Tccp=NULL #CCP timing
# Bfl=NULL  #FIML parameter storage
# Tfl=NULL  #FIML timing

if(hetero){
  Lccp=NULL
  LFl=NULL
  Iccp=NULL
  Binit=NULL
}


# Support for x1 and x2

zval=0.577216 # In Rust, x2 is deterministic : x2=0.577216 (see Rust 4.11)
zbin=length(zval) # cardinal of x2 support : only 1 value

xval=seq(0,max(df_rust$mileage),1) # support of x1
xbin=length(xval) # cardinal of x1 support

#Create transition matrices with x_grid

xtran=matrix(0,zbin*xbin,xbin) # to store the probability
xtranc=array(0,c(xbin,xbin,zbin)) # to store the cumulative distribution
for(z in 1:zbin){ # Here we have only 1 value for z
  temp=xgrid_Rust(zval[z],xval)  # call xgrid function
  xtran[(1+(z-1)*xbin):(z*xbin),] = temp$xtran # store probability transition
  xtranc[,,z] = temp$xtranc # store cumulated transition
}

# Check the transition matrices
length(xtran[1,])==xbin # xtran should have xbin columns
sum(xtran[1,])==1 # xtran row should sum to 1
xtranc[1,xbin,1]==1 # xtranc last column should be equal to 1

# Rcpp
xtrancRcpp = matrix(xtranc,xbin,xbin*zbin)
tbin=xbin*zbin # equal to xbin here because zbin=1

#z and x values for each state
zvalr=kronecker(zval,rep(1,xbin)) # we repeat xbin times xval
xvalr=kronecker(rep(1,zbin),xval)/10 # equal to xval because only 1 zval

#data for reduced form logits
#covers the state space
#RX1=cbind(rep(1,zbin*xbin),xvalr,zvalr,xvalr*zvalr,xvalr*xvalr,zvalr*zvalr) # corresponds to W1t in supplemental material

RX1=cbind(rep(1,xbin),xvalr,xvalr*xvalr) # corresponds to W1t in supplemental material


#starting values for FIML and CCP
alphaf= c(alpha[1:3],log(alpha[4])-log(1-alpha[4]))
alphac= alpha[1:4]

# Format our data as the output of genbus4 in the case of simulations

Y=dcast(data = df_rust,formula = Bus_ID~period,fun.aggregate = sum,value.var = "decision") # Pivot the dataframe for decision
Y= as.matrix(select(Y, -Bus_ID)) # Remove the column Bus_ID and format as matrix

X = dcast(data = df_rust,formula = Bus_ID~period,fun.aggregate = sum,value.var = "mileage") # Pivot the dataframe for mileage
X = as.matrix(select(X, -Bus_ID)) # Remove the column Bus_ID and format as matrix

State = dcast(data = df_rust,formula = Bus_ID~period,fun.aggregate = sum,value.var = "s") # Pivot the dataframe for decision
State= as.matrix(select(State, -Bus_ID))

N=length(X[,1]) # Number of buses in our data set
T=length(X[1,]) # Number of periods in our data set

Z=replicate(N,1)*zval # give value zval to each bus
Zstate=replicate(N,1) # the value of Z is alway zval[1] => position 1



Xstate= matrix(0,nrow=N,ncol=T)
for (i in 1:length(X[1,])){
  Xstate[,i]=match(X[,i],xval) # position of X value in the vector xval
}

# Vectorize 
y2=as.vector(Y)
x2=as.vector(X[,1:T])
z2=kronecker(rep(1,T),Z)
#s2=kronecker(rep(1,T),State)
s2=as.vector(State[,1:T])
t2=kronecker(1:T,rep(1,N)) # replicate the time vector for each bus

if(hetero){
  y2=c(y2,y2) # decision
  x2=c(x2,x2) # mileage
  z2=c(z2,z2) # x2
  s2=c(rep(0,N*T),rep(1,N*T)) # restated (because s is unobserved by the econometrician)
  t2=c(t2,t2) 
  stemp=c(rep(0,N),rep(1,N)) # temporary value for s
}

#estimating FIML----
if(FIML){
  
  tic = proc.time()[3] #start the timer
  
  if(!hetero){
    bfl=optim(alphaf,likebusML4,Y=Y,State=s2,N=N,T=T,X=X,Zstate=Zstate,Xstate=Xstate,xtran=xtran,tbin=tbin,zbin=zbin,xbin=xbin,xval=xval,Z=Z)
  } else {
    bfl=optim(alphaf,likebusML4,Y=y2,State=stemp,N=N,T=T,X=x2,Zstate=c(Zstate,Zstate),Xstate=c(Xstate,Xstate),xtran=xtran,tbin=tbin,zbin=zbin,xbin=xbin,xval=xval,Z=Z)
  }
  
  toc = proc.time()[3]-tic
  
  Tfl=c(Tfl,toc)
  Bfl=rbind(Bfl,bfl)
}


#estimating with data ccps----

tic = proc.time()[3] #start the timer

#setting up data for reduced form logit
# xx=cbind(rep(1,N*T),x2,z2,x2*z2,x2*x2,z2*z2,s2,s2*x2,s2*z2,s2*x2*z2,s2*x2*x2,s2*z2*z2) 
# xx=cbind(xx,matrix(rep(t2,12),ncol=12)*xx,matrix(rep(t2*t2,12),ncol=12)*xx) # corresponds the interaction W1t and W2t

xx=cbind(rep(1,N*T),x2,x2*x2,s2,s2*x2,s2*x2*x2) 
xx=cbind(xx,matrix(rep(t2,6),ncol=6)*xx,matrix(rep(t2*t2,6),ncol=6)*xx) # corresponds the interaction W1t and W2t

#estimating reduced form logit
if(!hetero){
  #   b1=rep(0,ncol(xx))
  #   b1 = optim(b1,wlogitd,Y=(y2==0),X=xx,P=rep(1,N*T),method="BFGS")$par #default fnscale = 1 = min which is what we want
  b1=glm(((y2==0)*1)~xx-1,family='binomial')$coef #MUCH faster, same result
} else {
  # PType=.5*rep(1,2*N*T)
  # oPType=rep(0,2*N*T)
  # Pi2=c(.5,.5) 
  
  PType=.35*rep(1,2*N*T)
  oPType=rep(0,2*N*T)
  Pi2=c(.35,.65)
  
  b1 = rep(0,ncol(xx))
  #b1 = optim(b1,wlogitd,Y=(y2==0),X=xx,P=PType,method="BFGS")$par
  b1 = optim(b1,wlogitd,Y=(y2==0),X=xx,P=PType,method="Nelder-Mead",hessian=FALSE)$par # when BFGS does not work
  
  #For a binomial GLM prior weights are used to give the number of trials when the response is the proportion of successes
  #So we cannot send them into the "weights" argument
}

#calculating fv terms
if(!hetero){
  fvt1 = fvdataRcpp_Rust(b1,RX1,tbin,xbin,Zstate,Xstate,xtran,N,T,State)  
} else {
  fvt1 = fvdataRcpp_Rust(b1,RX1,tbin,xbin,Zstate,Xstate,xtran,N,T,rep(1,N),hetero)  
}

#estimating the structural parameters
#xccp = cbind(rep(1,N*T),x2*10,s2)
xccp = cbind(rep(1,N*T),x2*10,s2) # we don't multiply x2 by 10

if(!hetero){
  #bccp = alphac  
  #bccp = optim(bccp,wlogit,Y=y2,X=cbind(xccp,fvt1),P=rep(1,N*T),method="BFGS")$par
  bccp = glm(y2~cbind(xccp,fvt1)-1,family='binomial')$coef #much faster than 'optim'
} else {
  
  #starting the EM algorithm
  
  j=0
  
  bccp = alphac
  
  #intcondX=cbind(rep(1,N), X[1:N,1],Z[1:N,1])
  intcondX=cbind(rep(1,N), X[1:N,1],Z[1:N]) # value of regressors at t=1
  binit=rep(0,3)
  
  cond=0 
  lp=NULL
  while(cond==0){
    
    #updating PType
    ##first getting the type-specific likelihoods
    
    oPType=PType
    
    #replaces call to "likeCPP"
    U1 = cbind(xccp,fvt1)%*%bccp # v2-v1
    Like = (y2*exp(U1)+(1-y2))/(1+exp(U1)) # N*T*2 vector
    
    Like2=array(Like,c(N,T,2)) # rewrite Like in a cube
    base=apply(Like2,c(1,3),prod) # Likelihood by bus and by brand
    
    #now getting the initial condition parameters 
    intcond_optim=optim(binit,intcond,like=base,X=intcondX,method="BFGS")
    binit = intcond_optim$par # initial of delta -> used to compute pi_s(m) in intcondP
    lp=c(lp,intcond_optim$value) # value of the likelihood
    
    #and the PType's
    PType=intcondP(binit,base,intcondX) # q_ns :86 bus *2 states
    PType=kronecker(rep(1,T),PType) # q_nst : 86 bus * 2 states * time
    PType=as.vector(PType) # q_nst as vector
    
    #estimating reduced form logit(to find U1 as a function of 18 regressors)
    #b1 = optim(b1,wlogitd,Y=(y2==0),X=xx,P=PType,method="BFGS")$par   
    b1 = optim(b1,wlogitd,Y=(y2==0),X=xx,P=PType,method="Nelder-Mead",hessian=FALSE)$par # y2==0 <=> d1t=1
    
    
    #calculating fv terms
    fvt1=fvdataRcpp_Rust(b1,RX1,tbin,xbin,Zstate,Xstate,xtran,N,T,rep(1,N),hetero)
    
    
    # Structural paremeters
    bccp = optim(bccp,wlogit,Y=y2,X=cbind(xccp,fvt1),P=PType,method="BFGS")$par 
    
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


if(hetero){
  Iccp=c(Iccp,j)
  Binit=rbind(Binit,binit)  
}
  

bccp
