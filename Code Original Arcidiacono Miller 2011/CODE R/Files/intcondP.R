intcondP=function(b,like,X){
  # like = base = Likelihood by bus and by brand
  # x = regressors at time t=1
  
  
  U1=X%*%b 
  p=exp(U1)/(1+exp(U1)) # pi_s(m)
  p=cbind(p,1-p) 
  
  PType = (like*p)/(rowSums(p*like)%*%matrix(1,1,2)) # eq.(2.17)
  
  PType
}