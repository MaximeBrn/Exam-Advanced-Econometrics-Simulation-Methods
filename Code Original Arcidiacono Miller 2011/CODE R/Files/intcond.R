intcond=function(b,like,X){
  # x = regressors at time t=1
  # like = base = Likelihood by bus and by brand
  
  U1=X%*%b # W0*delta : cost of keeping at t=1 (supp p.13)
  p=exp(U1)/(1+exp(U1)) # pi(2|x1) (supp p.13)
  p=cbind(p,1-p) # to weight the column of like
  
  
  Like = -sum(log(rowSums(p*like))) # objective of delta(m+1) (supp p. 13)
  
  Like
}