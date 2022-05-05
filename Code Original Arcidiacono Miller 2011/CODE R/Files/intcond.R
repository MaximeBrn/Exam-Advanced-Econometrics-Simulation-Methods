# Intcond is a flexible function of the first-period observatables
# See section B.1.4 of the supplemental material (pdf p.13)
intcond=function(b,like,X){
  # x = regressors at time t=1 (W0 in the supplemental)
  # like = base = Likelihood by bus and by brand
  
  
  U1=X%*%b # W0*delta : cost of keeping at t=1 
  p=exp(U1)/(1+exp(U1)) # pi(2|x1) 
  p=cbind(p,1-p) # Used to weight the column of brand contribution to the likelihood
  
  
  Like = -sum(log(rowSums(p*like))) # objective function to find delta(m+1)
  
  Like
}