# This function is used to update the q_ns
# See equation 2.17 in the paper
intcondP=function(b,like,X){
  # like = base = Likelihood by bus and by brand
  # x = regressors at time t=1
  
  
  U1=X%*%b 
  p=exp(U1)/(1+exp(U1)) # we use the (m+1)-th delta to update the m-th pi(s|x) (see supplemental B.1.4)
  p=cbind(p,1-p) # vectorize the use the (m+1)-th pi(s|x) as weights
  
  PType = (like*p)/(rowSums(p*like)%*%matrix(1,1,2)) # Update the m-th q_ns (see equation 2.17)
  
  PType
}