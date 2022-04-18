xgrid=function(theta,xval){

  n=length(xval)
  xub=c(xval[2:n],Inf) #
  xtran=matrix(0,n,n)
  xtranc=matrix(0,n,n)
  lcdf=0
  for(i in 1:n){
    xtran[,i]=(xub[i]>=xval)*(1-exp(-theta*(xub[i]-xval))-lcdf) # see eq. (7.2)
    lcdf=xtran[,i]+lcdf
    xtranc[,i]=xtranc[,i]+lcdf # cumulated distribution function
  }

  list(xtran=xtran,xtranc=xtranc)
}

