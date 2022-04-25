wlogit_nobeta=function(b,Y,X,P,Beta){
  
  U1=X[,1:3]%*%b+X[,4]*Beta
  Like=t(P)%*%(log(1+exp(U1))-Y*U1)
  
  Like
}

