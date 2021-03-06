// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List genbusRcpp(vec const& alpha,int N,int T,mat const& xtran,mat const& xtrancRcpp,
                int xbin,int zbin,vec const& xval,vec const& zval){
  
  double Beta = alpha[3];   //discount factor
  double Pi = alpha[4];     
  double eul = 0.577215665; // Euler constant (gamma in the paper)
  cube FV = zeros<cube>(xbin*zbin,2,T+10+1); // matrix of discounter value function with terminal value 0
  
  //Since Rcpp doesn't allow cubes to be sent in, we recreate the cube
  cube xtranc = zeros<cube>(xbin,xbin,zbin); 
  int inc;
  for(int z=0;z<zbin;z++){
    inc = z*xbin;
    xtranc.slice(z) = xtrancRcpp(span(),span(inc,inc+xbin-1));
  }
  
  int adj, adj0, z2, z3;
  double util1,util0,dem,p0;
  vec FVtsz;
  List out;
  
  int t=T+10;
  while(t>1){
    for(int s=0;s<=1;s++){        //0:1
      for(int z=1;z<=zbin;z++){   //1:zbin
        for(int x=1;x<=xbin;x++){ //1:xbin
          adj  = x+(z-1)*xbin-1; // index of x1=x and x2=z in xtran
          FVtsz = vectorise(FV(span((z-1)*xbin,z*xbin-1),span(s),span(t))); // dicounted value function for this (all x1,x2,s,t)
          util1=alpha[0]+alpha[1]*xval[x-1]+alpha[2]*s+as_scalar(xtran(span(adj),span())*FVtsz);
          util0=as_scalar(xtran(span((z-1)*xbin),span())*FVtsz);
          FV(adj,s,t-1)=Beta*(log(exp(util1)+exp(util0)) + eul);
        }
      }  
    }
    t=t-1;
  }

  vec Statedraw = randu<vec>(N);
  vec State = zeros<vec>(N); 
  State.elem(find(Statedraw>Pi)).ones(); //gives the states of each bus (1 or 2)
  mat Y = zeros<mat>(N,T+10);
  mat X = zeros<mat>(N,T+1+10);
  mat Xstate = ones<mat>(N,T+1+10);
  vec Zstate=ceil(zval.size()*runif(N));
   
  mat Draw1 = randu<mat>(N,T+10); // To be compared to ccp of replacing
  mat Draw2 = randu<mat>(N,T+10); // To be compared to the xtranc

  vec Z = zeros<vec>(N);
  mat FVT = zeros<mat>(N,T+1+10);
  
  for(int n = 0; n<N; n++){
    
    Z[n] = zval[Zstate[n]-1];
    adj0 = (Zstate[n]-1)*xbin;
    z2 = (Zstate[n]-1)*xbin;
    z3 = z2+xbin-1;
    
    for(int t = 0; t<(T+10);t++){
      
      adj = Xstate(n,t)+(Zstate[n]-1)*xbin-1;
      vec FVnt = vectorise(FV(span(z2,z3),span(State[n]),span(t+1)));
      FVT(n,t)=as_scalar(xtran(span(adj),span())*FVnt-xtran(span(adj0),span())*FVnt);
      util1 = alpha[0]+alpha[1]*X(n,t)+alpha[2]*State[n]+as_scalar(xtran(span(adj),span())*FVnt);
      util0 = as_scalar(xtran(span(adj0),span())*FVnt);
      dem=exp(util1)+exp(util0); // see eq. (2.3)
      p0=exp(util0)/dem; // see eq. (2.3)
      Y(n,t) = 1-(Draw1(n,t)<p0); // If draw < p0 => Y=0 => replace
      vec xbin1 = ones<vec>(xbin-1);
      Xstate(n,t+1)=1+(Y(n,t)==1)*sum((Draw2(n,t)*xbin1)>vectorise(xtranc(span(Xstate(n,t)-1),span(0,(xbin-2)),span(Zstate[n]-1)))); // find the rank of x_t+1 in xval
      X(n,t+1)=xval[Xstate(n,t+1)-1]; // assign value to x_t+1
    }
  }

  Y=Y(span(),span(10,T+10-1));
  X=X(span(),span(10,T+10-1));
  Xstate=Xstate(span(),span(10,T+10-1));
  FVT=FVT(span(),span(10,T+10-1));

  return List::create(
    Named("Y") = Y, // Return the decision (Y=0 if we replace)
    Named("X") = X, // X1_t
    Named("Z") = Z, // X2
    Named("Xstate") = Xstate,// Index of X in xval
    Named("Zstate") = Zstate, // Index of Z in zval
    Named("State") = State, // brand of the bus
    Named("FVT") = FVT);
}