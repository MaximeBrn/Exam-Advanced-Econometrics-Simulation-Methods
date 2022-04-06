function [Like]=likebusML4(alpha,Y,State,N,T,X,Zstate,Xstate,xtran,tbin,zbin,xbin,xval,Z)

Beta=exp(alpha(4))/(1+exp(alpha(4)));

FV=zeros(tbin,2,T+1);

t=T;
while t>1
    for s=0:1;
        for z=1:zbin;
            for x=1:xbin
              adj=x+(z-1)*xbin;  
              util1=alpha(1)+alpha(2)*xval(x)+alpha(3)*s+ xtran(adj,:)*FV((z-1)*xbin+1:z*xbin,s+1,t+1);
              util0=xtran(1+(z-1)*xbin,:)*FV((z-1)*xbin+1:z*xbin,s+1,t+1);
              FV(adj,s+1,t)=Beta*log(exp(util1)+exp(util0));
              
            end;   
        end;
    end;
    t=t-1;
end; 

Like=0;

for n=1:N;
    adj0=(Zstate(n)-1)*xbin+1;
    z2=(Zstate(n)-1)*xbin+1;
    z3=z2+xbin-1;
    for t=1:T;
        adj=Xstate(n,t)+(Zstate(n)-1)*xbin;
        util1=alpha(1)+alpha(2)*X(n,t)+alpha(3)*State(n)+((xtran(adj,:)-xtran(adj0,:))*FV(z2:z3,State(n)+1,t+1));
        dem=exp(util1)+1;
        Like=Like+log(dem)-((Y(n,t)==1).*util1);
    end;
end;








