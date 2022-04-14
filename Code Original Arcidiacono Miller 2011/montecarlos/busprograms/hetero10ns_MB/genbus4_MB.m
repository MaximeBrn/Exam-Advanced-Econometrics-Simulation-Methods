%bus data generation
%adds another variable
%change from orig folder to go with Jon's parameters for Z and X from 3_z

function [Y,X,Z,Xstate,Zstate,State,FVT,Adj]=genbus(alpha,N,T,xtran,xtranc,xbin,zbin,xval,zval,T2);

%% Describe the inputs
% alpha is the vector of true parameter values :
% alpha=[theta_0,theta_1,theta_2,beta,?]
% N is the number of buses
% T = 20
% xtran : transition matrix for mileages implied by equation (7.2)
% xtranc : 
% xbin : length of possible mileage
% zbin
% xval : possible values of x1t
% zval : support of x2
% T2 = 10

%% Describe the outputs


%% Code
Beta=alpha(4); % 
Pi=alpha(5);
eul=0.577215665;
FV=zeros(xbin*zbin,2,T+10+1);

Adj=zeros(T+11,1);

t=2;

while t<T+12
    
    Adj(t)=.7*Adj(t)+.5*randn(1);
    t=t+1;
end

t=T+10;
while t>1 % for each period
    for s=0:1; % for each brand
        for z=1:zbin; % for each value of x2
            for x=1:xbin % for each value of x1
              adj=x+(z-1)*xbin;  
              util1=alpha(1)+Adj(t)+alpha(2)*xval(x)+alpha(3)*s+ xtran(adj,:)*FV((z-1)*xbin+1:z*xbin,s+1,t+1);
              util0=xtran(1+(z-1)*xbin,:)*FV((z-1)*xbin+1:z*xbin,s+1,t+1);
              FV(adj,s+1,t)=Beta*(log(exp(util1)+exp(util0)) + eul);
              
            end;   
        end;
    end;
    t=t-1;
end; 

State=(rand(N,1)>Pi);
Y=zeros(N,T2+10);
Zstate=ceil(length(zval)*rand(N,1));
X=zeros(N,1);

Xstate=ones(N,T2+1+10);
Draw=rand(N,T2+10);
Draw2=rand(N,T2+10);
Draw3=rand(N,T2+10);
Xstate=X+1;
for n=1:N;
    
    Z(n,1)=zval(Zstate(n,1));
    adj0=(Zstate(n)-1)*xbin+1;
    z2=(Zstate(n)-1)*xbin+1;
    z3=z2+xbin-1;
    for t=1:T2+10;
        adj=Xstate(n,t)+(Zstate(n)-1)*xbin;
        FVT(n,t)=(xtran(adj,:)*FV(z2:z3,State(n)+1,t+1))-(xtran(adj0,:)*FV(z2:z3,State(n)+1,t+1));
        util1=alpha(1)+Adj(t)+alpha(2)*X(n,t)+alpha(3)*State(n) + (xtran(adj,:)*FV(z2:z3,State(n)+1,t+1));
        util0=(xtran(adj0,:)*FV(z2:z3,State(n)+1,t+1));
        dem=exp(util1)+exp(util0);
        p0=exp(util0)./dem;
        p1=1-p0;
        Y(n,t)=1-(Draw(n,t)<p0);
        Xstate(n,t+1)=1+(Y(n,t)==1).*sum((Draw2(n,t)*ones(1,xbin-1))>xtranc(Xstate(n,t),1:xbin-1,Zstate(n,1)),2);
        X(n,t+1)=xval(Xstate(n,t+1));
    end;
end;

Y=Y(:,11:T2+10);
X=X(:,11:T2+10);
Xstate=Xstate(:,11:T2+10);
FVT=FVT(:,11:T2+10);

        

    
