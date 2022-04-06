%shell for bus program
%throwing out first and last ten years of data

clear all;

Bccp=[];
Tccp=[];
Adj=[];
%parameter values
alpha=[2;-.15;1;.9;.4];
Beta=alpha(4);
eul=0.577215665;

T=20;
T2=10;
N=2000;

%optimization choices
%derivative is used in the minimization only for reduced form ccp's

o1=optimset('LargeScale','off','Display','off');
o2=optimset('LargeScale','off','Display','off','GradObj','on');

%Create transition matrices;

zval=(.25:.01:1.25)';
zbin=length(zval);
xval=(0:.125:25)';
xbin=length(xval);
for z=1:zbin;
    [xtran(1+(z-1)*xbin:z*xbin,:),xtranc(:,:,z)]=xgrid(zval(z),xval);
end;

tbin=xbin*zbin;

%z and x values for each state

zvalr=kron(zval,ones(xbin,1));
xvalr=kron(ones(zbin,1),xval)./10;

%data for reduced form logits

RX1=[ones(zbin*xbin,1) xvalr zvalr xvalr.*zvalr xvalr.*xvalr zvalr.*zvalr];

RX1f=[ones(zbin*xbin,1) xvalr zvalr];

%monte carlos

%starting values for FIML and CCP
alphaf=[alpha(1:3);log(alpha(4))-log(1-alpha(4));0;-.1;0];
alphac=[alpha(1:4)];
alphac2=[0;-.1;0];

MC=1;
for MC=1:50;
    
    %generating the data
    
    [Y,X,Z,Xstate,Zstate,State,TFV,adj]=genbus4(alpha,N,T,xtran,xtranc,xbin,zbin,xval,zval,T2);
    Adj=[Adj adj];
    y2=reshape(Y,N*T2,1);
    x2=reshape(X(:,1:T2),N*T2,1)./10;
    z2=kron(ones(T2,1),Z);
    s2=kron(ones(T2,1),State);
    t2=kron([1:T2]',ones(N,1))./10;

    td=zeros(N*T2,T2-1);
    t=1;
    while t<T2
        
        td(:,t)=t2*10-1==t;
        t=t+1;
    end
    
    %estimating with data ccps
        
    tic

    %setting up data for reduced form logit
    
    xx=[ones(N*T2,1) x2 z2 x2.*z2 x2.*x2 z2.*z2 s2 s2.*x2 s2.*z2 s2.*x2.*z2 s2.*x2.*x2 s2.*z2.*z2];
    xx=[xx td];

   
    %estimating reduced form logit
 
    b1=zeros(size(xx,2),1);
    
    PType=ones(N*T2,1);
    [b1]=fminunc('wlogitd',b1,o2,y2==0,xx,PType);

    %calculating fv terms
    
    [fvt1]=fvdata(b1,RX1,tbin,xbin,Zstate,Xstate,xtran,N,T2,State);
    
    %estimating
   
    bccp=[alphac(1:4);zeros(T2-2,1)];
           
    xccp=[ones(N*T2,1) x2*10 s2];

    index=t2<1;
      
    [bccp]=fminunc('wlogit',[bccp],o1,y2(index==1),[xccp(index==1,:) fvt1(index==1) td(index==1,1:T2-2)],PType(index==1));
    
    bccp(4)
    tccp=toc
    
    Tccp=[Tccp;tccp];
    Bccp=[Bccp;[bccp]'];
 
  
    MC 
    
    save busdata2 Bccp Tccp Adj
end


