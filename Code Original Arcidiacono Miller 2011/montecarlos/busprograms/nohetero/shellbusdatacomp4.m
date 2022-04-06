%shell for bus program


clear all;

Bccp=[];
Tccp=[];
Bf1=[];
Tf1=[];

%parameter values
alpha=[2;-.15;1;.9;.4];
Beta=alpha(4);
eul=0.577215665;

T=20;
N=1000;

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

%monte carlos

%starting values for FIML and CCP
alphaf=[alpha(1:3);log(alpha(4))-log(1-alpha(4))];
alphac=[alpha(1:4)];
alphac2=[0;-.1;0];

MC=1;
for MC=1:50;
    
    %generating the data
    
    [Y,X,Z,Xstate,Zstate,State,TFV]=genbus4(alpha,N,T,xtran,xtranc,xbin,zbin,xval,zval);
    y2=reshape(Y,N*T,1);
    x2=reshape(X(:,1:T),N*T,1)./10;
    z2=kron(ones(T,1),Z);
    s2=kron(ones(T,1),State);
    t2=kron([1:T]',ones(N,1))./10;
            
    %estimating FIML
    tic
    %bf=0;
    [bf]=fminunc('likebusML4',alphaf,o1,[Y],s2,N,T,[X],[Zstate],[Xstate],xtran,tbin,zbin,xbin,xval,Z);
    tf1=toc
    
    Tf1=[Tf1;toc];
    Bf1=[Bf1;bf']; 
    %estimating with data ccps
        
    tic

    %setting up data for reduced form logit
    
    xx=[ones(N*T,1) x2 z2 x2.*z2 x2.*x2 z2.*z2 s2 s2.*x2 s2.*z2 s2.*x2.*z2 s2.*x2.*x2 s2.*z2.*z2];
    xx=[xx ((t2)*ones(1,12)).*xx ((t2.*t2)*ones(1,12)).*xx];

   
    %estimating reduced form logit
 
    b1=zeros(size(xx,2),1);

    [b1]=fminunc('wlogitd',b1,o2,y2==0,xx,ones(N*T,1));

    %calculating fv terms
    
    [fvt1]=fvdata(b1,RX1,tbin,xbin,Zstate,Xstate,xtran,N,T,State);
    
    %estimating the structural parameters
   
    bccp=alphac;
    xccp=[ones(N*T,1) x2*10 s2];
    
    [bccp]=fminunc('wlogit',[bccp],o1,y2,[xccp fvt1],ones(N*T,1));

    tccp=toc
    
    Tccp=[Tccp;tccp];
    Bccp=[Bccp;bccp'];
  
    MC 
    
    save busdata Bccp Tccp Bf1 Tf1
end


