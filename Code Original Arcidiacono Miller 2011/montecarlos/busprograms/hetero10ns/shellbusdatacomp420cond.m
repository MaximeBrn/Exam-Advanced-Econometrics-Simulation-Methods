%shell for bus program

clear all;

Bccp=[];
Tccp=[];
Iccp=[];
Binit=[];
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
alphaf=[2.2577;-.1384;.4;2.656;0;-.1;0];
alphac=[2.233;-.1339;.4;.9115];
alphac2=[0;-.1;0];

MC=1;
for MC=1:50;
    
    %generating the data
    
    [Y,X,Z,Xstate,Zstate,State,TFV,adj]=genbus4(alpha,N,T,xtran,xtranc,xbin,zbin,xval,zval,T2);
    
    Adj=[Adj adj];
    y2=reshape(Y,N*T2,1);
    y2=[y2;y2];
    x2=reshape(X(:,1:T2),N*T2,1)./10;
    x2=[x2;x2];
    z2=kron(ones(T2,1),Z);
    z2=[z2;z2];
    s2=[zeros(N*T2,1);ones(N*T2,1)];
    t2=kron([1:T2]',ones(N,1))./10;
    t2=[t2;t2];
    stemp=[zeros(N,1);ones(N,1)];
    td=zeros(2*N*T2,T2-1);
    t=1;
    while t<T2
        
        td(:,t)=t2*10-1==t;
        t=t+1;
    end
    
    %estimating with data ccps
        
    tic

    %setting up data for reduced form logit
    
    xx=[ones(2*N*T2,1) x2 z2 x2.*z2 x2.*x2 z2.*z2 s2 s2.*x2 s2.*z2 s2.*x2.*z2 s2.*x2.*x2 s2.*z2.*z2];
    xx=[xx td];

    PType=.5*ones(2*N*T2,1);
    oPType=zeros(2*N*T2,1);
    Pi2=[.5 .5];
    
    %estimating reduced form logit
 
    b1=zeros(size(xx,2),1);

    [b1]=fminunc('wlogitd',b1,o2,y2==0,xx,PType);

    %calculating fv terms
    
    [fvt1]=fvdata(b1,RX1,tbin,xbin,Zstate,Xstate,xtran,N,T2);
    
    %starting the EM algorithm
    
    j=0;
   
    bccp=[alphac(1:4);zeros(T2-2,1)];
           
    xccp=[ones(N*T2*2,1) x2*10 s2];
    intcondX=[ones(N,1) X(1:N,1) Z(1:N,1)];
    binits=zeros(3,1);
    binit=binits;
    index=t2<1;
    lp=[];
    cond=0;
    tol=.0000001;
    while cond==0;
    
        %updating PType
        %first getting the type-specific likelihoods
        
        oPType=PType;
        Like=likeCCP(bccp,y2(index==1), [xccp(index==1,:) fvt1(index==1) td(index==1,1:T2-2)]);
        Like2=reshape(Like,N,T2-1,2);
        base=squeeze(prod(Like2,2));
        
        %now getting the initial condition parameters 
        if j>1
        if j<50
        [binit]=fminunc('intcond',binits,o1,base,intcondX);
        end
        if j>49
            [binit]=fminunc('intcond',binit,o1,base,intcondX);
        end
        end
        %and the PType's
        
        templ=intcond(binit,base,intcondX);
        
        lp=[lp;templ];
        
        [PType]=intcondP(binit,base,intcondX);
        
        PType=kron(ones(T2,1),PType);

        PType=reshape(PType,(N*T2*2),1);
        
        %estimating reduced form logit
            
        [b1]=fminunc('wlogitd',b1,o2,y2==0,xx,PType);
        
        %calculating fv terms
        
        [fvt1]=fvdata(b1,RX1,tbin,xbin,Zstate,Xstate,xtran,N,T2);
        
        
        [bccp]=fminunc('wlogit',[bccp],o1,y2(index==1),[xccp(index==1,:) fvt1(index==1) td(index==1,1:T2-2)],PType(index==1));
                
        if j>26;
            junk=abs((lp(j)-lp(j-25))./lp(j))<tol;
            junk2=abs((lp(j-1)-lp(j-26))./lp(j-1))<tol;
          
            cond=junk2.*junk;
        end
        j=j+1;
        
    end
    j
    tccp=toc
    
    Tccp=[Tccp;tccp];
    Bccp=[Bccp;[bccp]'];
    Iccp=[Iccp;j];
    Binit=[Binit;binit'];
    
  
    MC 
    
    save busdata0210 Bccp Tccp Iccp Binit Adj
end


