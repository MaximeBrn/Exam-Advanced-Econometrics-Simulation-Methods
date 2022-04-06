%this program is a shell program that allows for differences in intercept
%terms over time which are unobserved to the econmetrician.  The first set
%of lines allow the user to change
%1) the number of observations per time period,
%2) the number of other potential firms, 
%3) the number of unobserved states, 
%4) the coefficients on the unobserved states,
%5) the transitions on the unobserved states, and
%7) the number of times the user wants to run the simulation
%throughout the variance on the epsilon's is set to 1
clear

%SETTING PARAMETERS

Nsim=100;
Nm=3000;
Nf=5;
T=10;
Tl=10;
S=5;
Xn=10;
adj=1;
Beta=.9;

%SETTING COEFFICIENTS TO BE ESTIMATED
%structural entry coefficients
bx=-.05;
bx0=0;
bnf=-.2;
bss=.25;
be=-1.5;
%price coefficients (intercept, number of firms, x, state)

bp=[7; -.4 ;-.1 ;.3];

%SETTING PARAMETERS GOVERNING TRANSITIONS OF THE UNOBSERVABLES
ps=.7;
nps=(1-ps)/(S-1);
trans=ps*eye(S)+nps*(1-eye(S));

%SETTING THE STARTING VALUES FOR THE PARAMETERS

bpriors=[zeros(S-1,1);.9];
bs=[bx0+.2;bx;be;bnf];
bps=[bp(1)+.3;bp(2:3)];
sigs=1;
transs=trans;

%SPECIFYING THE CONVERGENCE CRITERIA
%TOL IS HOW SMALL THE DIFFERENCE IN LIKELIHOODS BETWEEN 25 ITERATIONS
%DIVIDED BY THE LAST LIKELIHOOD MUST BE TWO CONSECUTIVE TIMES FOR THE
%EM ALGORITHM TO CONVERGE

tol=.0000001;

%DECLARING VARIABLES, SETTING PATHS AND OPTIMIZATION PREFERENCES

addpath '/afs/econ.duke.edu/home/p/psarcidi/Papers Current/ccp/io1009'

o1=optimset('Display','off','LargeScale','off','GradObj','on');

o2=optimset('Display','off','LargeScale','off');

Bccp=[];
Byccp=[];
Tccp=[];
Iccp=[];
Lccp=[];
Binit=[];
Byfe=[];

%CREATING THE CDF OF THE UNOBSERVED STATE

ctrans=trans(:,1);

for s=1:S-2
    ctrans(:,s+1)=ctrans(:,s)+trans(:,s+1);
end

nf=[0:Nf]';

%getting the utilities for every state combination

Util=zeros(Nf+1,Xn,S,2);

j=0;

while j<Nf+1
    
    s=0;
    
    while s<S
        
        x=0;
        
        while x<Xn
        
            i=0;
        
            while i<2
            
                Util(j+1,x+1,s+1,i+1)=bx0+bx*(x/adj)+bnf*j+bss*s+be*(1-i);
            
                i=i+1;
            
            end
            x=x+1;
        end
        s=s+1;
    end
    j=j+1;
end

%getting the binomial coefficients needed to combine the individual choice
%probabilities to probabilities over numbers of firms

bine=zeros(Nf+1,Nf+1);
bine(1,:)=1;

n=1;

while n<Nf+1
    
    k=0;
    
    while k<n+1
        
        bine(n+1,k+1)=factorial(n)/(factorial(k)*factorial(n-k));
        
        k=k+1;
        
    end
    
    n=n+1;
    
end

%FINDING THE PROBABILITIES OF ENTERING (OR STAYING IN) GIVEN ALL POSSIBLE
%STATES. THIS IS USED FOR DATA CREATION

[p2]=prob1009(Util,trans,Nf,Xn,S,bine,Beta);

%RUNNING THE ESTIMATION

ns=1;

while ns<Nsim+1
    
    %calling data creation program
    %output is nmarketXntimeXnfirm number for entry/exit decisions of the firm
    %X is then nmarket
    %State is then nmarketXntime
    
    [Firm,X,State,Y,LFirm]=createdata1109(p2,ctrans,S,Xn,Nf+1,Nm,T,bp,Tl);
    
    %RESHAPING THE DATA 
    
    S2=reshape(State(:,1:T),Nm*T,1);
    
    NFirm=sum(Firm,3);
    LNFirm=sum(LFirm,3);
    
    NFirm2=kron(ones(Nf+1,1),reshape(NFirm,Nm*T,1));
    LNFirm2=kron(ones(Nf+1,1),reshape(LNFirm,Nm*T,1));
    
    X=kron(ones(T,1),X);
   
    Z2=kron(ones((Nf+1),1),[ones(Nm*T,1) X-1]);

    Y2=reshape(Y,Nm*T,1);
    
    ZZ2=[ones(Nm*T,1) reshape(NFirm,Nm*T,1) X-1];
     
    Firm2=[];
    LFirm2=[];
    for n=1:Nf+1;
        Firm2=[Firm2;reshape(Firm(:,:,n),Nm*T,1)];
        LFirm2=[LFirm2;reshape(LFirm(:,:,n),Nm*T,1)];
    end
    
    LNFirm2=LNFirm2-LFirm2;
    
    Z2=[Z2 1-LFirm2 LNFirm2];
          
    %now ready to estimate via CCP's
    
    %TAKES AS GIVEN THE TRANSITION PROBABILITIES ON THE UNOBSERVED STATES
    tic
    
    bprior=bpriors;
    b2=bs;
    sig=sigs;
    bp2=bps;
    trans2=transs;
    
    %GETTING PARAMETERS FOR REDUCED FORM CCP's
    Zred=[ones(Nm*T*(Nf+1),1) Z2(:,2) (Z2(:,2)./10).^2 LFirm2 LNFirm2 (LNFirm2./5).^2];

    [btemp]=fminunc('wlogitd',zeros(size(Zred,2),1),o1,Firm2,Zred,0,ones(Nm*(Nf+1)*T,1));
  
    btemp=[btemp;zeros(4,1)];
    %GETTING THE FUTURE VALUE TERMS
    
    p=uprobmred(btemp,Nf,S,Xn);
        
    [fv,bigp]=fv1109(p,trans2,Nf,Xn,S,bine);
 
    %MATCHING TO THE DATA AND GETTING THE FUTURE VALUE TERMS
    
    [Fv2,BigP2]=match1109(LNFirm2,LFirm2,Z2(:,2)+1,Z2(:,3)+1,fv,bigp,Nm*T*(Nf+1),Nf+1);

    Z2(:,4)=BigP2*nf;
        
    %ESTIMATING THE STRUCTURAL COEFFICIENTS 
        
    [b2]=fminunc('wlogitdadj',b2,o1,Firm2,Z2,Fv2,ones(Nm*T*(Nf+1),1));
  
    [bp2,sig]=wols(Y2,ZZ2,ones(Nm*T,1),Nm*T);
  
    t1=toc
    ns=ns+1
    
    Yhat=reshape(Y-mean(Y,2)*ones(1,10),Nm*T,1);
    [bpfe,sigfe]=wols(Yhat,ZZ2(:,1:2),ones(Nm*T,1),Nm*T);
    
    %SAVING THE PARAMETERS
    Bccp=[Bccp; b2'];
    Byccp=[Byccp;[bp2;sig]'];
    Byfe=[Byfe;[bpfe;sigfe]'];
    Tccp=[Tccp;t1];
    
    save bwrong Bccp Byccp Tccp Byfe
    
end
        
        