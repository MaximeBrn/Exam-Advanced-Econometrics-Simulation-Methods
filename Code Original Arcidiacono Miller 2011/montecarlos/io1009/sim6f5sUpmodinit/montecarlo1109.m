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

bpriors=[zeros(2*(S-1),1);.9];
bs=[bx0+.2;bx;bss-.1;be;bnf];
ub=10*ones(size(bpriors,1),1);
lb=-ub;
bps=[bp(1)+.3;bp(2:3);bp(4)-.1];
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
        
    btemp=regress(LNFirm(:,1),[ones(Nm,1) X X.^2]);
    Resid=LNFirm(:,1)-([ones(Nm,1) X X.^2])*btemp;
    
    NFirm2=kron(ones((Nf+1)*S,1),reshape(NFirm,Nm*T,1));
    LNFirm2=kron(ones((Nf+1)*S,1),reshape(LNFirm,Nm*T,1));
    
    
    X=kron(ones(T,1),X);
   
    Z2=kron(ones(S*(Nf+1),1),[ones(Nm*T,1) X-1]);
    
    temp=kron([0:S-1]',ones(Nm*T,1));
    
    Z2=[Z2 kron(ones(Nf+1,1),temp)];
    
    Y2=kron(ones(S,1),reshape(Y,Nm*T,1));
    
    ZZ2=[kron(ones(S,1),[ones(Nm*T,1) reshape(NFirm,Nm*T,1) X-1]) temp];
     
    Firm2=[];
    LFirm2=[];
    for n=1:Nf+1;
        Firm2=[Firm2;kron(ones(S,1),reshape(Firm(:,:,n),Nm*T,1))];
        LFirm2=[LFirm2;kron(ones(S,1),reshape(LFirm(:,:,n),Nm*T,1))];
    end
    
    LNFirm2=LNFirm2-LFirm2;
    
    Z2=[Z2 1-LFirm2 LNFirm2];
          
    %now ready for the EM algorithm
    
    %FIRST SETTING THE INITIAL GUESS ON THE PRIOR PROBABILITY OF BEING IN
    %EACH STATE AS WELL AS THE INITIAL GUESS ON THE COEFFS
    tic
    
    bprior=bpriors;
    b2=bs;
    sig=sigs;
    bp2=bps;
    trans2=transs;
    
    %GETTING THE INITIAL CCPS AND FUTURE VALUE TERMS
    [btemp]=fminunc('wlogitd',zeros(4,1),o1,Firm2,[ones(Nm*T*S*(Nf+1),1) Z2(:,2) LFirm2 LNFirm2],0,ones(Nm*(Nf+1)*T*S,1));
    btemp=[btemp;0];
    p=uprobmred010410(btemp,Nf,S,Xn);
    
    [fv,bigp]=fv1109(p,trans,Nf,Xn,S,bine);

    %MATHING THE FUTURE VALUE TERMS TO THE DATA
    
    [Fv2,BigP2]=match1109(LNFirm2,LFirm2,Z2(:,2)+1,Z2(:,3)+1,fv,bigp,Nm*T*S*(Nf+1),Nf+1);

    Zinit=[ones(Nm,1) Resid];
    
    ninit=2;
    
    cond=0;
    lp=[0];
    j=1;     
    
    o3=optimset('Display','off','LargeScale','off','MaxIter',1);

    while cond==0; 
       
        %CALCULATING THE LIKELIHOOD GIVEN THE PARAMETER GUESSES
        
        Like=likecalc1109(bp2,sig,b2,Y2,ZZ2,Z2,Fv2,Firm2,Nm*S*T,Nf+1);
        Like=reshape(Like,Nm,T,S)*T*S;

        %ESTIMATING THE INITIAL CONDITION PARAMETERS
        if j>1
            if j==11;
                o3=optimset('Display','off','LargeScale','off');
            end
        [bprior,lpx]=fmincon('likeunobstrans1209',bprior,[],[],[],[],lb,ub,[],o3,Nm,T,S,Zinit,Like,ninit);
       
        
        %KEEPING TRACK OF THE LIKELIHOOD AT EACH ITERATION
        lp=[lp;lpx];
        end
        %CALCULATING THE CONDITIONAL PROBABILITY OF BEING A PARTICULAR TYPE
        %AS WELL AS UPDATING THE TRANSITIONS ON THE UNOBSERVABLES
        
        [PType,trans2]=typeprob1209(bprior,Zinit,ninit,Nm,T,S,Like);
        
        PType2=kron(ones(Nf+1,1),PType);
        
        %UPDATING THE CCP'S WITH THE MODEL
        
        k=1;
        
        while k<2;
               
            p=uprobm1109(b2,fv,bigp,Nf,S,Xn);
        
            [fv,bigp]=fv1109(p,trans2,Nf,Xn,S,bine);
 
            k=k+1;
        end
        
        %MATCHING TO THE DATA AND GETTING THE FUTURE VALUE TERMS
        [Fv2,BigP2]=match1109(LNFirm2,LFirm2,Z2(:,2)+1,Z2(:,3)+1,fv,bigp,Nm*T*S*(Nf+1),Nf+1);

        Z2(:,5)=BigP2*nf;
        
        %ESTIMATING THE COEFFICIENTS 
        
        [b2]=fminunc('wlogitdadj',b2,o1,Firm2,Z2,Fv2,PType2);

        [bp2,sig]=wols(Y2,ZZ2,PType,Nm*T);

        %CHECKING CONVERGENCE
        
            if j>26;
                junk=abs((lp(j)-lp(j-25))./lp(j))<tol;
                junk2=abs((lp(j-1)-lp(j-26))./lp(j-1))<tol;
          
                cond=junk2.*junk;
                if j>1000;
                    cond=1;
                end
  
            end

        j=j+1;
        
    end
    
        
    t1=toc
       Bccp=[Bccp; b2'];
    Byccp=[Byccp;[bp2;sig]'];
    Tccp=[Tccp;t1];
    Iccp=[Iccp;j];
    Lccp=[Lccp;lpx];
    Binit=[Binit;bprior'];
    
    j
    ns=ns+1
    
    %SAVING THE PARAMETERS
 
   
    save bmodel Bccp Byccp Tccp Iccp Lccp Binit
    
    end
        
        