%calculates fv terms given the data 

function fvt1=fvdata(b1,RX1,tbin,xbin,Zstate,Xstate,xtran,N,T);


    FV1=zeros(tbin,2,T+1);
    
    for t=2:T;
        
            p0=plogit(b1,kron([1 t/10 (t/10)*(t/10)],[RX1]));
            FV1(:,t)=(-log(p0));
        
    end;
   

    for n=1:N;
        
        adj0=(Zstate(n)-1)*xbin+1;
        z2=(Zstate(n)-1)*xbin+1;
        z3=z2+xbin-1;
        
        for t=1:T;
            adj=z2+Xstate(n,t)-1; 
                                  
            FVT1(n,t)=(xtran(adj,:)-xtran(adj0,:))*FV1(z2:z3,t+1);
    
        end;
        
    end;
    fvt1=reshape(FVT1,N*T,1);