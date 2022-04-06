
function [p,BigP,fv]=prob1009(Util,trans,N,Xn,S,bine,Beta)
BigP=zeros(N+1,N+1,Xn,S,2);
fv=zeros(N+1,Xn,S,2);
eul=.5772;

p=exp(Util)./(1+exp(Util));

p0=zeros((N+1)*Xn*S*2,1);
p2=reshape(p,(N+1)*Xn*S*2,1);

%calculating the probability of making various choices given the current
%state
%first row is the case when all other competitors are entrants
%second row has all other competitors as entrants but one, and so on 
while max(abs(p0-p2))>.0000000001;
    

    p0=p2;
    
    %loop over whether you're an incumbent
    i=0;
    
    while i<2
        
        %loop over how many other firms were incumbents
        
        n1=0;
        
        while n1<N+1
            
            %loop over transitory states
            
            s=1;
            
            while s<S+1
                
                %loop over market characterisitics
                
                xn=1;
                
                while xn<Xn+1;
                    
                    %loop over future transitory state
                    
                    BigP(n1+1,:,xn,s,i+1)=nfirms(p(min(n1+1+i,N+1),xn,s,1),p(max(n1+i,1),xn,s,2),bine(N-n1+1,:),bine(n1+1,:),N-n1,N);
                    
                    v=0;
                    
                    s2=1;
                    
                    while s2<S+1
                        
                        v=v+trans(s,s2)*(BigP(n1+1,:,xn,s,i+1)*log(1-p(:,xn,s2,2)));
                        
                        s2=s2+1;
                        
                    end
                    fv(n1+1,xn,s,i+1)=-v;
                    tu=BigP(n1+1,:,xn,s,i+1)*Util(:,xn,s,i+1)-Beta*v+Beta*eul;
                    p(n1+1,xn,s,i+1)=exp(tu)./(1+exp(tu));

                    xn=xn+1;
                end
                
                s=s+1;
                
            end
            
            n1=n1+1;
        end
        
        i=i+1;
        
    end

    p2=reshape(p,(N+1)*Xn*S*2,1);

end
fv=Beta*(fv+eul);

