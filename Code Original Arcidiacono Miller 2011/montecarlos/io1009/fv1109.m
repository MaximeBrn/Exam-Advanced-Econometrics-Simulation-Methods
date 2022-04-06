
function [fv,BigP]=fv1109(p,trans,N,Xn,S,bine)

eul=.5772;
beta=.9;

fv=zeros(N+1,Xn,S,2);
BigP=zeros(N+1,N+1,Xn,S,2);

%calculating the probability of making various choices given the current
%state
%first row is the case when all other competitors are entrants
%second row has all other competitors as entrants but one, and so on 

   
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
                        
                        fv(n1+1,xn,s,i+1)=fv(n1+1,xn,s,i+1)-trans(s,s2)*(BigP(n1+1,:,xn,s,i+1)*log(1-p(:,xn,s2,2)));
                        
                        s2=s2+1;
                        
                    end

                    xn=xn+1;
                end
                
                s=s+1;
                
            end
            
            n1=n1+1;
        end
        
        i=i+1;
        
    end
fv=.9*(fv+.5772);
end

 
