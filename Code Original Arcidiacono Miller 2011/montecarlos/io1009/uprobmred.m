
function [p]=uprobm1109(b,N,S,Xn)


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
                    
                    v=[1 xn-1 ((xn-1)./10).^2 i n1 (n1./5).^2 s-1 (s-1).*(xn-1)./10 i.*(s-1) n1.*(s-1)./10]*b;
                
                    p(n1+1,xn,s,i+1)=exp(v)./(1+exp(v));
                    
                    xn=xn+1;
                end
                
                s=s+1;
                
            end
            
            n1=n1+1;
        end
        
        i=i+1;
        
    end

    
end
