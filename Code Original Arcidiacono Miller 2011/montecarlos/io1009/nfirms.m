%transition probs on the number of firms

function BigP=nfirms(pe1,pi1,bine,bini,ne,N)

%first getting the probability associated with each possible non-incumbent
%outcome

%first row gives the probability that all the non-incumbents stay out
%next row gives the probability that all the non-incumbents but one stay out,
%and so on

Pe=zeros(ne+1,1);

j=0;

while j<ne+1

    Pe(j+1)=bine(j+1)*(pe1^j)*((1-pe1)^(ne-j));
    
    j=j+1;
end

%now getting the probability associated with each possible incumbent
%outcome

%first row gives the probability that all the incumbents exit
%next row gives the probability that all the incumbents but one exit,
%and so on

Pi=zeros(N-ne+1,1);

k=0;

while k<N-ne+1

    Pi(k+1)=bini(k+1)*(pi1^k)*((1-pi1)^(N-ne-k));
    
    k=k+1;
end

%Now aggreating Pe and Pi to form the probability of different numbers of
%firms in each period.

BigP=zeros(1,N+1);

j=0;

while j<ne+1
    
    k=0;
    
    while k<N-ne+1
        
        %j+k gives different ways of having the same number of incumbent
        %firms
        
        BigP(j+k+1)=Pe(j+1)*Pi(k+1)+BigP(j+k+1);
        
        k=k+1;
    end
    
    j=j+1;
    
end




