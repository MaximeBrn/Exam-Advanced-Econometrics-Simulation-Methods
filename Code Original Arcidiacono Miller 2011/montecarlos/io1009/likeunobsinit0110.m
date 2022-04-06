

function [fullike]=typeprob0806(b,N,S,X,Like,nx);

num=ones(N,S);

s=2;

while s<S+1

    num(:,s)=exp(X*b((s-2)*nx+1:(s-1)*nx));
    
    s=s+1;
    
end

prior=num./(sum(num,2)*ones(1,S));


fullike=-sum(log(sum(prior.*Like,2)));
