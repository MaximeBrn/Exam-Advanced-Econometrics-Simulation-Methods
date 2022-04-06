

function [fullike]=typeprob0806(b,N,T,S,X,Like,nx);

%parse the b vector
binit=b(1:nx*(S-1),1);

num=ones(N,S);

s=2;

while s<S+1

    num(:,s)=exp(X*binit((s-2)*nx+1:(s-1)*nx));
    
    s=s+1;
    
end

prior=num./(sum(num,2)*ones(1,S));

per=exp(b(nx*(S-1)+1))./(1+exp(b(nx*(S-1)+1)));

trans=eye(S)*per+(1-eye(S))*((1-per)/(S-1));

s=1;
while s<S+1
    pfor(:,1,s)=prior(:,s).*Like(:,1,s);
    s=s+1;
end

t=2;
while t<T+1
    s2=1;
    while s2<S+1;
        temp=0;
        s1=1;
        while s1<S+1;
            temp=temp+pfor(:,t-1,s1)*trans(s1,s2);
            s1=s1+1;
        end
        pfor(:,t,s2)=temp.*Like(:,t,s2);
        s2=s2+1;
    end
    t=t+1;
end

fullike=-sum(log(sum(pfor(:,T,:),3)));
