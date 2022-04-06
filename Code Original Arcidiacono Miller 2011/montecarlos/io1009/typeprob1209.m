
function [PType,trans]=typeprob1209(b,X,nx,N,T,S,Like);

%the program takes as inputs:
%1) the coefficients
%2) the variables that affect the initial conditions, X
%3) the number of columns of X, the number of observations in the data set (N), the number of time
%periods (T), and the number of unobserved states (S)
%4) the individual components of the likelihood (Like)


%
%the program returns:
%1)the conditional probability of each observation being in one of the
%observed states, PType (N x T x S) 
%2) the updated initial conditions and the updated unobservable transition
%matrix



%In order to update the type probabilities at time t, need to be able to
%incorporate information from periods after t and from periods before t
%
%pback gives the contribution from periods after t

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


pback=ones(N,T,S);
pback(:,T,:)=Like(:,T,:);


t=T-1;
while t>0
    s1=1;
    while s1<S+1
        temp=0;
        s2=1;
        while s2<S+1
            temp=temp+pback(:,t+1,s2)*trans(s1,s2);
            s2=s2+1;
        end
        pback(:,t,s1)=temp.*Like(:,t,s1);
        s1=s1+1;
    end
    t=t-1;
end

%pfor gives the contribution from periods before t

pfor=ones(size(Like));

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

%Using pfor and pack we can now calculate the conditional probability of
%being in each unobserved state for each observation

PType=zeros(size(Like));

dem=sum(pfor.*pback./Like,3);

s=1;
while s<S+1
    PType(:,:,s)=(pfor(:,:,s).*pback(:,:,s)./Like(:,:,s))./dem;
    s=s+1;
end

PType=reshape(PType,(N*T*S),1);
