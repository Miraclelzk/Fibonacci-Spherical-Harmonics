function [weight_f]=makeFibonacciWeight(bw_,n_)

    bw2_=2*bw_;
    % 每个点的球谐函数
    [theta, phi] = getFibonacci(n_);
    YF=SH(bw2_,theta,phi,'real');
    YF=conj(YF');

    %a是未知数的个数
    a=size(YF,1);
    AF=zeros(a,n_);
    for i=1:a 
        for j=1:n_    
            AF(i,j)=YF(i,j)*sqrt(2*pi);
        end 
    end

    %解权
    b=zeros(a,1);
    b(1)=1;
    xxF2 = lsqminnorm(AF,b,'warn');
    w1=xxF2*2*pi*sqrt(2);
    weight_f=real(w1);
end