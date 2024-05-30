function [theta, phi] = getFibonacci(n)

%斐波那契采样采样的每个点的位置
golden=(sqrt(5)-1)/2;
theta=zeros(n,1);
phi=zeros(n,1);

for i=0:n-1
    d1=(1-n)/2+i;
    theta(i+1)= asin(2.*d1./n) + pi/2.0 ;  
    phi(i+1) = 2.0 * pi*(d1*golden - floor(d1*golden));
end  

end