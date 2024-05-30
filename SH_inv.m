function [data]=SH_inv(coeffs,bw,theta,phi)

Y=SH(bw,theta,phi,'complex');
Y=conj(Y);

%计算正变换
data=Y*coeffs;
end