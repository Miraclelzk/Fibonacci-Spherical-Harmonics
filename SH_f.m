%data  采样点数据
%theta phi采样点点位
%weight 权值
%输出的bw

function [coeffs]=SH_f(data, theta, phi, weight, bw)
% 每个点的球谐函数 
Y=SH(bw, theta, phi,'complex');
Y=Y.';

%计算正变换
coeffs=Y*(data.*weight);
end