clear
%系数的阶数
bw_f=32;
nf=4*(bw_f+2)*(bw_f+2);
disp(nf)
w_f=makeFibonacciWeight(bw_f,nf);

[theta, phi] = getFibonacci(nf);


%%
%循环

%循环次数
num=40;

%MeanSE_e=zeros(40,1);
%MaxSE_e=zeros(40,1);
RMSE_e=zeros(40,1);
MAE_e=zeros(40,1);
for i=1:num
    %生成随机系数
    coeff =zeros(bw_f*bw_f,1);
    for l=1:bw_f
        coeff(l*l-l+1)=complex(randn(),0);
        for m=1:l-1
            t1=randn();
            t2=randn();
            coeff(l*l-l+1+m)=complex(t1,t2);
            coeff(l*l-l+1-m)=complex((-1)^m*t1,(-1)^(m+1)*t2);
        end
    end

    %逆变换

    data=SH_inv(coeff,bw_f,theta,phi);
    new_coeff=SH_f(data,theta, phi,w_f,bw_f);

    error=new_coeff-coeff;
    %均方误差
    %MeanSE_e(i)=sum(abs(error).^2)/num;
    %MaxSE_e(i)=max(abs(error)).^2;
      
    RMSE_e(i)=sqrt(sum(abs(error).^2)/num);
    MAE_e(i)=sum(abs(error))/num;
end

%% draw

x=(1:40);

figure
plot(x,RMSE_e,'-o');
hold on;
plot(x,MAE_e,'-*');


%legend('经纬网均方误差','经纬网最大方误差','斐波那契均方误差','斐波那契最大方误差');
legend('RMSE','MAE');
axis([-inf inf 0  3e-13]);

ylabel('Deviation');
xlabel('Number of Iterations');
set(gca, 'FontSize', 18)
set(gcf, 'Position', [10 10 500 400]);