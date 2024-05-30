addpath 'code'
clear

%% 参数
bw_f=32;
bw_e=32;
%% 生成随机系数
coeff =zeros(bw_e*bw_e,1);
for l=1:bw_e
    coeff(l*l-l+1)=complex(randn(),0);
    for m=1:l-1
        t1=randn();
        t2=randn();
        coeff(l*l-l+1+m)=complex(t1,t2);
        coeff(l*l-l+1-m)=complex((-1)^m*t1,(-1)^(m+1)*t2);
    end
end


%% 斐波那契采样
nf=4*(bw_f+2)*(bw_f+2);
disp(nf)
w_f=makeFibonacciWeight(bw_f,nf);

[theta1, phi1] = getFibonacci(nf);
[t(:,1),t(:,2),t(:,3)]= sph2cart(phi1,pi/2-theta1,1.0);
% plot3(t(:,1),t(:,2),t(:,3),'k.');
data_f=SH_inv(coeff,bw_f,theta1,phi1);
new_f_coeff =SH_f(data_f,theta1, phi1,w_f,bw_f);


error_f=new_f_coeff-coeff;
tmpmag_f=abs(error_f);
curmax_f=max(tmpmag_f);
cursum_f=sum(tmpmag_f);
curmean_f=mean(tmpmag_f);
disp(curmax_f);
disp(cursum_f);
disp(curmean_f);

%% 斐波那契采样  
nf=4*(bw_f+2)*(bw_f+2);
disp(nf)
w_f=makeEqualWeight(nf);

[theta1, phi1] = getFibonacci(nf);
[t(:,1),t(:,2),t(:,3)]= sph2cart(phi1,pi/2-theta1,1.0);
% plot3(t(:,1),t(:,2),t(:,3),'k.');
data_f=SH_inv(coeff,bw_f,theta1,phi1);
new_f_coeff =SH_f(data_f,theta1, phi1,w_f,bw_f);


error_f=new_f_coeff-coeff;
tmpmag_f=abs(error_f);
curmax_f=max(tmpmag_f);
cursum_f=sum(tmpmag_f);
curmean_f=mean(tmpmag_f);
disp(curmax_f);
disp(cursum_f);
disp(curmean_f);



%% Icosahedral  二十面体
k = 4;
ico_N = 10*4^k+2;
disp(ico_N)
[xyz, ~] = getIcosNodes(k,0);
% plot3(xyz(:,1),xyz(:,2),xyz(:,3),'k.');
[phi2,theta2,~] = cart2sph(xyz(:,1),xyz(:,2),xyz(:,3));
theta2 = pi/2-theta2;
w_ico=makeEqualWeight(ico_N);

data_ico=SH_inv(coeff,bw_f,theta2,phi2);
new_ico_coeff =SH_f(data_ico,theta2, phi2,w_ico,bw_f);


error_ico=new_ico_coeff-coeff;
tmpmag_f=abs(error_ico);
curmax_f=max(tmpmag_f);
cursum_f=sum(tmpmag_f);
curmean_f=mean(tmpmag_f);
disp(curmax_f);
disp(cursum_f);
disp(curmean_f);

%% 经纬度采样

lat = bw_e*2;
lon = bw_e*2;
n=lon*lat;
disp(n)
thetastep = pi/lat/2;
phistep=2*pi/lon/2;

theta3=zeros(n,1);
phi3=zeros(n,1);
for i=0:lat-1  
    for j=0:lon-1  
        theta3(i*lon+j+1)=(2*i+1)*thetastep;
    	phi3(i*lon+j+1)=(2*j+1)*phistep;   
    end   
end

w_e=makeEquiangularWeight(bw_e,lon,lat);
data_e=SH_inv(coeff,bw_e,theta3,phi3);
new_e_coeff=SH_f(data_e,theta3, phi3,w_e,bw_e);

error_e=new_e_coeff-coeff;
tmpmag_e=abs(error_e);
curmax_e=max(tmpmag_e);
cursum_e=sum(tmpmag_e);
curmean_e=mean(tmpmag_e);
disp(curmax_e);
disp(cursum_e);
disp(curmean_e);


%% HEALPix

S = 20;
HEALPix_N = 12*S^2;
disp(HEALPix_N)
HEALPix = getHEALPixNodes(S);
% plot3(HEALPix(:,1),HEALPix(:,2),HEALPix(:,3),'k.');
[phi4,theta4,~] = cart2sph(HEALPix(:,1),HEALPix(:,2),HEALPix(:,3));
theta4 = pi/2-theta4;
w_HEALPix=makeEqualWeight(HEALPix_N);

data_HEALPix=SH_inv(coeff,bw_f,theta4,phi4);
new_HEALPix_coeff =SH_f(data_HEALPix,theta4, phi4,w_HEALPix,bw_f);


error_HEALPix=new_HEALPix_coeff-coeff;
tmpmag_f=abs(error_HEALPix);
curmax_f=max(tmpmag_f);
cursum_f=sum(tmpmag_f);
curmean_f=mean(tmpmag_f);
disp(curmax_f);
disp(cursum_f);
disp(curmean_f);
