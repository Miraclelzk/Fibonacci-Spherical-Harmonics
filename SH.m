function [basis] = SH(degree, theta, phi, real_or_complex)
    %degree 阶数
    %theta phi采样点点位
    %返回 球谐函数*采样位置  的矩阵

    %theta是个数组  phi也是

    k =degree*degree;  %% k is the # of coefficients
    n = size(theta,1);
    Y = zeros(n,k);

    for l = 0:degree-1
        % calculate the spherical harmonics
        Pm = legendre(l,cos(theta')); % legendre part
        Pm = Pm';
        lconstant = sqrt((2*l + 1)/(4*pi));

        center = (l+1)*(l+1) - l;

        Y(:,center) = lconstant*Pm(:,1);
        for m=1:l
            precoeff = lconstant * sqrt(factorial(l - m)/factorial(l + m));

            switch lower(real_or_complex)
                case 'real'
                     if mod(m,2) == 1
                         precoeff = -precoeff;
                     end
                    Y(:,center + m) = sqrt(2)*precoeff*Pm(:,m+1).*cos(m*phi);
                    Y(:,center - m) = sqrt(2)*precoeff*Pm(:,m+1).*sin(m*phi);
                case 'complex'
                    if mod(m,2) == 1
                        Y(:,center + m) = precoeff*Pm(:,m+1).*exp(1i*m*phi);
                        Y(:,center - m) = -precoeff*Pm(:,m+1).*exp(-1i*m*phi);
                    else
                        Y(:,center + m) = precoeff*Pm(:,m+1).*exp(1i*m*phi);
                        Y(:,center - m) = precoeff*Pm(:,m+1).*exp(-1i*m*phi);
                    end
                otherwise
                    error('The last argument must be either \"real\" (default) or \"complex\".');
            end
        end

    end

    basis = Y;
end