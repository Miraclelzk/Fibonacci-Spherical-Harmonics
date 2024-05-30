% 计算等角网格的权值
function [weight_e]=makeEquiangularWeight(bw_,lon_,lat_)
    %解析权
    fudge = pi/(2*lat_);
    weight=ones(lat_,1);
    for i=0:lat_-1
       tmpsum=0;
       for j=0:bw_-1
            tmpsum=tmpsum+1/(2*j+1)*sin((2*j+1)*(2*i+1)*fudge);
       end

       tmpsum=tmpsum*sin((2*i+1)*fudge);
       tmpsum=tmpsum*2*pi*4/lat_/lon_;
       weight(i+1)=tmpsum; 
    end
    weight_e=weight(:,ones(lon_,1))';
    weight_e=reshape(weight_e,lon_*lat_,1);    
end