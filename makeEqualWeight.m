function [weight]=makeEqualWeight(n_)

    weight=ones(n_,1)*sqrt(2)*sqrt(2*pi)/n_*2*sqrt(pi);
end