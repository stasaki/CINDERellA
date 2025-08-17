function [nlogL] = my_normlike(res,mu,sigma)

zscore = (res - mu) ./ sigma;
L = -.5.*zscore.*zscore - log(sqrt(2.*pi).*sigma);
nlogL = -sum(L);
