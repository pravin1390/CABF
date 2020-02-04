function [sigma_r,sal] = sigmaMap_detail(f,sigma_min,sigma_max,smoothness)
%SIGMAMAP_DETAIL Summary of this function goes here
%   Detailed explanation goes here

if(~exist('sigma_min','var') || isempty(sigma_min))
    sigma_min = 1;
end
if(~exist('sigma_max','var') || isempty(sigma_max))
    sigma_max = 80;
end
if(~exist('smoothness','var') || isempty(smoothness))
    smoothness = 10;
end

sal = saliencyIG(f);

center = 0.5;
sigma_r = sigmoidMap(sal,center,smoothness,sigma_min,sigma_max);
sigma_r = round(sigma_r);

end

function [ y ] = sigmoidMap( x,x0,lambda,ymin,ymax )
% Transform input values by a sigmoid function of given parameters
% x = Arrayof input values
% x0 = Center of sigmoid
% lambda = Exponential parameter
% ymin = Infimum of range
% ymax = Supremum of range

if(~exist('x0','var')), x0 = 0.5; end
if(~exist('lambda','var')), lambda = 2; end
if(~exist('ymin','var')), ymin = 0; end
if(~exist('ymax','var')), ymax = 1; end

a = ymax - ymin;
y = a./(1 + exp(-lambda*(x-x0))) + ymin;

end

