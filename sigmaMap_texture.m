function [sigma_r] = sigmaMap_texture(f,sigma_min,sigma_max,smoothness,winsize)
% Compute sigma_r map for texture filtering
% f = Input image (grayscale/RGB)
% sigma_min = Min. value of sigma_r (default = 25)
% sigma_max = Max. value of sigma_r (default = 80)
% smoothness = Parameter (>=0) to control smoothness of sigma values. If
%   smoothness=0, the map is constant. If smoothness = inf, the map is
%   binary (either min or max value). Default = 8.
% winsize = Window size to compute ralative total variation (RTV) (optional, default = 7)

if(~exist('sigma_min','var') || isempty(sigma_min))
    sigma_min = 25;
end
if(~exist('sigma_max','var') || isempty(sigma_max))
    sigma_max = 80;
end
if(~exist('smoothness','var') || isempty(smoothness))
    smoothness = 8;
end
if(~exist('winsize','var') || isempty(winsize))
    winsize = 7;
end
if(size(f,3)==3)
    f = double(rgb2gray(uint8(f)));
end

[~,E] = fastSRTV(f,winsize);
E = E/max(E(:));

center = 0.7;
sigma_r = sigmoidMap(1-E,center,smoothness,sigma_min,sigma_max);
sigma_r = imdilate(sigma_r,strel('disk',2,4));  % Clean up the fine noise

end

function [ srtv,mrtv ] = fastSRTV( I,k )
% Smoothed modified relative total variation
% I = Input image (double)
% k = Window length
% srtv = Smoothed MRTV
% mrtv = MRTV

if(size(I,3)==3)
    I = double(rgb2gray(uint8(I)));
end

I = imfilter(I,ones(5)/25,'symmetric');

use_delta = true;   % Set to false if tonal range should not be used in MRTV
if(use_delta)
    Imin = minFilter(double(I),k);
    Imax = maxFilter(double(I),k);
    Delta = Imax - Imin;
else
    Delta = 1;
end

% Compute MRTV
[gradI,~] = imgradient(I,'centralDifference');
num = maxFilter(gradI,k);
den = imfilter(gradI,ones(k),'symmetric');
mrtv = Delta.*num./(den + 1e-9);
mrtv(den==0) = 0;

% Apply smoothing to MRTV
srtv = minFilter(mrtv,k);

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

function [ Rmin ] = minFilter( fin,w )
% Filter to find local minimum in every pixel neighborhood

sym    = (w - 1)/2;
[minput, ninput] = size(fin);
rowpad=(ceil(minput/w)*w)-minput;
columnpad=(ceil(ninput/w)*w)-ninput;
template=padarray(fin,[rowpad,columnpad],Inf,'post');
[m,n]=size(template);
Rmin = nan(minput,ninput);

% scan along rows
for ii = 1 : minput
    L     = zeros(n, 1);
    R     = zeros(n, 1);
    L(1)  = template(ii, 1);
    R(n)  = template(ii, n);
    for k = 2 : n
        if  mod(k - 1, w) == 0
            L(k)          = template(ii ,  k);
            R(n - k + 1)  = template(ii ,  n - k + 1);
        else
            L(k)          = min( L(k-1) , template(ii, k) );
            R(n - k + 1)  = min( R(n - k + 2), template(ii, n - k + 1) );
        end
    end
    for k = 1 : n
        p = k - sym;
        q = k + sym;
        if p < 1
            r = Inf;
        else
            r = R(p);
        end
        if q > n
            l = Inf;
        else
            l = L(q);
        end
        template(ii, k) = min(r,l);
    end
end
% scan along columns
for jj = 1 : ninput
    L    = zeros(m, 1);
    R    = zeros(m, 1);
    L(1) = template(1, jj);
    R(m) = template(m, jj);
    for k = 2 : m
        if  mod(k - 1, w) == 0
            L(k)          = template(k, jj);
            R(m - k + 1)  = template(m - k + 1, jj);
        else
            L(k)          = min( L(k - 1), template(k, jj) );
            R(m - k + 1)  = min( R(m - k + 2), template(m - k + 1, jj));
        end
    end
    for k = 1 : m
        p = k - sym;
        q = k + sym;
        if p < 1
            r = Inf;
        else
            r = R(p);
        end
        if q > m
            l = Inf;
        else
            l = L(q);
        end
        if (k<=minput)
            Rmin(k,jj) = min(r,l);
        end
    end
end

end

function [ Rmax ] = maxFilter( fin,w )
% Filter to find local maximum in every pixel neighborhood

sym    = (w - 1)/2;
[minput, ninput] = size(fin);
rowpad=(ceil(minput/w)*w)-minput;
columnpad=(ceil(ninput/w)*w)-ninput;
template=padarray(fin,[rowpad,columnpad],'post','replicate');
[m,n]=size(template);
Rmax = nan(minput,ninput);

% scan along rows
for ii = 1 : minput
    L     = zeros(n, 1);
    R     = zeros(n, 1);
    L(1)  = template(ii, 1);
    R(n)  = template(ii, n);
    for k = 2 : n
        if  mod(k - 1, w) == 0
            L(k)          = template(ii ,  k);
            R(n - k + 1)  = template(ii ,  n - k + 1);
        else
            L(k)          = max( L(k-1) , template(ii, k) );
            R(n - k + 1)  = max( R(n - k + 2), template(ii, n - k + 1) );
        end
    end
    for k = 1 : n
        p = k - sym;
        q = k + sym;
        if p < 1
            r = -1;
        else
            r = R(p);
        end
        if q > n
            l = -1;
        else
            l = L(q);
        end
        template(ii, k) = max(r,l);
    end
end

% scan along columns
for jj = 1 : ninput
    L    = zeros(m, 1);
    R    = zeros(m, 1);
    L(1) = template(1, jj);
    R(m) = template(m, jj);
    for k = 2 : m
        if  mod(k - 1, w) == 0
            L(k)          = template(k, jj);
            R(m - k + 1)  = template(m - k + 1, jj);
        else
            L(k)          = max( L(k - 1), template(k, jj) );
            R(m - k + 1)  = max( R(m - k + 2), template(m - k + 1, jj));
        end
    end
    for k = 1 : m
        p = k - sym;
        q = k + sym;
        if p < 1
            r = -1;
        else
            r = R(p);
        end
        if q > m
            l = -1;
        else
            l = L(q);
        end
        if (k<=minput)
            Rmax(k,jj) = max(r,l);
        end
    end
end

end



