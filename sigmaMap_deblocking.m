function [ sigma_r ] = sigmaMap_deblocking( f,sigma_min,k0 )
% Compute sigma_r map for JPEG deblocking
% f = Input image (grayscale/RGB)
% sigma_min = Min. value of sigma (optional, default = 25)
% k0 = Gain (optional, default = 1)

if(~exist('sigma_min','var') || isempty(sigma_min))
    sigma_min = 25;
end
if(~exist('k0','var') || isempty(k0))
    k0 = 1;
end

if(size(f,3)==3)
    f = double(rgb2gray(uint8(f)));
end

[M,N] = size(f);

[bmask,vmask,hmask,cornermask,centermask] = edgeMask(M,N);
[R,C] = find(bmask);

% Find discontinuities at block boundaries
bV = imfilter(f,[-1,0,1],'symmetric');
bH = imfilter(f,[-1,0,1]','symmetric');

V = nan(M,N);
V(cornermask) = max(abs(bV(cornermask)),abs(bH(cornermask)));
V(vmask) = abs(bV(vmask));
V(hmask) = abs(bH(hmask));
V(centermask) = 0;

Vdata = V(bmask);
F = scatteredInterpolant(R,C,Vdata);
[RR,CC] = find(true(size(f)));
BDmap = F(RR,CC);   % Block discontinuity map
BDmap = reshape(BDmap,size(f));
sigma_r = max(sigma_min,k0*BDmap);

end

function [ bmask,vmask,hmask,cornermask,centermask ] = edgeMask(M,N)
Br = floor(M/8);
Bc = floor(N/8);

vtile = false(8);
vtile(2:7,[1,8]) = true;

htile = false(8);
htile([1,8],2:7) = true;

crtile = false(8);
crtile(1,1) = true;
crtile(1,8) = true;
crtile(8,1) = true;
crtile(8,8) = true;

cntile = false(8);
cntile(4:5,4:5) = true;

btile = vtile | htile | crtile | cntile;

vmask = repmat(vtile,Br+1,Bc+1);
vmask = vmask(1:M,1:N);
hmask = repmat(htile,Br+1,Bc+1);
hmask = hmask(1:M,1:N);
cornermask = repmat(crtile,Br+1,Bc+1);
cornermask = cornermask(1:M,1:N);
centermask = repmat(cntile,Br+1,Bc+1);
centermask = centermask(1:M,1:N);
bmask = repmat(btile,Br+1,Bc+1);
bmask = bmask(1:M,1:N);

end



