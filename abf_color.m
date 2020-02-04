function [ Bf ] = abf_color( f,sigma_s,sigma_r,filttype,padtype )
%ABF_COLOR Adaptve bilateral filter for color images

f = double(f);
[rr, cc, ~] = size(f);
if(isscalar(sigma_r))
    sigma_r = sigma_r * ones(rr,cc);
end
if(~exist('filttype','var'))
    filttype = 'gaussian';
end
if(~exist('padtype','var'))
    padtype = 'symmetric';
end
if(strcmp(filttype,'gaussian'))
    kerrad = 3*sigma_s;
elseif(strcmp(filttype,'box'))
    kerrad = sigma_s;
end

% Initialize
W = zeros([rr,cc,3]);
Z = zeros([rr,cc]);
sigma_r2 = sigma_r.*sigma_r;

if(strcmp(padtype,'zeros'))
    f = padarray(f,[kerrad, kerrad, 0]);
elseif(strcmp(padtype,'symmetric'))
    f = padarray(f,[kerrad, kerrad, 0],'symmetric');  % Pad image
end

% Gaussian spatial kernel
if(strcmp(filttype,'gaussian'))
    spker = fspecial('gaussian',2*kerrad+1,sigma_s);
elseif(strcmp(filttype,'box'))
    spker = ones(2*kerrad+1);
end

% Start implementation
for j1 = kerrad+1:kerrad+rr
    for j2 = kerrad+1:kerrad+cc
        nb = f(j1-kerrad:j1+kerrad,j2-kerrad:j2+kerrad,:);
        r_arg = sum((nb-f(j1,j2,:)).^2,3);
        rker = exp(-0.5*r_arg/(sigma_r2(j1-kerrad,j2-kerrad)));
        W(j1-kerrad,j2-kerrad,:) = sum(sum(repmat(spker.*rker,1,1,3) .* nb,1),2);
        Z(j1-kerrad,j2-kerrad) = sum(sum(spker .* rker));
    end
end

Bf = nan(size(W));
for k = 1:3
    Bfk = nan(size(Z));
    Wk = W(:,:,k);
    mask = (Z==0);
    Bfk(mask) = Wk(mask);
    Bfk(~mask) = Wk(~mask)./Z(~mask);
    Bf(:,:,k) = Bfk;
end

end


