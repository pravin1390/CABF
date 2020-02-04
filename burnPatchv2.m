function [ hfig ] = burnPatchv2( f,cen,width,height,zoomflag,scale,placement,color )
%BURNPATCHV2
% Version 2
% f = Input image
% cen = N-by-2 array containing patch centers like [r1,c1;r2,c2;r3,c3]
% width = Scalar OR 1-by-N array containing patch widths
% height = Scalar OR 1-by-N array containing patch heights
% zoomflag = FALSE (default) if zoomed version of patches not required,
%            TRUE otherwise
% scale = Scalar OR 1-by-N array containing zoom factor of patches
% placement = 1-by-N cell array containing entries 'bl','br','tl','tr',
%             e.g. {'bl','tr',tl'}
% color = 1-by-N cell array containing color names of highlighting
%         rectangles, e.g. {'r','g','y'}

N = size(cen,1);
if(~exist('zoomflag','var'))
    zoomflag = false;
end
if(isscalar(width))
    width = width*ones(1,N);
end
if(isscalar(height))
    height = height*ones(1,N);
end
if(isscalar(scale))
    scale = scale*ones(1,N);
end
if(length(width)~=N || length(height)~=N)
    error('Specify %d widths and heights\n',N);
end
if(zoomflag)
    if(N==1 && ~iscell(placement))
        placement = {placement};
    end
    if(N>1 && ~iscell(placement))
        error('Specify %d placement locations in a cell array\n',N);
    end
end
if(~iscell(color))
    color = repmat({color},1,N);
end


if(zoomflag)    % Zoom and highlight patches by rectangles
    [RR,CC,~] = size(f);
    offset = 2;
    box_patch = cell(1,N);
    box_zoom = cell(1,N);
    for k = 1:N
        [P,left,top] = getPatch(f,cen(k,:),width(k),height(k));
        box_patch{k} = [left-1,top-1,width(k),height(k)];
        Pbig = imresize(P,scale(k));
        [rr,cc,~] = size(Pbig);
        if(strcmp(placement{k},'bl'))
            r2 = RR-offset; r1 = r2-rr+1;
            c1 = 1+offset; c2 = c1+cc-1;
            f(r1:r2,c1:c2,:) = Pbig;
        elseif(strcmp(placement{k},'br'))
            r2 = RR-offset; r1 = r2-rr+1;
            c2 = CC-offset; c1 = c2-cc+1;
            f(r1:r2,c1:c2,:) = Pbig;
        elseif(strcmp(placement{k},'tl'))
            r1 = 1+offset; r2 = r1+rr-1;
            c1 = 1+offset; c2 = c1+cc-1;
            f(r1:r2,c1:c2,:) = Pbig;
        elseif(strcmp(placement{k},'tr'))
            r1 = 1+offset; r2 = r1+rr-1;
            c2 = CC-offset; c1 = c2-cc+1;
            f(r1:r2,c1:c2,:) = Pbig;
        end
        box_zoom{k} = [c1-1,r1-1,cc,rr];
    end
    hfig = figure;
    imshow(f);
    for k = 1:N
        rectangle('Position',box_patch{k},'EdgeColor',color{k},'LineWidth',1.2);
        rectangle('Position',box_zoom{k},'EdgeColor',color{k},'LineWidth',1.2);
    end
else    % Only highlight patches by rectangles
    box_patch = cell(1,N);
    for k = 1:N
        [~,left,top] = getPatch(f,cen(k,:),width(k),height(k));
        box_patch{k} = [left-1,top-1,width(k),height(k)];
    end
    hfig = figure;
    imshow(f);
    for k = 1:N
        rectangle('Position',box_patch{k},'EdgeColor',color{k},'LineWidth',1.2);
    end
end

end

function [ P,left,top ] = getPatch(f,cen,w,h)

hw = (w-1)/2;
hh = (h-1)/2;
cen_r = cen(1);
cen_c = cen(2);
left = cen_c - hw;
top = cen_r - hh;
P = f(cen_r-hh:cen_r+hh,cen_c-hw:cen_c+hw,:);

end







