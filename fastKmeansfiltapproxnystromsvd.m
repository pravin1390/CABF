function outimg=fastKmeansfiltapproxnystromsvd(img,S,sigmacent,sigmamap,Centre,spatialkernel,imgguide)
if ~exist('imgguide','var')
     % guideimage is same as inputimage
     imgguide=img;
end
[m,n,~]=size(img);
guided=size(imgguide,3);
Cluster=size(Centre,1);
L=length(sigmacent);
num=zeros(size(img));
den=zeros(m,n);
imgguide=reshape(imgguide,m*n,guided);
sigmamap=reshape(sigmamap,m*n,1);
    
if strcmp(spatialkernel,'box')
    filt     = ones(2*S+1,2*S+1);       
elseif strcmp(spatialkernel,'gaussian')       
    w  = round(6*S); if (mod(w,2) == 0); w  = w+1; end
    filt     = fspecial('gaussian', [w w], S);
else
end  
B=zeros(Cluster,Cluster);
for i=1:Cluster
    B(i,i)=0;
    for k=i+1:Cluster
        B(i,k)=sum((Centre(i,:)-Centre(k,:)).^2,2);
        B(k,i)=B(i,k);
    end
end
A=[];
for j=1:L
    A=[A;exp(-B/(2*sigmacent(j)^2))];
end
[U,D,V]=svd(A,'econ');
Dhat=diag(D);    
    Eigmaty=zeros(m,n,Cluster);
    Wx=zeros(m*n,Cluster);
for i=1:Cluster    
    Wx(:,i)=sum((imgguide-Centre(i,:)).^2,2);
end
for i=1:L
    Eigmaty=Eigmaty+reshape(exp(-Wx/(2*(sigmacent(i)^2)))*U(Cluster*(i-1)+1:Cluster*(i-1)+Cluster,:),m,n,Cluster);
end
Wx=exp(-bsxfun(@rdivide,Wx,(2*(sigmamap.^2))));

%% Calculating Bilateral image for each cluster centre as index pixel
    for i=1:Cluster
        Eigmatx=(1/Dhat(i,1))*reshape(Wx*V(:,i),m,n);
        den=den+bsxfun(@times,Eigmatx,imfilter(Eigmaty(:,:,i),filt));
        num=num+bsxfun(@times,Eigmatx,imfilter(bsxfun(@times,Eigmaty(:,:,i),img),filt));  
    end
    outimg=bsxfun(@rdivide,num,den);
end


