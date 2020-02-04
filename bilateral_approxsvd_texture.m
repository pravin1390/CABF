% clc;
clear;
close all;
%% input
filename = 'input.jpg';
I  =  double(imread(filename));

%% Construction of sigmamap and non uniform quantization of sigma values
sigmamap=sigmaMap_texture(I);
[m,n,d]=size(I);
L=4; %% Number of quantized levels 
[~,sigmacent] = kmeans(reshape(round(sigmamap),m*n,1),L);
sigs=3;
Iact=I./255;
fast_flag=1;

% Three iterations of filtering are done here

%% direct bilateral filtering
tic,
Idirectbf=abf_color(I,sigs,sigmamap,'gaussian','zeros');    %exact bilateral
Idirectbf=abf_color(Idirectbf,sigs,sigmamap,'gaussian','zeros');    %exact bilateral
Idirectbf=abf_color(Idirectbf,sigs,sigmamap,'gaussian','zeros');    %exact bilateral
Tdirect=toc;
fprintf('Direct adaptive bilateral filter : \n');
fprintf('time for direct adaptive bilateral(ms)=%3.0f \n',Tdirect*1000);

%% Kmeans filtering 
% Done in two steps : Clustering and Filtering
%% Filtering
Cluster=16;
%% ICIP algorithm with change in clustering, weights for approximation for range kernel and interpolation between clusters
tic,
%% Bisecting K-means clustering
for i=1:d
    Iact2(:,:,i)=imresize(Iact(:,:,i),[256 256],'nearest');
end
Ares=reshape(Iact2,size(Iact2,1)*size(Iact2,2),d);
Centre=kmeans_recursive(Ares,Cluster);

spatialtype='gaussian';     
convmethod='matlab'; % Change convmethod to 'O1' for O(1) convolutions
Ikmean=fastKmeansfiltapproxnystromsvd(Iact,sigs,sigmacent./255,sigmamap./255,Centre,spatialtype);      % bilateral kmeans
for i=1:d
    Iact2(:,:,i)=imresize(Ikmean(:,:,i),[256 256],'nearest');
end
Ares=reshape(Iact2,size(Iact2,1)*size(Iact2,2),d);
Centre=kmeans_recursive(Ares,Cluster);
Ikmean=fastKmeansfiltapproxnystromsvd(Ikmean,sigs,sigmacent./255,sigmamap./255,Centre,spatialtype);      % bilateral kmeans
for i=1:d
    Iact2(:,:,i)=imresize(Ikmean(:,:,i),[256 256],'nearest');
end
Ares=reshape(Iact2,size(Iact2,1)*size(Iact2,2),d);
Centre=kmeans_recursive(Ares,Cluster);
Ikmean=fastKmeansfiltapproxnystromsvd(Ikmean,sigs,sigmacent./255,sigmamap./255,Centre,spatialtype);      % bilateral kmeans
Ikmean=Ikmean.*255;
Ikmean(Ikmean<0)=0;
Ikmean(Ikmean>255)=255;
Tkmeans=toc;
fprintf('Fast adaptive bilateral filter by Kmeans complete with %d clusters \n',size(Centre,1));
fprintf('time for fast adaptive bilateral(ms)=%3.0f \n',Tkmeans*1000);
error2 = reshape(Idirectbf-Ikmean, [d*m*n,1]);
MSE_mcbf2 = sqrt(sum(error2.^2)/(d*m*n));
PSNR2=20*log10(255/(MSE_mcbf2));
fprintf('mean sq error=%f, PSNR = %f db  \n',MSE_mcbf2,PSNR2);


% % %% output
figure;
imshow(uint8(I));title('Original image');
figure;
imshow(uint8(Idirectbf));title('direct bilateral filter');
figure;
imshow(uint8(Ikmean));title('fast bilateral filter');
figure; imagesc(sigmamap); axis image; axis off; colorbar;
