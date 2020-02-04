clc;
clear;
close all;
%% input
filename = 'flower20.jpg';
I  =  imread(filename);
I = imresize(I,0.5);
I = double(I);

%% Construction of sigmamap and non uniform quantization of sigma values
sigmamap=sigmaMap_detail(I,5,80);
[m,n,d]=size(I);
L=4; %% Number of quantized levels 
[~,sigmacent] = kmeans(reshape(round(sigmamap),m*n,1),L);
sigs=5;
Iact=I./255;
fast_flag=1;

% 
%% direct bilateral filtering
tic,
Idirectbf=abf_color(I,sigs,sigmamap,'gaussian','zeros');    %exact bilateral
Tdirect=toc;
fprintf('Direct adaptive bilateral filter : \n');
fprintf('time for direct adaptive bilateral(ms)=%3.0f \n',Tdirect*1000);
%% Kmeans filtering
% Done in two steps : Clustering and Filtering
%% Filtering
Cluster=64;
%% ICIP algorithm with change in clustering, weights for approximation for range kernel and interpolation between clusters
tic,
%% Matlab clustering code     
[~,Centre] = rgb2ind(uint8(I(1:4:end,1:4:end,:)),Cluster,'nodither');

spatialtype='gaussian';     
convmethod='matlab'; % Change convmethod to 'O1' for O(1) convolutions
Ikmean=fastKmeansfiltapproxnystromsvd(Iact,sigs,sigmacent./255,sigmamap./255,Centre,spatialtype);      % bilateral kmeans
Ikmean=Ikmean.*255;
Ikmean(Ikmean<0)=0;
Ikmean(Ikmean>255)=255;
Tkmeans=toc;
fprintf('Fast adaptive bilateral filter by Kmeans complete with %d clusters \n',size(Centre,1));
fprintf('time for fast adaptive bilateral(ms)=%3.0f \n',Tkmeans*1000);

%% Detail enhancement
K = 3;      % Gain applied to detail layer (>1)
D_direct = I - Idirectbf;   % Detail layer computed using direct BF
D_kmean = I - Ikmean;       % Detail layer comppted using fast algorithm
E_direct = Idirectbf + K*D_direct;      % Enhanced using direct BF
E_kmean = Ikmean + K*D_kmean;
p_enhanced = psnr(E_kmean,E_direct,255)

%% Display enhanced images
figure; imshow(uint8(I)); pause(0.2); drawnow;
figure; imshow(uint8(E_direct));  pause(0.2); drawnow;
figure; imshow(uint8(E_kmean));  pause(0.2); drawnow;

