% clc;
clear;
close all;
%% input
filename = 'pirate_q10.jpg';
I  =  double(imread(filename));

%% Construction of sigmamap and non uniform quantization of sigma values
sigmamap=sigmaMap_deblocking(I,15,1);
[m,n,d]=size(I);
L=4; %% Number of quantized levels 
[~,sigmacent] = kmeans(reshape(round(sigmamap),m*n,1),L);
sigs=3;
Iact=I./255;
fast_flag=1;


%% direct bilateral filtering
tic,
Idirectbf=abf_bruteforce(I,sigs,sigmamap,I,'gaussian');    %exact bilateral
Tdirect=toc;
fprintf('Direct adaptive bilateral filter : \n');
fprintf('time for direct adaptive bilateral(ms)=%3.0f \n',Tdirect*1000);
%% Kmeans filtering
% Done in two steps : Clustering and Filtering
%% Filtering
Cluster=8;
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
Ikmean=Ikmean.*255;
Ikmean(Ikmean<0)=0;
Ikmean(Ikmean>255)=255;
Tkmeans=toc;
fprintf('Fast adaptive bilateral filter by Kmeans complete with %d clusters \n',size(Centre,1));
fprintf('time for fast adaptive  bilateral(ms)=%3.0f \n',Tkmeans*1000);

% % %% output
figure;
imshow(uint8(I));%title('Original image');
figure;
imshow(uint8(Idirectbf));%title('direct bilateral filter');
figure;
imshow(uint8(Ikmean));%title('fast bilateral filter');
