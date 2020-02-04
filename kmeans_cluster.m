function [gIdx,c,dist,clust]=kmeans_cluster(X)
% function [gIdx,c,dist]=kmeans_cluster(X)
%   kmeans_modified    k-means clustring
%   [IDX, C, diatance] = kmeans_cluster(X, K) partititions the N x P data matrix X into K
%   clusters through a fully vectorized algorithm, where N is the number of
%   data points and P is the number of dimensions (variables). The
%   partition minimizes the sum of point-to-cluster-centroid Euclidean
%   distances of all clusters. 
%   IDX is the returned N x 1 vector contains the cluster indices of each point.
%   C is the K cluster centroid locations in the K x P matrix.
%   distance is maximum distance of data point in the cluster
%

[n,m]=size(X);
c=zeros(2,m);
Y=sum(X.^2,2);
[~,minindex]=min(Y);
[~,maxindex]=max(Y);
c(1,:)=X(minindex,:);
c(2,:)=X(maxindex,:);

% allocating variables
g0=ones(n,1);
gIdx=zeros(n,1);
D=zeros(n,2);

% Main loop converge if previous partition is the same as current
while any(g0~=gIdx)
    g0=gIdx;
    % Loop for each centroid
    for t=1:2
          D(:,t)=sum((X-c(t,:)).^2,2);
    end
    % Partition data to closest centroids
    [z,gIdx]=min(D,[],2);
    % Update centroids using means of partitions
    for t=1:2
        c(t,:)=mean(X(gIdx==t,:)); 
        clust(t)=size(X(gIdx==t,:),1); 
    end
end
    for t=1:2
          dist(t)=sum(z(gIdx==t));
    end
end
