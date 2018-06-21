function [ strong_users, users_noma, users_oma ] = pair_kmeans_A()

global P;
global Users;
global Pairs;
global Cluster;

% 随机选择 1 个用户作为初始聚类的质心
seed = randperm(P.nums, 1);
Cluster(1).centroid = Users(seed).h.';

% 计算每个用户到所有质心的相似度
for u = 1:P.nums
    hu = Users(u).h.';
    cos_vec = zeros(1, length(Cluster));
    selected_rate = zeros(1, P.nums);
    for c = 1:1:length(Cluster)
        hc = Cluster(c).centroid;
        cos_vec(c) = norm(hu*hc') / norm(hu) / norm(hc);
    end
    [~, index] = max(cos_vec(u,:));
    selected_rate(u) = cos_vec(index)^2;
    
    % 添加用户到相似度最高的质心
    Cluster(index).fans = [Cluster(index).fans, u];
end

[~, next] = max(selected_rate);



% 根据 Cluster 进行 pair


end

