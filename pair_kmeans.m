function [ strong_users, users_noma, users_oma ] = pair_kmeans( K )

global P;
global Users;
global Pairs;
global Cluster;

% 随机选择 K 个用户作为 K 个簇的质心
seed = randperm(P.nums, K);

for i = 1:1:length(seed)
    Cluster(i).centroid = Users(seed(i)).h.';
end

last_centroid = [Cluster(1:length(Cluster)).centroid];

while 1
    % 清空结果，方便重新计算
    for i = 1:1:length(Cluster)
       Cluster(i).fans = []; 
    end
    
    % 计算每个用户到 K 个质心的相似度
    for u = 1:P.nums
        hu = Users(u).h.';
        cos_vec = zeros(1, length(Cluster));
        for c = 1:1:length(Cluster)
            hc = Cluster(c).centroid;
            cos_vec(c) = norm(hu*hc') / norm(hu) / norm(hc);
        end
        [~, index] = max(cos_vec);
        
        % 添加用户到相似度最高的质心
        Cluster(index).fans = [Cluster(index).fans, u];
    end
    
    % 重新计算每个簇的质心
    for c = 1:length(Cluster)
        fans = Cluster(c).fans;
        mu = zeros(1, P.Nt);
        for f = fans
            mu = mu + Users(f).h.';
        end
        mu = mu / length(fans);
        Cluster(c).centroid = mu;
    end
    
    % 如果分簇结果不再发生变化，就提前退出
    curr_centroid = [Cluster(1:length(Cluster)).centroid];
    if all(last_centroid == curr_centroid)
       break; 
    end
    last_centroid = curr_centroid;
end

% 根据 Cluster 进行 pair
p = 1;
for c = 1:1:length(Cluster)
    fans = Cluster(c).fans;
    pls = [];
    for u = fans
        pls = [pls, Users(u).pathloss];
    end
    [seed, index] = sort(pls);
    span = length(seed)/2;
    for i = 1:1:span
        strong = fans(index(length(seed)-i+1));
        weak = fans(index(i));
        Pairs(p).pair = [strong, weak];
        Pairs(p).strong_user = strong;
        Users(strong).partner = weak;
        Users(weak).partner = strong;
        p = p + 1;
    end
end

% 将未配对的用户单独一个pair
users_all = 1:P.nums;
users_noma = [Pairs(1:p-1).pair];
users_oma = setdiff(users_all, users_noma);
for i = 1:length(users_oma)
    Pairs(p).pair = [users_oma(i), 0];
    Pairs(p).strong_user = users_oma(i);
    p = p + 1;
end

strong_users = [Pairs(1:p-1).strong_user];

end

