function [ strong_users, users_noma, users_oma ] = pair_rho( rho )

global P;
global Users;
global Pairs;

% 计算信道相似度矩阵和信道增益差矩阵
Diff = zeros(P.nums, P.nums);
Corr = zeros(P.nums, P.nums);

for i = 1:P.nums
    hi = Users(i).h.';
    pli = Users(i).pathloss;
    for j = 1:P.nums
        hj = Users(j).h.';
        plj = Users(j).pathloss;
        Diff(i,j) = pli - plj;
        Corr(i,j) = norm(hi*hj') / norm(hi) / norm(hj);
    end
end

% 筛选相似度满足要求(rho)的可能配对组合
potential = zeros(P.nums,P.nums);
[Users(1:P.nums).candidates] = deal([]);
[Users(1:P.nums).partner] = deal([]);
for i = 1:P.nums
    for j = 1:P.nums
        if Corr(i,j) > rho && Diff(i,j) > 0
            Users(i).candidates = [Users(i).candidates,j];
            potential(i,j) = Diff(i,j);
        end
    end
end

% 将信道增益差最大的组合依此提取
tmp = potential;
p = 1;   % index for pairs_noma
while any(tmp(:))
    index = find(tmp == max(tmp(:)));
    [strong,weak] = ind2sub([P.nums,P.nums], index);
    Pairs(p).pair = [strong, weak];
    Pairs(p).strong_user = strong;
    Users(strong).partner = weak;
    Users(weak).partner = strong;
    p = p + 1;
    tmp(strong,:) = 0;
    tmp(weak,:) = 0;
    tmp(:,strong) = 0;
    tmp(:,weak) = 0;
end

% 将剩下没有配对的用户单独一组
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