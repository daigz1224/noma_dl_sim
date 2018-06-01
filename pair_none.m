function [ strong_users, users_noma, users_oma ] = pair_none()

global P;
global Pairs;

k = 1;  % index for cluster_noma
for u = 1:P.nums
    strong = u;
    Pairs(k).pair = [strong, 0];
    Pairs(k).strong_user = strong;
    k = k + 1;
end

users_all = 1:P.nums;
users_noma = [Pairs(1:k-1).pair];
users_oma = setdiff(users_all, users_noma);

strong_users = [Pairs(1:k-1).strong_user];

end

