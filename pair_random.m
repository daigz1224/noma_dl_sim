function [ strong_users, users_noma, users_oma ] = pair_random()

global P;
global Users;
global Pairs;

p = 1;  % index for cluster_noma
for u = 1:P.nums/2
    strong = u;
    weak = P.nums - u + 1;
    Pairs(p).pair = [strong, weak];
    Pairs(p).strong_user = strong;
    Users(strong).partner = weak;
    Users(weak).partner = strong;
    p = p + 1;
end

users_all = 1:P.nums;
users_noma = [Pairs(1:p-1).pair];
users_oma = setdiff(users_all, users_noma);

strong_users = [Pairs(1:p-1).strong_user];

end

