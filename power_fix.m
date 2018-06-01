function [] = power_fix()

global P;
global Users;
global Pairs;

% 固定功率分配
for p = 1:length(Pairs)
    if Pairs(p).pair == 0
        break;
    end
    u1 = Pairs(p).pair(1);  % strong user
    u2 = Pairs(p).pair(2);  % weak user
    if u2 == 0  % oma
        Users(u1).a = 1;
    else  % noma
        Users(u1).a = P.alpha;
        Users(u2).a = 1 - P.alpha;
    end
end

end