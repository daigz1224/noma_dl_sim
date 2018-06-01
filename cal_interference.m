function [ I ] = cal_interference( u, p )

global P;
global Users;
global Pairs;

h = Users(u).h.';
pl = Users(u).pathloss;
I = 0;
for k = 1:length(Pairs)
    if Pairs(k).pair == 0
        break;
    end
    if k == p
        continue
    end
    w = Users(Pairs(k).pair(1)).w;
    tmp = (pl*norm(h*w))^2*P.tx_power;
    I = I + tmp;
end

end

