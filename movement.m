function [] = movement()

global P;
global Users;

speed_m_s = P.speed * 1000 / 60 / 60 / 1000;

for u = 1:P.nums
    curr_coor = Users(u).coor;
    
    Users(u).ang = a;
    Users(u).h = h.';  % Nt x Nr
    Users(u).dist = dist;
    Users(u).coor = dist*exp(1i*a);  % зјБъ
    Users(u).pathloss = 10^(-0.1*pathloss);  % dB2lin
end


end

