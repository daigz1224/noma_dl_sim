function [] = drop_rand()

global P;
global Users;

for u = 1:P.nums
    a = rand()*2*pi;  % 随机生成角度 a∈(0,2*pi)
    if (a > pi/3 && a < 2*pi/3) || (a > 4*pi/3 && a < 5*pi/3)
        angle = a - pi/3;
    elseif (a > 2*pi/3 && a < pi) || (a > 5*pi/3 && a < 2*pi)
        angle = a - 2*pi/3;
    else
        angle = a;
    end
    border = sqrt(3)*P.cell_radius / (2*sin(pi/3 + angle));  % 六边形边界
    
    while 1
        dist = rand()*abs(border);  % 距离 m
        if dist > 35  % 避免灯下黑
            break
        end
    end
    
    pathloss = 128.1 + 37.6 * log10(dist / 1000);  % 路损 dB
    h = deal(0);  % Nr x Nt
    for i = 1:P.Nt
        h(i) = exp(1i*sin(a - 2*pi*(i-1)/P.Nt));
    end
    Users(u).ang = a;
    Users(u).h = h.';  % Nt x Nr
    Users(u).dist = dist;
    Users(u).coor = dist*exp(1i*a);  % 坐标
    Users(u).pathloss = 10^(-0.1*pathloss);  % dB2lin
end

end

