function [] = drop_pcp()

global P;
global Users;

lambda = log2(P.nums-2);  % 生成M点的密度
M = 1;
U = rand();
while U >= exp(-lambda)
    U = U * rand();
    M = M + 1;
end

M_x = zeros(1, M);
M_y = zeros(1, M);
for i = 1:M
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
        if dist > 35 && dist < 400 % 避免灯下黑
            break
        end
    end
    
    coor = dist*exp(1i*a);  % 坐标
    M_x(i) = real(coor);
    M_y(i) = imag(coor);
end

% 每个M点周围会生成点的个数
tmp = rand(1, M);
tmp_norm = tmp(1:M)/sum(tmp);
num_per_M = ceil(tmp_norm * P.nums);
m_x = [];
m_y = [];
for i = 1:M
    N = num_per_M(i);
    R = 30*(sqrt(N)+1);  % 用户间距尽量大于20m
    theta = unifrnd(-1,1,[1, N]);
    x = zeros(1, N);
    y = zeros(1, N);
    for j = 1:N
        while 1
            r = R * rand();
            x(j) = M_x(i) + r * cos(theta(j));
            y(j) = M_y(i) + r * sin(theta(j));
            a = atan2(y(j), x(j));
            dist = sqrt(y(j)^2 + x(j)^2);
            if (a > pi/3 && a < 2*pi/3) || (a > 4*pi/3 && a < 5*pi/3)
                angle = a - pi/3;
            elseif (a > 2*pi/3 && a < pi) || (a > 5*pi/3 && a < 2*pi)
                angle = a - 2*pi/3;
            else
                angle = a;
            end
            border = sqrt(3)*P.cell_radius / (2*sin(pi/3 + angle));  % 六边形边界
            if dist < abs(border) && dist > 35
                break;
            end
        end
    end
    m_x = [m_x, x];
    m_y = [m_y, y];
end

m_index = 1;
for u = 1:P.nums
    a = atan2(m_y(m_index), m_x(m_index));
    dist = sqrt(m_y(m_index)^2 + m_x(m_index)^2);

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
    
    m_index = m_index + 1;
end





% figure; hold on;
% axis square;
% plot(P.cell_radius*exp(1i*(pi/3*(0:6))),'-.k','linewidth',2);
% plot(0,'h', 'MarkerEdgeColor', 'k',...
%     'MarkerFaceColor', 'k', 'MarkerSize', 12);
% 
% for i = 1:M
%    plot(M_x(i),M_y(i),'ob'); 
% end
% 
% 
% for u = 1:P.nums
%     coor = Users(u).coor;
%     plot(coor, 'x', 'MarkerEdgeColor', 'k',...
%             'MarkerFaceColor', 'k', 'MarkerSize', 8)
%     text(real(coor),imag(coor),num2str(u));
% end
% 
% axis off;

end