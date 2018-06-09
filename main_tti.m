% Title: 不同分簇与功率分配算法仿真对比
clear -global; clearvars; close all; clc;
dbstop if error;
digits(64);  % 64位运算精度，默认32

global P;
global Users;
global Pairs;

%% 1.设置参数

% 数值属性
P.cell_radius      = 500;        % m
P.sys_bandwidth    = 4.32*10^6;  % Hz
P.tx_power_dB      = 43;         % dBm
P.noise_density_dB = -169;       % dBm/Hz
P.Nr               = 1;
P.Nt               = 128;
P.nums             = 60;
P.pair_mode        = 2;  % 0 - none 1 - random 2 - rho 3 - kmeans
P.rho              = 0.95;  % for pair_rho
P.K                = ceil(0.4 * P.nums);  % for pair_kmeans
P.alpha = 0.2;  % 功率分配因子 (0, 0.5) for strong user
P.ber   = 0.001;
P.time   = 100;  % ms

% 计算部分参数
P.tx_power      = 10^(P.tx_power_dB / 10);  % dB2lin mW
P.noise_power   = P.sys_bandwidth*10^(0.1*P.noise_density_dB);  % dB2lin mW

% 检查参数是否符合要求
% assert(P.Nt >= P.nums*P.Nr);

%% 2.初始化

rng(987654321);

[Users(1:P.nums).ang]        = deal(0);
[Users(1:P.nums).h]          = deal(0);
[Users(1:P.nums).dist]       = deal(0);
[Users(1:P.nums).coor]       = deal(0);
[Users(1:P.nums).pathloss]   = deal(0);
[Users(1:P.nums).w]          = deal(0);
[Users(1:P.nums).test_hW]    = deal(0);
[Users(1:P.nums).candidates] = deal(0);
[Users(1:P.nums).partner]    = deal(0);
[Users(1:P.nums).rate]       = deal(zeros(1, P.time));

sum_rate = deal(0);

%% 3. 撒点
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

%% 4.分簇配对
switch P.pair_mode
    case 0
        [ strong_users, users_noma, users_oma ] = pair_none();
    case 1
        [ strong_users, users_noma, users_oma ] = pair_random();
    case 2
        [ strong_users, users_noma, users_oma ] = pair_rho(P.rho);
    case 3
        [ strong_users, users_noma, users_oma ] = pair_kmeans(P.K);
end

for tti = 1:1:P.time
    
    %% 7.计算本次调度的配对用户
    
    schedule_pairs = cal_priority();
    
    %% 5.波束赋形
    zfbf_users = zeors(1, length(schedule_pairs));
    
    for i = 1:1:length(schedule_pairs)
        p = schedule_pairs(i);
        strong = Pairs(p).pair(1);
        zfbf_users(i) = strong;
    end
    
    H = [Users(zfbf_users).h].';  % #strong_users 个波束赋形向量 #strong_users x Nt
    W = pinv(H);                   % hm*wn/|hm| = 0 if m≠n else 1
    for p = schedule_pairs
        strong = Pairs(p).pair(1);  % strong user
        Users(strong).w = W(:,p);
    end
    
    % 测试波束赋形向量的有效性
    for u = 1:P.nums
        h = (Users(u).h).';  % Nr x Nt
        Users(u).test_hW = (h/norm(h))*W;
    end
    
    %% 6.功率分配
    
    power_fix();
    
    %% 8.计算速率
    
    for p = schedule_pairs
        u1 = Pairs(p).pair(1);  % strong user
        u2 = Pairs(p).pair(2);  % weak user
        
        if u2 == 0  %% oma
            h1  = Users(u1).h.';  % Nr x Nt
            pl1 = Users(u1).pathloss;
            w   = Users(u1).w;
            a1  = Users(u1).a;
            I   = cal_interference(u1, p);
            Gammma = pl1*(norm(h1*w))^2*a1*P.tx_power / ...
                (P.noise_power + I);
            SINR = 10*log10(Gammma);
            Users(u1).rate = 0.5*P.sys_bandwidth*log2(1+Gammma);
            sum_rate(p) = Users(u1).rate;
        else  %% noma
            h1  = Users(u1).h.';  % Nr x Nt
            pl1 = Users(u1).pathloss;
            w   = Users(u1).w;
            a1  = Users(u1).a;
            h2  = Users(u2).h.';  % Nr x Nt
            pl2 = Users(u2).pathloss;
            a2  = Users(u2).a;
            
            % cal strong user's rate
            I1 = cal_interference(u1, p);
            Gamma1_x2 = (pl1*norm(h1*w)^2*a2*P.tx_power) / ...
                (P.noise_power + pl1*norm(h1*w)^2*a1*P.tx_power + I1);
            SINR1_x2 = 10*log10(Gamma1_x2);  % dB
            Gamma1_x1_sic = (pl1*norm(h1*w)^2*a1*P.tx_power) / ...
                (P.noise_power + I1);
            SINR1_x1_sic = 10*log10(Gamma1_x1_sic);  % dB
            Gamma1_x1_err = (pl1*norm(h1*w)^2*a1*P.tx_power) / ...
                (P.noise_power + pl1*norm(h1*w)^2*a2*P.tx_power + I1);
            SINR1_x1_err = 10*log10(Gamma1_x1_err);  % dB
            Users(u1).rate = P.sys_bandwidth*((1-P.ber)*log2(1+Gamma1_x1_sic) + ...
                P.ber*log2(1+Gamma1_x1_err));
            
            % cal weak user's rate
            I2 = cal_interference(u2, p);
            Gamma2_x2 = pl2*norm(h2*w)^2*a2*P.tx_power / ...
                (P.noise_power + pl2*norm(h2*w)^2*a1*P.tx_power + I2);
            SINR2 = 10*log10(Gamma2_x2);  % dB
            Users(u2).rate = P.sys_bandwidth*log2(1+Gamma2_x2);
            
            sum_rate(p) = Users(u1).rate + Users(u2).rate;
        end
    end  % end for c = 1:length(Pairs)
end


% disp('Sum Rate = M/Hz/s');
% disp(sum(sum_rate) / 10^6);

%% 9.画图

% 撒点与配对图
figure(1); hold on;
axis square;
plot(P.cell_radius*exp(1i*(pi/3*(0:6))),'-.k','linewidth',2);
plot(0,'h', 'MarkerEdgeColor', 'k',...
    'MarkerFaceColor', 'k', 'MarkerSize', 12);
colors = colormap(hsv(length(users_noma)/2));
shift = -10-10*1i;

for p = 1:length(Pairs)
    u1 = Pairs(p).pair(1);
    u2 = Pairs(p).pair(2);
    if u2 == 0
        coor = Users(u1).coor;
        plot(coor+shift, 'x', 'MarkerEdgeColor', 'k',...
            'MarkerFaceColor', 'k', 'MarkerSize', 8)
        text(real(coor),imag(coor),num2str(u1));
    else
        coor1 = Users(u1).coor;
        coor2 = Users(u2).coor;
        plot(coor1+shift, 'o', 'MarkerEdgeColor', 'k',...
            'MarkerFaceColor', colors(p,:), 'MarkerSize', 8)
        text(real(coor1),imag(coor1),num2str(u1));
        plot(coor2+shift, 'o', 'MarkerEdgeColor', 'k',...
            'MarkerFaceColor', colors(p,:), 'MarkerSize', 8)
        text(real(coor2),imag(coor2),num2str(u2));
    end
end

for u = 1:P.nums
    coor = Users(u).coor;
    text(real(coor),imag(coor),num2str(u));
end

axis off;