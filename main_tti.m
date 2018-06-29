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
P.Nt               = 8;   % 基站天线数
P.nums             = 30;  % 用户数
% drop_mode
P.drop_mode        = 1;    % | 0 - rand | 1 - pcp |
% move_mode
P.move_mode        = 0;    % | 0 - still | 1 - move |
P.speed            = 3;    % km/h for move
% pair_mode
P.pair_mode        = 2;    % | 0 - none | 1 - random | 2 - rho | 3 - kmeans |
P.rho              = 0.99;                % for pair_rho
P.K                = floor(sqrt(P.nums));  % for pair_kmeans
% schedule_mode
P.schedule_mode    = 0;    % | 0 - round | 1 - PF  |
% power_mode
P.power_mode       = 0;    % | 0 - fix | 1 - fair |
P.alpha            = 0.2;   % 功率分配因子 (0, 0.5) for strong user

P.ber              = 0.00001;
P.time             = 300;  % 仿真时间 ms

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

sum_rate = deal(zeros(1, P.time));

%% 3.撒点
switch P.drop_mode
    case 0
        drop_rand();
    case 1
        drop_pcp();
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

%% 5.功率分配

switch P.power_mode
    case 0
        power_fix();
    case 1
        power_fair();
end

for tti = 1:1:P.time
    
    %% 6.调度
    switch P.schedule_mode
        case 0
            schedule_pairs = schedule_round(tti);
        case 1
            schedule_pairs = schedule_PF(tti);
    end
    
    %% 7.波束赋形
    zfbf_users = zeros(1, length(schedule_pairs));
    
    for i = 1:1:length(schedule_pairs)
        p = schedule_pairs(i);
        strong = Pairs(p).pair(1);
        zfbf_users(i) = strong;
    end
    
    H = [Users(zfbf_users).h].';  % #strong_users 个波束赋形向量 #strong_users x Nt
    W = pinv(H);                  % hm*wn/|hm| = 0 if m≠n else 1
    
    for i = 1:1:length(schedule_pairs)
        strong = Pairs(schedule_pairs(i)).pair(1);
        Users(strong).w = W(:,i);
    end
    
    % 测试波束赋形向量的有效性
    for u = 1:P.nums
        h = (Users(u).h).';  % Nr x Nt
        Users(u).test_hW = (h/norm(h))*W;
    end
    
    %% 8.计算速率
    
    for p = schedule_pairs
        u1 = Pairs(p).pair(1);  % strong user
        u2 = Pairs(p).pair(2);  % weak user
        
        if u2 == 0  %% oma
            h  = Users(u1).h.';  % Nr x Nt
            pl = Users(u1).pathloss;
            w   = Users(u1).w;
            a  = Users(u1).a;
            
            I   = cal_interference(u, p);
            
            Gamma = pl*(norm(h*w))^2*a*P.tx_power / ...
                (P.noise_power + I);
            
            SINR = 10*log10(Gamma);
            Users(u1).SINR(tti) = SINR;
            Users(u1).rate = 0.5*P.sys_bandwidth*log2(1+Gamma);
            
            sum_rate(tti) = Users(u1).rate;
        else  %% noma
            h1  = Users(u1).h.';  % Nr x Nt
            pl1 = Users(u1).pathloss;
            w   = Users(u1).w;
            a1  = Users(u1).a;
            h2  = Users(u2).h.';  % Nr x Nt
            pl2 = Users(u2).pathloss;
            a2  = Users(u2).a;
            
            % 计算强用户的速率
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
            
            Users(u1).SINR(tti) = SINR1_x1_sic;
            Users(u1).rate(tti) = P.sys_bandwidth*((1-P.ber)*log2(1+Gamma1_x1_sic) + ...
                P.ber*log2(1+Gamma1_x1_err));
            
            % 计算弱用户的速率
            I2 = cal_interference(u2, p);
            
            Gamma2_x2 = pl2*norm(h2*w)^2*a2*P.tx_power / ...
                (P.noise_power + pl2*norm(h2*w)^2*a1*P.tx_power + I2);
            
            SINR2 = 10*log10(Gamma2_x2);  % dB
            
            Users(u2).SINR(tti) = SINR2;
            Users(u2).rate(tti) = P.sys_bandwidth*log2(1+Gamma2_x2);
            
            sum_rate(tti) = sum_rate(tti) + Users(u1).rate(tti) + Users(u2).rate(tti);  % 每个时刻的和速率
        end
    end  % end for p = schedule_pairs
end  % end for tti = 1:1:P.time

fprintf('Average Sum Rate = %f M/Hz/s\n', sum(sum_rate) / P.time / 10^6);

%% 9.画图

% 撒点与配对图
figure; hold on;
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
axis off;

% figure;
% scatter(1:1:P.time, sum_rate, '.k');
% 
% figure;
% r_vec = zeros(1, length(Users));
% for i = 1:1:length(r_vec)
%     r_vec(i) = sum(Users(i).rate);
% end
% scatter(1:1:length(Users), r_vec, '*k');