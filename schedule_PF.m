function [ schedule_pairs ] = schedule_PF(tti)

global P;
global Users;
global Pairs;

schedule_pairs = zeros(1, P.Nt);
[Pairs(1:length(Pairs)).priority] = deal(0);

for p = 1:1:length(Pairs)
    u1 = Pairs(p).pair(1);  % strong user
    u2 = Pairs(p).pair(2);  % weak user
    
    if u2 == 0
        pl = Users(u1).pathloss;
        a = Users(u1).a;
        
        % 计算瞬时速率（假设没有干扰）
        Gamma = pl * a * P.tx_power / P.noise_power;
        SINR = 10 * log10(Gamma);
        rate = P.sys_bandwidth * log2(1 + Gamma);
        
        % 统计历史速率
        ave_rate = (sum(Users(u1).rate)+10^-6) / tti;
        
        % 计算用户组的优先级
        Pairs(p).priority = rate / ave_rate;
    else
        pl1 = Users(u1).pathloss;
        pl2 = Users(u2).pathloss;
        a1 = Users(u1).a;
        a2 = Users(u2).a;
        
        % 计算瞬时速率（假设没有干扰）
        Gamma_u1 = pl1 * a1 * P.tx_power / P.noise_power;
        SINR_u1 = 10 * log10(Gamma_u1);
        rate_u1 = P.sys_bandwidth * log2(1 + Gamma_u1);
        
        Gamma_u2 = pl2 * a2 * P.tx_power / P.noise_power;
        SINR_u2 = 10 * log10(Gamma_u2);
        rate_u2 = P.sys_bandwidth * log2(1 + Gamma_u2);
        
        % 统计历史速率
        ave_rate_u1 = (sum(Users(u1).rate)+10^-6) / tti;
        ave_rate_u2 = (sum(Users(u2).rate)+10^-6) / tti;
        
        % 计算用户组的优先级
        % Pairs(p).priority = (rate_u1 + rate_u2) / (ave_rate_u1 + ave_rate_u2);
        % disp(rate_u1 / ave_rate_u1);
        % disp(rate_u2 / ave_rate_u2);
        Pairs(p).priority = rate_u1 / ave_rate_u1 + rate_u2 / ave_rate_u2;
    end
end

priority_vec = [Pairs(1:1:length(Pairs)).priority];

for i = 1:1:length(schedule_pairs)
    [~, index] = max(priority_vec);
    priority_vec(index) = 0;
    schedule_pairs(i) = index;
end

end