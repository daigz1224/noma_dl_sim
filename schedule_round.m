function [ schedule_pairs ] = schedule_round(tti)

global P;
global Pairs;

schedule_pairs = zeros(1, P.Nt);

index_start = mod((tti - 1) * P.Nt + 1, length(Pairs));
index_end = mod(tti * P.Nt, length(Pairs));

if index_start < index_end
    range_pairs = index_start:1:index_end;
else
    range_pairs = [index_start:1:length(Pairs) , 1:1:index_end];
end

for i = 1:1:length(schedule_pairs)
    schedule_pairs(i) = range_pairs(i);
end

end