clear -global; clearvars; close all; clc;

alpha_vector = 0.001:0.001:0.5;

R1 = load('R1.mat');
R2 = load('R2.mat');
R3 = load('R3.mat');
R4 = load('R4.mat');

plot(alpha_vector, sum(R1.R_vector,2), '-', 'DisplayName', 'R1');hold on;
plot(alpha_vector, sum(R2.R_vector,2), '-', 'DisplayName', 'R2');hold on;
plot(alpha_vector, sum(R3.R_vector,2), '-', 'DisplayName', 'R3');hold on;
plot(alpha_vector, sum(R4.R_vector,2), '-', 'DisplayName', 'R4');hold on;

grid on;
legend('show')