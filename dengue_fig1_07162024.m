% Dengue Control Strategy using Ross-Macdonald Model
% Reproducing Figure 1: Nullclines

clear;
clc;
close all;

% Parameters
a = 0.36065; % Mosquito biting rate on humans
mu = 1/30; % Natural mortality rate of the mosquito
p = 0.22687; % Probability of infection of a susceptible human
q = 0.08058; % Probability of infection of a susceptible mosquito
alpha = 1/10; % Recovery rate of humans
m = 2.59691; % Number of female mosquitoes per person adjsuted from 1.5969

% Reproduction number R0
R0 = (p * q * a^2 * m) / (alpha * mu);
disp(['Reproduction Number R0: ', num2str(R0)]);

% Nullclines
x = linspace(0, 1, 1000);
y1_nullcline = mu * x ./ (p * a * (1 - x));
y2_nullcline = q * a * m * x ./ (q * a * m * x + alpha);

% Equilibrium points
xE = (p * q * a^2 * m - alpha * mu) / (p * q * a^2 * m + q * a * mu * m);
yE = (p * q * a^2 * m - alpha * mu) / (p * q * a^2 * m + p * a * alpha);

% Plotting the nullclines
figure;
plot(x, y1_nullcline, 'k', 'LineWidth', 2); % eq. 21
hold on;
plot(x, y2_nullcline, 'r', 'LineWidth', 2); % eq. 22
plot(xE, yE, 'ko', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('Proportion of Infected Vectors');
ylabel('Proportion of Infected Humans');
title(['Reproduction Number R0=', num2str(R0)]);
legend('eq. 21', 'eq. 22');
grid on;

% Set axes limits for better visualization
xlim([0 1]);
ylim([0 1]);

% Add annotations to indicate equilibrium points
annotation('textarrow', [0.3 0.37], [0.35 0.3], 'String', 'P_E');
annotation('textarrow', [0.05 0.1], [0.05 0.1], 'String', 'P_0');
