% COVID_poster plotting

%% First Younger Population (Population 1)
% Initial conditions: [S_1, I_1, R_1, D_1, I_2]
Y01 = [S_1(1); I_1(1); R_1(1); D_1(1); I_2(1)]; % Example initial conditions, assuming I_1 is initially 10

% Time span for the solution
tspan = [0 100];

% Define the function handle including the parameters
odeFunc1 = @(t, Y) myODESystem1(t, Y, mu, N_1, beta_11, beta_12, alpha1_plus_kappa1, kappa_1);

% Solve the ODE
[t1, Y1] = ode45(odeFunc1, tspan, Y01);

% Plot the results for Population 1
figure;
subplot(2, 2, 1);
hold on;
plot(t1, Y1(:,1), '-', 'DisplayName', 'S_1 (ODE)');
plot(S_1, '-o', 'DisplayName', 'S_1');
title('S_1 Population');
xlabel('Time');
ylabel('Population');
legend show;
hold off;

subplot(2, 2, 2);
hold on;
plot(t1, Y1(:,2), '-', 'DisplayName', 'I_1 (ODE)');
plot(I_1, '-o', 'DisplayName', 'I_1');
title('I_1 Population');
xlabel('Time');
ylabel('Population');
legend show;
hold off;

subplot(2, 2, 3);
hold on;
plot(t1, Y1(:,3), '-', 'DisplayName', 'R_1 (ODE)');
plot(R_1, '-o', 'DisplayName', 'R_1');
title('R_1 Population');
xlabel('Time');
ylabel('Population');
legend show;
hold off;

subplot(2, 2, 4);
hold on;
plot(t1, Y1(:,4), '-', 'DisplayName', 'D_1 (ODE)');
plot(D_1, '-o', 'DisplayName', 'D_1');
title('D_1 Population');
xlabel('Time');
ylabel('Population');
legend show;
hold off;

% Save as myODESystem1.m
function dYdt = myODESystem1(t, Y, mu, N_1, beta_11, beta_12, alpha1_plus_kappa1, kappa_1)
    S_1 = Y(1);
    I_1 = Y(2);
    R_1 = Y(3);
    D_1 = Y(4);
    I_2 = Y(5); % Assuming I_2 is provided as part of the state vector
    
    dS_1dt = mu * N_1 - beta_11 * (S_1 * I_1 / N_1) - beta_12 * (S_1 * I_2 / N_1) - mu * S_1;
    dI_1dt = beta_11 * (S_1 * I_1 / N_1) + beta_12 * (S_1 * I_2 / N_1) - (mu + alpha1_plus_kappa1) * I_1;
    dR_1dt = alpha1_plus_kappa1 * I_1 - mu * R_1;
    dD_1dt = kappa_1 * I_1;
    
    dYdt = [dS_1dt; dI_1dt; dR_1dt; dD_1dt; 1];
end

%% Second Older Population (Population 2)
% Initial conditions: [S_2, I_2, R_2, D_2, I_1]
Y02 = [S_2(1); I_2(1); R_2(1); D_2(1); I_1(1)]; % Example initial conditions, assuming I_1 is initially 10

% Time span for the solution
tspan = [0 100];

% Define the function handle including the parameters
odeFunc2 = @(t, Y) myODESystem2(t, Y, mu, N_2, beta_22, beta_21, alpha2_plus_kappa2, kappa_2);

% Solve the ODE
[t2, Y2] = ode45(odeFunc2, tspan, Y02);

% Plot the results for Population 2
figure;
subplot(2, 2, 1);
hold on;
plot(t2, Y2(:,1), '-', 'DisplayName', 'S_2 (ODE)');
plot(S_2, '-o', 'DisplayName', 'S_2');
title('S_2 Population');
xlabel('Time');
ylabel('Population');
legend show;
hold off;

subplot(2, 2, 2);
hold on;
plot(t2, Y2(:,2), '-', 'DisplayName', 'I_2 (ODE)');
plot(I_2, '-o', 'DisplayName', 'I_2');
title('I_2 Population');
xlabel('Time');
ylabel('Population');
legend show;
hold off;

subplot(2, 2, 3);
hold on;
plot(t2, Y2(:,3), '-', 'DisplayName', 'R_2 (ODE)');
plot(R_2, '-o', 'DisplayName', 'R_2');
title('R_2 Population');
xlabel('Time');
ylabel('Population');
legend show;
hold off;

subplot(2, 2, 4);
hold on;
plot(t2, Y2(:,4), '-', 'DisplayName', 'D_2 (ODE)');
plot(D_2, '-o', 'DisplayName', 'D_2');
title('D_2 Population');
xlabel('Time');
ylabel('Population');
legend show;
hold off;

% Save as myODESystem2.m
function dYdt = myODESystem2(t, Y, mu, N_2, beta_22, beta_21, alpha2_plus_kappa2, kappa_2)
    S_2 = Y(1);
    I_2 = Y(2);
    R_2 = Y(3);
    D_2 = Y(4);
    I_1 = Y(5); % Assuming I_1 is provided as part of the state vector
    
    dS_2dt = mu * N_2 - beta_22 * (S_2 * I_2 / N_2) - beta_21 * (S_2 * I_1 / N_2) - mu * S_2;
    dI_2dt = beta_22 * (S_2 * I_2 / N_2) + beta_21 * (S_2 * I_1 / N_2) - (mu + alpha2_plus_kappa2) * I_2;
    dR_2dt = alpha2_plus_kappa2 * I_2 - mu * R_2;
    dD_2dt = kappa_2 * I_2;
    
    dYdt = [dS_2dt; dI_2dt; dR_2dt; dD_2dt; 1];
end

%% Plotting both populations in subplots
figure;
subplot(2, 2, 1);
hold on;
plot(t1, Y1(:,1), '-o', 'DisplayName', 'S_1 (ODE)');
plot(S_1, '-x', 'DisplayName', 'S_1');
plot(t2, Y2(:,1), '-x', 'DisplayName', 'S_2 (ODE)');
plot(S_2, '-s', 'DisplayName', 'S_2');
hold off;
title('S_1 and S_2 Population');
xlabel('Time');
ylabel('Population');
legend('show');

subplot(2, 2, 2);
hold on;
plot(t1, Y1(:,2), '-o', 'DisplayName', 'I_1 (ODE)');
plot(I_1, '-x', 'DisplayName', 'I_1');
plot(t2, Y2(:,2), '-x', 'DisplayName', 'I_2 (ODE)');
plot(I_2, '-s', 'DisplayName', 'I_2');
hold off;
title('I_1 and I_2 Population');
xlabel('Time');
ylabel('Population');
legend('show');

subplot(2, 2, 3);
hold on;
plot(t1, Y1(:,3), '-o', 'DisplayName', 'R_1 (ODE)');
plot(R_1, '-x', 'DisplayName', 'R_1');
plot(t2, Y2(:,3), '-x', 'DisplayName', 'R_2 (ODE)');
plot(R_2, '-s', 'DisplayName', 'R_2');
hold off;
title('R_1 and R_2 Population');
xlabel('Time');
ylabel('Population');
legend('show');

subplot(2, 2, 4);
hold on;
plot(t1, Y1(:,4), '-o', 'DisplayName', 'D_1 (ODE)');
plot(D_1, '-x', 'DisplayName', 'D_1');
plot(t2, Y2(:,4), '-x', 'DisplayName', 'D_2 (ODE)');
plot(D_2, '-s', 'DisplayName', 'D_2');
hold off;
title('D_1 and D_2 Population');
xlabel('Time');
ylabel('Population');
legend('show');

%% Plotting all populations on a single plot
figure;
hold on;
plot(t1, Y1(:,1), '-o', 'DisplayName', 'S_1 (ODE)');
plot(S_1, '-x', 'DisplayName', 'S_1');
plot(t2, Y2(:,1), '-x', 'DisplayName', 'S_2 (ODE)');
plot(S_2, '-s', 'DisplayName', 'S_2');
plot(t1, Y1(:,2), '-o', 'DisplayName', 'I_1 (ODE)');
plot(I_1, '-x', 'DisplayName', 'I_1');
plot(t2, Y2(:,2), '-x', 'DisplayName', 'I_2 (ODE)');
plot(I_2, '-s', 'DisplayName', 'I_2');
plot(t1, Y1(:,3), '-o', 'DisplayName', 'R_1 (ODE)');
plot(R_1, '-x', 'DisplayName', 'R_1');
plot(t2, Y2(:,3), '-x', 'DisplayName', 'R_2 (ODE)');
plot(R_2, '-s', 'DisplayName', 'R_2');
plot(t1, Y1(:,4), '-o', 'DisplayName', 'D_1 (ODE)');
plot(D_1, '-x', 'DisplayName', 'D_1');
plot(t2, Y2(:,4), '-x', 'DisplayName', 'D_2 (ODE)');
plot(D_2, '-s', 'DisplayName', 'D_2');
hold off;
title('Population Over Time');
xlabel('Time');
ylabel('Population');
legend('show');
