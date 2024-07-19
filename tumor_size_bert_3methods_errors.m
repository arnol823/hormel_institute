% Matrices for Kuznetsov Tumor Modeling for Classic von Bertalanffy Model
% Approximate alpha and beta using the 3 methods and get errors
% Matrices for Kuznetsov Tumor Modeling for Classic von Bertalanffy Model
% Approximate alpha and beta using the Euler method

% Load the dataset
data = readtable('tumor_size_example_study1_pt2.csv'); 

% Extract the time vector and population data as numerical arrays
t_data = data{:, 1}; % Time vector
V = data{:, 2}; % Population data V

% Calculate time step differences
h = diff(t_data);

%% Population Data
num_points = length(V) - 1;

% Initialize matrices for Forward Difference
V_matrixA_fwd = zeros(num_points, 2);
V_matrixB_fwd = zeros(num_points, 1);

% Populate the matrices for Forward Difference
for i = 1:num_points
    V_matrixA_fwd(i, :) = [h(i) * V(i)^(2/3), -h(i) * V(i)]; 
    V_matrixB_fwd(i) = V(i+1) - V(i); 
end

% Solve the linear system using non-negative least squares for Forward Difference
A = V_matrixA_fwd;
b = V_matrixB_fwd;
xhat = lsqnonneg(A, b);  % Use lsqnonneg to ensure non-negative constraints
alpha_fwd = xhat(1)
beta_fwd = xhat(2)

% Define the differential equation for Forward Difference
dVdt_fwd = @(t, V) alpha_fwd * V^(2/3) - beta_fwd * V;

% Use ode45 to solve the differential equation for Forward Difference
[t_ode_fwd, V_ode_fwd] = ode45(dVdt_fwd, t_data, V(1));

% Calculate errors between the predicted points and actual data points for Forward Difference
predicted_V_fwd = interp1(t_ode_fwd, V_ode_fwd, t_data); % Interpolate to find predicted values at actual time points

% Error measures for Forward Difference
absolute_error_fwd = abs(predicted_V_fwd - V);
squared_error_fwd = (predicted_V_fwd - V).^2;
mean_absolute_error_fwd = mean(absolute_error_fwd);
mean_squared_error_fwd = mean(squared_error_fwd);
root_mean_squared_error_fwd = sqrt(mean_squared_error_fwd);

% Display error measures for Forward Difference
fprintf('Forward Difference Method:\n');
fprintf('Mean Absolute Error: %.4f\n', mean_absolute_error_fwd);
fprintf('Mean Squared Error: %.4f\n', mean_squared_error_fwd);
fprintf('Root Mean Squared Error: %.4f\n\n', root_mean_squared_error_fwd);

%% Initialize matrices for Backward Difference
V_matrixA_bwd = zeros(num_points, 2);
V_matrixB_bwd = zeros(num_points, 1);

% Populate the matrices for Backward Difference
for i = 2:num_points+1
    V_matrixA_bwd(i-1, :) = [h(i-1) * V(i-1)^(2/3), -h(i-1) * V(i-1)]; 
    V_matrixB_bwd(i-1) = V(i) - V(i-1); 
end

% Solve the linear system using non-negative least squares for Backward Difference
A = V_matrixA_bwd;
b = V_matrixB_bwd;
xhat = lsqnonneg(A, b);  % Use lsqnonneg to ensure non-negative constraints
alpha_bwd = xhat(1)
beta_bwd = xhat(2)

% Define the differential equation for Backward Difference
dVdt_bwd = @(t, V) alpha_bwd * V^(2/3) - beta_bwd * V;

% Use ode45 to solve the differential equation for Backward Difference
[t_ode_bwd, V_ode_bwd] = ode45(dVdt_bwd, t_data, V(1));

% Calculate errors between the predicted points and actual data points for Backward Difference
predicted_V_bwd = interp1(t_ode_bwd, V_ode_bwd, t_data); % Interpolate to find predicted values at actual time points

% Error measures for Backward Difference
absolute_error_bwd = abs(predicted_V_bwd - V);
squared_error_bwd = (predicted_V_bwd - V).^2;
mean_absolute_error_bwd = mean(absolute_error_bwd);
mean_squared_error_bwd = mean(squared_error_bwd);
root_mean_squared_error_bwd = sqrt(mean_squared_error_bwd);

% Display error measures for Backward Difference
fprintf('Backward Difference Method:\n');
fprintf('Mean Absolute Error: %.4f\n', mean_absolute_error_bwd);
fprintf('Mean Squared Error: %.4f\n', mean_squared_error_bwd);
fprintf('Root Mean Squared Error: %.4f\n\n', root_mean_squared_error_bwd);

%% Initialize matrices for Central Difference
num_points_central = length(V) - 2;
V_matrixA_central = zeros(num_points_central, 2);
V_matrixB_central = zeros(num_points_central, 1);

% Populate the matrices for Central Difference
for i = 2:num_points_central+1
    V_matrixA_central(i-1, :) = [h(i) * (V(i+1)^(2/3) + V(i-1)^(2/3)) / 2, -h(i) * (V(i+1) + V(i-1)) / 2]; 
    V_matrixB_central(i-1) = V(i+1) - V(i-1); 
end

% Solve the linear system using non-negative least squares for Central Difference
A = V_matrixA_central;
b = V_matrixB_central;
xhat = lsqnonneg(A, b);  % Use lsqnonneg to ensure non-negative constraints
alpha_central = xhat(1)
beta_central = xhat(2)

% Define the differential equation for Central Difference
dVdt_central = @(t, V) alpha_central * V^(2/3) - beta_central * V;

% Use ode45 to solve the differential equation for Central Difference
[t_ode_central, V_ode_central] = ode45(dVdt_central, t_data, V(1));

% Calculate errors between the predicted points and actual data points for Central Difference
predicted_V_central = interp1(t_ode_central, V_ode_central, t_data); % Interpolate to find predicted values at actual time points

% Error measures for Central Difference
absolute_error_central = abs(predicted_V_central - V);
squared_error_central = (predicted_V_central - V).^2;
mean_absolute_error_central = mean(absolute_error_central);
mean_squared_error_central = mean(squared_error_central);
root_mean_squared_error_central = sqrt(mean_squared_error_central);

% Display error measures for Central Difference
fprintf('Central Difference Method:\n');
fprintf('Mean Absolute Error: %.4f\n', mean_absolute_error_central);
fprintf('Mean Squared Error: %.4f\n', mean_squared_error_central);
fprintf('Root Mean Squared Error: %.4f\n\n', root_mean_squared_error_central);

% Plot the results for all methods
figure;
hold on;
plot(t_data, predicted_V_fwd, '-b', 'DisplayName', 'Predicted V (Forward Difference)');
plot(t_data, predicted_V_bwd, '-g', 'DisplayName', 'Predicted V (Backward Difference)');
plot(t_data, predicted_V_central, '-m', 'DisplayName', 'Predicted V (Central Difference)');
plot(t_data, V, 'or', 'MarkerSize', 3, 'DisplayName', 'Actual V Data Points');
hold off;

%title('Solution of Differential Equation using ode45 with Different Finite Difference Methods');
xlabel('Time', 'FontSize', 24);
ylabel('Tumor Volume', 'FontSize', 24);
legend('show');

set(legend, 'FontSize', 16, 'FontWeight', 'bold');
