% applies all three matlab methods (lsqcurvefit, fmincon, and fitnlm) for parameter estimation, calculates the error values, and stores them in a table (dataframe).
% 3 Finite difference approaches 
% use file with multiple patients
% Matrices for Kuznetsov Tumor Modeling for Classic von Bertalanffy Model
% Approximate alpha and beta using the 3 methods and get errors

% Load the dataset
data = readtable('combined_filtered_tumor_data_example01_02.csv'); 

% Extract the unique patient IDs
patientIDs = unique(data.patientID);

% Initialize the results table
results = table('Size', [0, 16], 'VariableTypes', repmat("double", 1, 16), ...
                'VariableNames', {'PatientID', 'Alpha_FWD', 'Beta_FWD', 'MAE_FWD', 'MSE_FWD', 'RMSE_FWD', ...
                                  'Alpha_BWD', 'Beta_BWD', 'MAE_BWD', 'MSE_BWD', 'RMSE_BWD', ...
                                  'Alpha_CEN', 'Beta_CEN', 'MAE_CEN', 'MSE_CEN', 'RMSE_CEN'});

% Loop over each patient ID
for pid = 1:length(patientIDs)
    patientID = patientIDs(pid);
    patient_data = data(data.patientID == patientID, :);
    
    % Extract the time vector and population data as numerical arrays
    t_data = patient_data{:, 'Treatment_Day'}; % Time vector
    V = patient_data{:, 'TargetLesionLongDiam_mm'}; % Population data V

    % Calculate time step differences
    h = diff(t_data);

    %% Population Data
    num_points = length(V) - 1;

    %% Forward Difference Method
    V_matrixA_fwd = zeros(num_points, 2);
    V_matrixB_fwd = zeros(num_points, 1);

    for i = 1:num_points
        V_matrixA_fwd(i, :) = [h(i) * V(i)^(2/3), -h(i) * V(i)]; 
        V_matrixB_fwd(i) = V(i+1) - V(i); 
    end

    A = V_matrixA_fwd;
    b = V_matrixB_fwd;
    xhat = lsqnonneg(A, b);
    alpha_fwd = xhat(1);
    beta_fwd = xhat(2);

    dVdt_fwd = @(t, V) alpha_fwd * V^(2/3) - beta_fwd * V;
    [t_ode_fwd, V_ode_fwd] = ode45(dVdt_fwd, t_data, V(1));
    predicted_V_fwd = interp1(t_ode_fwd, V_ode_fwd, t_data);

    absolute_error_fwd = abs(predicted_V_fwd - V);
    squared_error_fwd = (predicted_V_fwd - V).^2;
    mean_absolute_error_fwd = mean(absolute_error_fwd);
    mean_squared_error_fwd = mean(squared_error_fwd);
    root_mean_squared_error_fwd = sqrt(mean_squared_error_fwd);

    %% Backward Difference Method
    V_matrixA_bwd = zeros(num_points, 2);
    V_matrixB_bwd = zeros(num_points, 1);

    for i = 2:num_points+1
        V_matrixA_bwd(i-1, :) = [h(i-1) * V(i-1)^(2/3), -h(i-1) * V(i-1)]; 
        V_matrixB_bwd(i-1) = V(i) - V(i-1); 
    end

    A = V_matrixA_bwd;
    b = V_matrixB_bwd;
    xhat = lsqnonneg(A, b);
    alpha_bwd = xhat(1);
    beta_bwd = xhat(2);

    dVdt_bwd = @(t, V) alpha_bwd * V^(2/3) - beta_bwd * V;
    [t_ode_bwd, V_ode_bwd] = ode45(dVdt_bwd, t_data, V(1));
    predicted_V_bwd = interp1(t_ode_bwd, V_ode_bwd, t_data);

    absolute_error_bwd = abs(predicted_V_bwd - V);
    squared_error_bwd = (predicted_V_bwd - V).^2;
    mean_absolute_error_bwd = mean(absolute_error_bwd);
    mean_squared_error_bwd = mean(squared_error_bwd);
    root_mean_squared_error_bwd = sqrt(mean_squared_error_bwd);

    %% Central Difference Method
    num_points_central = length(V) - 2;
    V_matrixA_central = zeros(num_points_central, 2);
    V_matrixB_central = zeros(num_points_central, 1);

    for i = 2:num_points_central+1
        V_matrixA_central(i-1, :) = [h(i) * (V(i+1)^(2/3) + V(i-1)^(2/3)) / 2, -h(i) * (V(i+1) + V(i-1)) / 2]; 
        V_matrixB_central(i-1) = V(i+1) - V(i-1); 
    end

    A = V_matrixA_central;
    b = V_matrixB_central;
    xhat = lsqnonneg(A, b);
    alpha_central = xhat(1);
    beta_central = xhat(2);

    dVdt_central = @(t, V) alpha_central * V^(2/3) - beta_central * V;
    [t_ode_central, V_ode_central] = ode45(dVdt_central, t_data, V(1));
    predicted_V_central = interp1(t_ode_central, V_ode_central, t_data);

    absolute_error_central = abs(predicted_V_central - V);
    squared_error_central = (predicted_V_central - V).^2;
    mean_absolute_error_central = mean(absolute_error_central);
    mean_squared_error_central = mean(squared_error_central);
    root_mean_squared_error_central = sqrt(mean_squared_error_central);

    % Store the results in the results table
    results = [results; {patientID, alpha_fwd, beta_fwd, mean_absolute_error_fwd, mean_squared_error_fwd, root_mean_squared_error_fwd, ...
                         alpha_bwd, beta_bwd, mean_absolute_error_bwd, mean_squared_error_bwd, root_mean_squared_error_bwd, ...
                         alpha_central, beta_central, mean_absolute_error_central, mean_squared_error_central, root_mean_squared_error_central}];

    % Plot the results for all methods
    figure;
    hold on;
    plot(t_data, predicted_V_fwd, '-b', 'DisplayName', 'Predicted V (Forward Difference)');
    plot(t_data, predicted_V_bwd, '-g', 'DisplayName', 'Predicted V (Backward Difference)');
    plot(t_data, predicted_V_central, '-m', 'DisplayName', 'Predicted V (Central Difference)');
    plot(t_data, V, 'or', 'MarkerSize', 3, 'DisplayName', 'Actual V Data Points');
    hold off;

    title(['Solution of Differential Equation using ode45 for Patient ' num2str(patientID)]);
    xlabel('Time');
    ylabel('Tumor Volume');
    legend('show');
end

% Save the results table to a file
writetable(results, 'results.csv');

% Display the results table
disp(results);
