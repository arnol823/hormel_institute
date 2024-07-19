% Combined script to apply all three MATLAB methods and three finite difference methods
% Load the dataset
data = readtable('combined_filtered_tumor_data_with_volume.csv'); 

% Extract the unique patient IDs
patientIDs = unique(data.patientID);

% Initialize the results table
results = table('Size', [0, 37], 'VariableTypes', repmat("double", 1, 37), ...
                'VariableNames', {'PatientID', ...
                                  'Alpha_LSQ', 'Beta_LSQ', 'MAE_LSQ', 'MSE_LSQ', 'RMSE_LSQ', 'RunTime_LSQ', ...
                                  'Alpha_FMIN', 'Beta_FMIN', 'MAE_FMIN', 'MSE_FMIN', 'RMSE_FMIN', 'RunTime_FMIN', ...
                                  'Alpha_FITNLM', 'Beta_FITNLM', 'MAE_FITNLM', 'MSE_FITNLM', 'RMSE_FITNLM', 'RunTime_FITNLM', ...
                                  'Alpha_FWD', 'Beta_FWD', 'MAE_FWD', 'MSE_FWD', 'RMSE_FWD', 'RunTime_FWD', ...
                                  'Alpha_BWD', 'Beta_BWD', 'MAE_BWD', 'MSE_BWD', 'RMSE_BWD', 'RunTime_BWD', ...
                                  'Alpha_CEN', 'Beta_CEN', 'MAE_CEN', 'MSE_CEN', 'RMSE_CEN', 'RunTime_CEN'});

% Loop over each patient ID
for pid = 1:length(patientIDs)
    patientID = patientIDs(pid);
    patient_data = data(data.patientID == patientID, :);
    
    % Extract the time vector and population data as numerical arrays
    t_data = patient_data{:, 'Treatment_Day'}; % Time vector
    V = patient_data{:, 'TargetLesionLongDiam_mm'}; % Population data V
    
    % Ensure the time vector is strictly increasing
    [t_data, sortIdx] = sort(t_data);
    V = V(sortIdx);
    
    % Calculate time step differences
    h = diff(t_data);

    %% LSQCURVEFIT Method
    tic;
    lsq_model = @(params, V) params(1) * V.^(2/3) - params(2) * V;
    params0 = [0.1, 0.1]; % Initial guess for [alpha, beta]
    options = optimoptions('lsqcurvefit', 'Display', 'off');
    [params_lsq, ~, ~] = lsqcurvefit(lsq_model, params0, V, V, [], [], options);
    run_time_lsq = toc;
    alpha_lsq = params_lsq(1);
    beta_lsq = params_lsq(2);
    
    predicted_V_lsq = lsq_model(params_lsq, V);
    absolute_error_lsq = abs(predicted_V_lsq - V);
    squared_error_lsq = (predicted_V_lsq - V).^2;
    mean_absolute_error_lsq = mean(absolute_error_lsq);
    mean_squared_error_lsq = mean(squared_error_lsq);
    root_mean_squared_error_lsq = sqrt(mean_squared_error_lsq);
    
    %% FMINCON Method
    tic;
    fmin_model = @(params) sum((params(1) * V.^(2/3) - params(2) * V - V).^2);
    options = optimoptions('fmincon', 'Display', 'off');
    params_fmin = fmincon(fmin_model, params0, [], [], [], [], [], [], [], options);
    run_time_fmin = toc;
    alpha_fmin = params_fmin(1);
    beta_fmin = params_fmin(2);
    
    predicted_V_fmin = alpha_fmin * V.^(2/3) - beta_fmin * V;
    absolute_error_fmin = abs(predicted_V_fmin - V);
    squared_error_fmin = (predicted_V_fmin - V).^2;
    mean_absolute_error_fmin = mean(absolute_error_fmin);
    mean_squared_error_fmin = mean(squared_error_fmin);
    root_mean_squared_error_fmin = sqrt(mean_squared_error_fmin);
    
    %% FITNLM Method
    tic;
    fitnlm_model = @(b, V) b(1) * V.^(2/3) - b(2) * V;
    mdl = fitnlm(V, V, fitnlm_model, params0);
    run_time_fitnlm = toc;
    params_fitnlm = mdl.Coefficients.Estimate';
    alpha_fitnlm = params_fitnlm(1);
    beta_fitnlm = params_fitnlm(2);
    
    predicted_V_fitnlm = fitnlm_model(params_fitnlm, V);
    absolute_error_fitnlm = abs(predicted_V_fitnlm - V);
    squared_error_fitnlm = (predicted_V_fitnlm - V).^2;
    mean_absolute_error_fitnlm = mean(absolute_error_fitnlm);
    mean_squared_error_fitnlm = mean(squared_error_fitnlm);
    root_mean_squared_error_fitnlm = sqrt(mean_squared_error_fitnlm);
    
    %% Forward Difference Method
    tic;
    num_points = length(V) - 1;
    V_matrixA_fwd = zeros(num_points, 2);
    V_matrixB_fwd = zeros(num_points, 1);

    for i = 1:num_points
        V_matrixA_fwd(i, :) = [h(i) * V(i)^(2/3), -h(i) * V(i)]; 
        V_matrixB_fwd(i) = V(i+1) - V(i); 
    end

    A = V_matrixA_fwd;
    b = V_matrixB_fwd;
    xhat = lsqnonneg(A, b);
    run_time_fwd = toc;
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
    tic;
    V_matrixA_bwd = zeros(num_points, 2);
    V_matrixB_bwd = zeros(num_points, 1);

    for i = 2:num_points+1
        V_matrixA_bwd(i-1, :) = [h(i-1) * V(i-1)^(2/3), -h(i-1) * V(i-1)]; 
        V_matrixB_bwd(i-1) = V(i) - V(i-1); 
    end

    A = V_matrixA_bwd;
    b = V_matrixB_bwd;
    xhat = lsqnonneg(A, b);
    run_time_bwd = toc;
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
    tic;
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
    run_time_central = toc;
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
    results = [results; {patientID, ...
                         alpha_lsq, beta_lsq, mean_absolute_error_lsq, mean_squared_error_lsq, root_mean_squared_error_lsq, run_time_lsq, ...
                         alpha_fmin, beta_fmin, mean_absolute_error_fmin, mean_squared_error_fmin, root_mean_squared_error_fmin, run_time_fmin, ...
                         alpha_fitnlm, beta_fitnlm, mean_absolute_error_fitnlm, mean_squared_error_fitnlm, root_mean_squared_error_fitnlm, run_time_fitnlm, ...
                         alpha_fwd, beta_fwd, mean_absolute_error_fwd, mean_squared_error_fwd, root_mean_squared_error_fwd, run_time_fwd, ...
                         alpha_bwd, beta_bwd, mean_absolute_error_bwd, mean_squared_error_bwd, root_mean_squared_error_bwd, run_time_bwd, ...
                         alpha_central, beta_central, mean_absolute_error_central, mean_squared_error_central, root_mean_squared_error_central, run_time_central}];

    % Plot the results for all methods
    figure;
    hold on;
    plot(t_data, predicted_V_lsq, '-k', 'DisplayName', 'Predicted V (LSQ)');
    plot(t_data, predicted_V_fmin, '-r', 'DisplayName', 'Predicted V (FMINCON)');
    plot(t_data, predicted_V_fitnlm, '-c', 'DisplayName', 'Predicted V (FITNLM)');
    plot(t_data, predicted_V_fwd, '-b', 'DisplayName', 'Predicted V (Forward Difference)');
    plot(t_data, predicted_V_bwd, '-g', 'DisplayName', 'Predicted V (Backward Difference)');
    plot(t_data, predicted_V_central, '-m', 'DisplayName', 'Predicted V (Central Difference)');
    plot(t_data, V, 'or', 'MarkerSize', 3, 'DisplayName', 'Actual V Data Points');
    set(legend, 'FontSize', 16, 'FontWeight', 'bold');
    hold off;

    title(['Solution of Differential Equation using Various Methods for Patient ' num2str(patientID)]);
    xlabel('Time', 'FontSize', 24);
    ylabel('Tumor Volume', 'FontSize', 24);
    legend('show');
end

% Save the results table to a file
writetable(results, 'combined_results.csv');

% Display the results table
disp(results);
