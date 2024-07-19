% generates table with RMSE and run times compared across parameter estimation methods
% Load the results data from combined_results.csv
results = readtable('combined_results.csv');

% Initialize arrays to store error metrics and runtimes for each method
methods = {'LSQ', 'FMIN', 'FITNLM', 'FWD', 'BWD', 'CEN'};
mean_absolute_errors = zeros(1, 6);
mean_squared_errors = zeros(1, 6);
root_mean_squared_errors = zeros(1, 6);
run_times = zeros(1, 6);

% Extract error metrics and runtimes for each method
mean_absolute_errors(1) = mean(results.MAE_LSQ);
mean_squared_errors(1) = mean(results.MSE_LSQ);
root_mean_squared_errors(1) = mean(results.RMSE_LSQ);
run_times(1) = mean(results.RunTime_LSQ);

mean_absolute_errors(2) = mean(results.MAE_FMIN);
mean_squared_errors(2) = mean(results.MSE_FMIN);
root_mean_squared_errors(2) = mean(results.RMSE_FMIN);
run_times(2) = mean(results.RunTime_FMIN);

mean_absolute_errors(3) = mean(results.MAE_FITNLM);
mean_squared_errors(3) = mean(results.MSE_FITNLM);
root_mean_squared_errors(3) = mean(results.RMSE_FITNLM);
run_times(3) = mean(results.RunTime_FITNLM);

mean_absolute_errors(4) = mean(results.MAE_FWD);
mean_squared_errors(4) = mean(results.MSE_FWD);
root_mean_squared_errors(4) = mean(results.RMSE_FWD);
run_times(4) = mean(results.RunTime_FWD);

mean_absolute_errors(5) = mean(results.MAE_BWD);
mean_squared_errors(5) = mean(results.MSE_BWD);
root_mean_squared_errors(5) = mean(results.RMSE_BWD);
run_times(5) = mean(results.RunTime_BWD);

mean_absolute_errors(6) = mean(results.MAE_CEN);
mean_squared_errors(6) = mean(results.MSE_CEN);
root_mean_squared_errors(6) = mean(results.RMSE_CEN);
run_times(6) = mean(results.RunTime_CEN);

% Create a table for RMSE and Run Times
comparison_table = table(methods', root_mean_squared_errors', run_times', ...
    'VariableNames', {'Method', 'RMSE', 'RunTime'});

% Display the comparison table
disp('Comparison Table:');
disp(comparison_table);

% Save the comparison table to a new CSV file
writetable(comparison_table, 'comparison_results_table.csv');

disp('Comparison table has been saved to comparison_results_table.csv');

% Create a figure for the table
f = figure('Position', [100, 100, 600, 200]);
uitable('Parent', f, 'Data', [comparison_table.Method, num2cell(comparison_table.RMSE), num2cell(comparison_table.RunTime)], ...
        'ColumnName', comparison_table.Properties.VariableNames, ...
        'RowName', {}, ...
        'Units', 'Normalized', ...
        'Position', [0, 0, 1, 1]);

% Save the table as an image
saveas(f, 'comparison_results_table.png');

disp('Comparison table image has been saved as comparison_results_table.png');
