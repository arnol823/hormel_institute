% Generates two graphs: one for RMSE and one for run time based on different methods of parameter approximation
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

% Display the error metrics
disp('Mean Absolute Errors for each method:'); 
disp(array2table(mean_absolute_errors, 'VariableNames', methods));

disp('Mean Squared Errors for each method:');
disp(array2table(mean_squared_errors, 'VariableNames', methods));

disp('Root Mean Squared Errors for each method:');
disp(array2table(root_mean_squared_errors, 'VariableNames', methods));

disp('Run Times for each method:');
disp(array2table(run_times, 'VariableNames', methods));

% Create a bar chart for RMSE
figure;
bar(root_mean_squared_errors, 'FaceColor', 'b');
set(gca, 'XTickLabel', methods, 'FontSize', 14);
% title('Comparison of RMSE for Each Method', 'FontSize', 24);
xlabel('Methods', 'FontSize', 24);
ylabel('RMSE', 'FontSize', 24);
ylim([0 max(root_mean_squared_errors)*1.1]);
xtips1 = 1:length(root_mean_squared_errors);
ytips1 = root_mean_squared_errors;
labels1 = arrayfun(@(x) sprintf('%.2e', x), root_mean_squared_errors, 'UniformOutput', false);
text(xtips1, ytips1, labels1, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 14);

% Create a bar chart for Run Times
figure;
bar(run_times, 'FaceColor', 'r', 'FaceAlpha', 0.5);
set(gca, 'XTickLabel', methods, 'FontSize', 14);
% title('Comparison of Run Times for Each Method', 'FontSize', 24);
xlabel('Methods', 'FontSize', 24);
ylabel('Run Time (s)', 'FontSize', 24);
ylim([0 max(run_times)*1.1]);
xtips2 = 1:length(run_times);
ytips2 = run_times;
labels2 = arrayfun(@(x) sprintf('%.2e', x), run_times, 'UniformOutput', false);
text(xtips2, ytips2, labels2, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 14);

% Save the comparison results to a new CSV file
comparison_results = table(mean_absolute_errors', mean_squared_errors', root_mean_squared_errors', run_times', ...
    'VariableNames', {'MAE', 'MSE', 'RMSE', 'RunTime'}, 'RowNames', methods);
writetable(comparison_results, 'comparison_results.csv', 'WriteRowNames', true);

disp('Comparison results have been saved to comparison_results.csv');
