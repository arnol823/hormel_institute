% number of unique values and the number of their data entries for tumor
% size dataset
% 
% data = readtable('combined_filtered_tumor_data.csv');%
% columnData = data.PatientID; 
% [uniqueValues, ~, idx] = unique(columnData);
% 
% % Count the number of entries for each unique value
% counts = accumarray(idx, 1);
% 
% % Display the unique values and their counts
% disp(table(uniqueValues, counts, 'VariableNames', {columnData, 'Count'}));

% Function to count unique values in a specific column of a CSV file
function uniqueValueCounts = countUniqueValues(filename, columnName)
    % Read the table from the CSV file
    data = readtable(filename);
    
    % Extract the column data
    columnData = data.(columnName);
    
    % Find the unique values and their counts
    [uniqueValues, ~, idx] = unique(columnData);
    counts = accumarray(idx, 1);
    
    % Create a table of unique values and their counts
    uniqueValueCounts = table(uniqueValues, counts, 'VariableNames', {columnName, 'Count'});
    
    % Display the table
    disp(uniqueValueCounts);
end

% Example usage
uniqueValueCounts = countUniqueValues('combined_filtered_tumor_data.csv', 'PatientID');

% Verify that all counts are above 5
    for i = 1:height(uniqueValueCounts)
        if uniqueValueCounts.Count(i) <= 5
            disp(['Count for ' char(uniqueValueCounts.(columnName)(i)) ' is 5 or less: ' num2str(uniqueValueCounts.Count(i))]);
        end

    end
