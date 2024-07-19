% Concatenate filtered tumor size and patient data into one file
% List of CSV files to concatenate
fileList = {'filtered_tumor_data_s1.csv', 'filtered_tumor_data_s2.csv', ...
            'filtered_tumor_data_s3.csv', 'filtered_tumor_data_s4.csv', ...
            'filtered_tumor_data_s5.csv'};

% Initialize an empty table to store concatenated data
combinedData = table();

% Loop through each file and concatenate the data
for i = 1:length(fileList)
    % Read the current CSV file
    currentData = readtable(fileList{i});
    
    % Concatenate the current data with the combined data
    combinedData = [combinedData; currentData];
end

% Write the combined data to a new CSV file
outputFileName = 'combined_filtered_tumor_data.csv';
writetable(combinedData, outputFileName);

% Display a message indicating that the concatenation is complete
disp(['Data concatenation complete. Combined data saved to ' outputFileName '.']);


%%
data = readmatrix('filtered_tumor_data_s1.csv');
x1 = size(data)

data = readmatrix('filtered_tumor_data_s2.csv');
x2 = size(data)

data = readmatrix('filtered_tumor_data_s3.csv');
x3 = size(data)

data = readmatrix('filtered_tumor_data_s4.csv');
x4 = size(data)

data = readmatrix('filtered_tumor_data_s5.csv');
x5 = size(data)

z = x1+x2 + x3 + x4 + x5

data = readmatrix('combined_filtered_tumor_data.csv');
z_true = size(data)

z == size(data)

% file has correct number of rows and three columns, so good!
