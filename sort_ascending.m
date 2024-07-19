% sort data filet to ensure ascending values for each patient in the
% Treatment_day field

% Load the dataset
data = readtable('combined_filtered_tumor_data_example01.csv'); 

% Sort the data by patientID and Treatment_Day
sorted_data = sortrows(data, {'patientID', 'Treatment_Day'});

% Save the sorted data to a new file
writetable(sorted_data, 'sorted_tumor_size_example_study1_pt2_02.csv');

% Display a message
disp('Data has been sorted and saved to sorted_tumor_size_example_study1_pt2_02.csv');
