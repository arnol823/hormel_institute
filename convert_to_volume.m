% Load the dataset
data = readtable('combined_filtered_tumor_data_example01.csv');

% Extract the diameter column
diameter = data.TargetLesionLongDiam_mm;

% Calculate the volume of the sphere
volume = (4/3) * pi * (diameter / 2).^3;

% Add the volume column to the table
data.TargetLesionLongDiam_mm = volume;

% Save the updated table to a new CSV file
writetable(data, 'combined_filtered_tumor_data_with_volume.csv');

disp('The updated data with volumes has been saved to combined_filtered_tumor_data_with_volume.csv');
