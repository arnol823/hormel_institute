% filter down tumor dataset to include only patients with 6 points or more
% Load the data from the CSV file
data = readtable('tumor_size02_study05.csv');

% Extract relevant columns
patientID = data.Patient_Anonmyized; % Assuming the first column is named 'PatientID'
time = data.Treatment_Day; % Assuming the second column is named 'Treatment_Day'
tumor_size = data.TargetLesionLongDiam_mm; % Assuming the third column is named 'Tumor_Size'

% Find unique patient IDs
uniquePatients = unique(patientID);

% Filter out patients with less than six data points
minDataPoints = 5;
patientsToInclude = uniquePatients(arrayfun(@(id) sum(patientID == id) > minDataPoints, uniquePatients));

% Initialize arrays for filtered data
filteredPatientID = [];
filteredTime = [];
filteredTumorSize = [];

% Filter data
for i = 1:length(patientsToInclude)
    patientData = patientID == patientsToInclude(i);
    filteredPatientID = [filteredPatientID; patientID(patientData)];
    filteredTime = [filteredTime; time(patientData)];
    filteredTumorSize = [filteredTumorSize; tumor_size(patientData)];
end

% Create a new table with the filtered data
filteredData = table(filteredPatientID, filteredTime, filteredTumorSize, ...
                     'VariableNames', {'PatientID', 'Treatment_Day', 'Tumor_Size'});

% Save the filtered data to a new CSV file
writetable(filteredData, 'filtered_tumor_data_s5.csv');

% Display a message indicating that the filtering is complete
disp('Data filtering complete. Filtered data saved to filtered_tumor_data.csv.');
