%% takes filepaths and gives labels

% Read the CSV file

function labels = label_creater(filepath)


data = readtable(filepath); % Replace with the actual file path

% Extract the 'TargetSpatialSuffix' column
targetSpatialSuffix = data.TargetSpatialSuffix;

% Initialize an array to store the converted values
outputArray = zeros(length(targetSpatialSuffix), 1);

% Loop through the column and assign values
for i = 1:length(targetSpatialSuffix)
    if strcmp(targetSpatialSuffix{i}, '_30')
        outputArray(i) = 1;
    elseif strcmp(targetSpatialSuffix{i}, '_-30')
        outputArray(i) = 0;
    else
        % Handle other cases if necessary, or throw an error
        error('Unexpected value in TargetSpatialSuffix at row %d', i);
    end
end

labels = outputArray';

end

