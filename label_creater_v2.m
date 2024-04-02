function labels = label_creater_v2(filepath)

data = readtable(filepath); % Read the CSV file

% Check for the presence of 'TargetSpatialSuffix' or 'TargetScreen' column
if any(strcmp(data.Properties.VariableNames, 'TargetSpatialSuffix'))
    targetColumn = data.TargetSpatialSuffix;
    isCharType = true;
elseif any(strcmp(data.Properties.VariableNames, 'TargetScreen'))
    targetColumn = data.TargetScreen;
    isCharType = false;
else
    error('Expected column not found in the file');
end

% Initialize an array to store the converted values
outputArray = zeros(length(targetColumn), 1);

% Loop through the column and assign values
for i = 1:length(targetColumn)
    if isCharType % For 'TargetSpatialSuffix' type
        if strcmp(targetColumn{i}, '_30')
            outputArray(i) = 1;
        elseif strcmp(targetColumn{i}, '_-30')
            outputArray(i) = 0;
        else
            error('Unexpected value in TargetSpatialSuffix at row %d', i);
        end
    else % For 'TargetScreen' type
        if targetColumn(i) == 1
            outputArray(i) = 1; % Right
        elseif targetColumn(i) == 0
            outputArray(i) = 0; % Left
        else
            error('Unexpected value in TargetScreen at row %d', i);
        end
    end
end

labels = outputArray';

end
