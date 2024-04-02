function [performanceLDALedoitHbO,performanceLDALedoitHbO_valtrain] = vector_feature(trials_HbOM_Tr,trials_HbRM_Tr,indexMoviesTrainMultiple,trials_HbOM_Tst,trials_HbRM_Tst,indexMoviesTestMultiple)


    fs = 8.988;

    timeLen = [2].*fs; %how long the response is

    cue_start = ceil(2*fs);

    %%

    windowStartT = 2;
    windowEndT = 10;

    %%

    startT = cue_start + ceil(windowStartT*fs);
    endT = cue_start + ceil(windowEndT*fs);


    % offset trials
    % trials_HbOM_Tr = offsetTrials(trials_HbOM_Tr,cue_start);
    % trials_HbRM_Tr = offsetTrials(trials_HbRM_Tr,cue_start);
    % 
    % trials_HbOM_Tst = offsetTrials(trials_HbOM_Tst,cue_start);
    % trials_HbRM_Tst = offsetTrials(trials_HbRM_Tst,cue_start);

    %% vector phase extraction

    deltaHbO_Tr = trials_HbOM_Tr(:,startT:endT,:);
    deltaHbR_Tr = trials_HbRM_Tr(:,startT:endT,:);

    deltaHbO_Tst = trials_HbOM_Tst(:,startT:endT,:);
    deltaHbR_Tst = trials_HbRM_Tst(:,startT:endT,:);

    deltaCBV_Tr = zeros(size(deltaHbO_Tr));
    deltaCOE_Tr = zeros(size(deltaHbO_Tr));
    magnitudeL_Tr = zeros(size(deltaHbO_Tr));
    angleK_Tr = zeros(size(deltaHbO_Tr));

    deltaCBV_Tst = zeros(size(deltaHbO_Tst));
    deltaCOE_Tst = zeros(size(deltaHbO_Tst));
    magnitudeL_Tst = zeros(size(deltaHbO_Tst));
    angleK_Tst = zeros(size(deltaHbO_Tst));


    for trial = 1:size(deltaHbO_Tr, 3)
        for channel = 1:size(deltaHbO_Tr, 1)
            for timePoint = 1:size(deltaHbO_Tr, 2)
                % Extract the values for the current channel and time point
                currentDeltaHbO_Tr = deltaHbO_Tr(channel, timePoint, trial);
                currentDeltaHbR_Tr = deltaHbR_Tr(channel, timePoint, trial);
    
                % Corrected calculation of each feature
                deltaCBV_Tr(channel, timePoint, trial) = (currentDeltaHbO_Tr + currentDeltaHbR_Tr) ./ sqrt(2);
                deltaCOE_Tr(channel, timePoint, trial) = (currentDeltaHbR_Tr - currentDeltaHbO_Tr) ./ sqrt(2);
                magnitudeL_Tr(channel, timePoint, trial) = sqrt((currentDeltaHbO_Tr.^2 + currentDeltaHbR_Tr.^2) ./ 2);
                angleK_Tr(channel, timePoint, trial) = atan2(deltaCOE_Tr(channel, timePoint, trial), deltaCBV_Tr(channel, timePoint, trial));
            end
        end
    end

    for trial = 1:size(deltaHbO_Tst, 3)
        for channel = 1:size(deltaHbO_Tst, 1)
            for timePoint = 1:size(deltaHbO_Tst, 2)
                % ExTstact the values for the current channel and time point
                currentDeltaHbO_Tst = deltaHbO_Tst(channel, timePoint, trial);
                currentDeltaHbR_Tst = deltaHbR_Tst(channel, timePoint, trial);
    
                % Corrected calculation of each feature
                deltaCBV_Tst(channel, timePoint, trial) = (currentDeltaHbO_Tst + currentDeltaHbR_Tst) ./ sqrt(2);
                deltaCOE_Tst(channel, timePoint, trial) = (currentDeltaHbR_Tst - currentDeltaHbO_Tst) ./ sqrt(2);
                magnitudeL_Tst(channel, timePoint, trial) = sqrt((currentDeltaHbO_Tst.^2 + currentDeltaHbR_Tst.^2) ./ 2);
                angleK_Tst(channel, timePoint, trial) = atan2(deltaCOE_Tst(channel, timePoint, trial), deltaCBV_Tst(channel, timePoint, trial));
            end
        end
    end

    % active_channels_logical = logical(active_channels);
    % 
    % deltaCBV = deltaCBV(active_channels_logical,:,:);
    % deltaCOE = deltaCOE(active_channels_logical,:,:);
    % magnitudeL = magnitudeL(active_channels_logical,:,:);
    % angleK = angleK(active_channels_logical,:,:);

    %% concatenate features

        % Number of features per channel
    numFeatures = 4; 
    
    % Total number of channels
    numChannels = size(deltaCBV_Tr,1);
    
    % Initialize the feature vector
    % Total features = numChannels * numFeatures
    % Total data points = number of trials * number of time points
    featureVector_Tr = zeros(size(deltaCOE_Tr,3) * size(deltaCOE_Tr,2), numChannels * numFeatures);
    featureVector_Tst = zeros(size(deltaCOE_Tst,3) * size(deltaCOE_Tst,2), numChannels * numFeatures);

    
    % Loop through each trial and time point
    for trial = 1:size(deltaCOE_Tr,3)
        for timePoint = 1:size(deltaCOE_Tr,2)
            % Index for rows in the feature vector
            rowIndex = (trial - 1) * size(deltaCOE_Tr,2) + timePoint;
    
            % Extract features for all channels and concatenate
            for channel = 1:numChannels
                colIndex = (channel - 1) * numFeatures;
                featureVector_Tr(rowIndex, colIndex + 1) = deltaCBV_Tr(channel, timePoint, trial);
                featureVector_Tr(rowIndex, colIndex + 2) = deltaCOE_Tr(channel, timePoint, trial);
                featureVector_Tr(rowIndex, colIndex + 3) = magnitudeL_Tr(channel, timePoint, trial);
                featureVector_Tr(rowIndex, colIndex + 4) = angleK_Tr(channel, timePoint, trial);
            end
        end
    end

    for trial = 1:size(deltaCOE_Tst,3)
        for timePoint = 1:size(deltaCOE_Tst,2)
            % Index for rows in the feature vector
            rowIndex = (trial - 1) * size(deltaCOE_Tst,2) + timePoint;
    
            % Extract features for all channels and concatenate
            for channel = 1:numChannels
                colIndex = (channel - 1) * numFeatures;
                featureVector_Tst(rowIndex, colIndex + 1) = deltaCBV_Tst(channel, timePoint, trial);
                featureVector_Tst(rowIndex, colIndex + 2) = deltaCOE_Tst(channel, timePoint, trial);
                featureVector_Tst(rowIndex, colIndex + 3) = magnitudeL_Tst(channel, timePoint, trial);
                featureVector_Tst(rowIndex, colIndex + 4) = angleK_Tst(channel, timePoint, trial);
            end
        end
    end

        % Finding columns that have all zeros
    colsToRemove_Tr = all(featureVector_Tr == 0);
    colsToRemove_Tst = all(featureVector_Tst == 0);

    
    % Removing these columns
    featureVector_Tr(:, colsToRemove_Tr) = [];
    featureVector_Tst(:, colsToRemove_Tst) = [];


    expandedLabels_Tr = zeros(size(deltaCOE_Tr,3) * size(deltaCOE_Tr,2), 1); % Initialize expanded label vector
    expandedLabels_Tst = zeros(size(deltaCOE_Tst,3) * size(deltaCOE_Tst,2), 1); % Initialize expanded label vector

    response_data_Tr = indexMoviesTrainMultiple';
    response_data_Tst = indexMoviesTestMultiple';


    %%
    
    for trial = 1:size(deltaCOE_Tr,3)
        for timePoint = 1:size(deltaCOE_Tr,2)
            % Calculate the row index for the expanded labels
            rowIndex = (trial - 1) * size(deltaCOE_Tr,2) + timePoint;
    
            % Assign the label of the current trial to the corresponding time points
            expandedLabels_Tr(rowIndex) = response_data_Tr(trial);
        end
    end

    for trial = 1:size(deltaCOE_Tst,3)
        for timePoint = 1:size(deltaCOE_Tst,2)
            % Calculate the row index for the expanded labels
            rowIndex = (trial - 1) * size(deltaCOE_Tst,2) + timePoint;
    
            % Assign the label of the current trial to the corresponding time points
            expandedLabels_Tst(rowIndex) = response_data_Tst(trial);
        end
    end

    performanceLDALedoitHbO = zeros(1,length(timeLen));
    performanceLDALedoitHbO_valtrain = zeros(1,length(timeLen));


    for i2 = 1:length(timeLen)

      [performanceLDALedoitHbO(1,i2),performanceLDALedoitHbO_valtrain(1,i2) ] = trainClassifier_qda_sudan_v3(featureVector_Tr,expandedLabels_Tr,featureVector_Tst,expandedLabels_Tst);
       
    end  
end