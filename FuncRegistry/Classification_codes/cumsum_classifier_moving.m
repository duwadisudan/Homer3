function [performanceLDALedoitHbO] = cumsum_classifier_moving(trialsTr, trialsTst, indexMoviesTrain, indexMoviesTest,timeLen, mlList, chOI, fs)

    % Define window parameters
    windowDuration = 4 * fs; % 4 seconds
    stepSize = 2 * fs; % shift of 2 seconds for each window
    numWindows = 5; % Total number of windows

    performanceLDALedoitHbO = zeros(1, numWindows);

    chOI = chOI(chOI <= length(mlList));
    updatedMlList = zeros(size(mlList));
    updatedMlList(chOI) = mlList(chOI) & 1;
    mlList = logical(updatedMlList);

    for winIdx = 1:numWindows
        windowStart = 1 + (winIdx-1) * stepSize;
        windowEnd = windowStart + windowDuration - 1;

        % Cumulative sum for the training data in the current window
        temp = cumsum(trialsTr(mlList, windowStart:windowEnd, :), 2);
        cumsumTr = squeeze(temp(:, end, :));

        % Cumulative sum for the testing data in the current window
        temp = cumsum(trialsTst(mlList, windowStart:windowEnd, :), 2);
        cumsumTst = squeeze(temp(:, end, :));

        % Classification for the current window
        [performanceLDALedoitHbO(1,winIdx)] = train_RLDA_Ledoit_TrTst_sudan(...
        cumsumTr, cumsumTst, 2, indexMoviesTrain, indexMoviesTest);
    end
end
