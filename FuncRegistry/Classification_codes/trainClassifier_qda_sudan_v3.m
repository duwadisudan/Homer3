function [validationAccuracy,TrainvalidationAccuracy] = trainClassifier_qda_sudan_v3(trainData, trainResponse, testData, testResponse)
    % Transpose trainData and testData
    %trainData = trainData';  % Now 48 x 567
    %testData = testData';    % Now 12 x 567

    % Apply PCA to transposed training data
    [coefficients, score, ~, ~, explained, mu] = pca(trainData);
    % Keep components that explain a specified variance
    explainedVarianceToKeep = 95;
    idx = find(cumsum(explained) >= explainedVarianceToKeep, 1, 'first');
    pcaTrainData = score(:, 1:idx);

    %scoreTest95 = (XTest-mu)*coeff(:,1:idx);
    classificationDiscriminant = fitcdiscr(pcaTrainData, trainResponse, 'DiscrimType', 'linear','Gamma', 0.5, 'Delta', 0.5);
    %classificationDiscriminant = fitcdiscr(pcaTrainData, trainResponse, 'DiscrimType', 'quadratic');

    % Transform transposed testing data using PCA parameters from training data
    %pcaTestData = bsxfun(@minus, testData, mu) * coefficients(:, 1:idx);

    pcaTestData = (testData - mu) * coefficients(:, 1:idx);

    trainPredictions = predict(classificationDiscriminant, pcaTrainData);

    % Make predictions on PCA-transformed testing data
    testPredictions = predict(classificationDiscriminant, pcaTestData);

    % Calculate validation accuracy
    TraincorrectPredictions = (trainPredictions == trainResponse);
    TrainvalidationAccuracy = sum(TraincorrectPredictions) / numel(trainResponse);

    % Calculate validation accuracy
    correctPredictions = (testPredictions == testResponse);
    validationAccuracy = sum(correctPredictions) / numel(testResponse);

    % Create the result struct
    trainedClassifier.ClassificationDiscriminant = classificationDiscriminant;
    trainedClassifier.PCACoefficients = coefficients;
    trainedClassifier.PCACenters = mu;
    trainedClassifier.About = 'This struct is a trained model using custom train-test split.';
end
