function [validationScore] = train_RLDA_Ledoit_TrTst_sudan(...
    trialsTr,trialsTst,numClasses,indexMoviesTrain,indexMoviesTest)

    trainingPredictors = trialsTr;
    trainingPredictors(isinf(trainingPredictors)) = NaN;
    trainingResponse = indexMoviesTrain;
    trainingResponse(trainingResponse == 0) = 2;

    testPredictors = trialsTst;
    testPredictors(isinf(testPredictors)) = NaN;
    testResponse = indexMoviesTest;  
    testResponse(testResponse == 0) = 2;
    

    numChn = size(testPredictors,1);

    
    if size(trainingPredictors,2) == 0 || size(testPredictors,2) == 0
        validationScore = 0;
        label = zeros(size(testResponse));
    else

        numTr = size(testPredictors,2);        
        % Use a priori uniform distribution instead of empirical prior.
        prior = 1/numClasses;
        
        discriminantScore = zeros(numClasses,numTr);
        
        innerSigma = zeros(numClasses,numChn,numChn);
        
        for i3 = 1:numClasses
            thisClassIdx = trainingResponse == i3;
            thisClassX = trainingPredictors(:,thisClassIdx);
            thisClassMu = mean(thisClassX,2);
            innerSigma(i3,:,:) = (thisClassX-thisClassMu)*(thisClassX-thisClassMu)';
        end
        
        pooledS = squeeze(sum(innerSigma,1))/(size(trainingPredictors,2)-numClasses);
        
        sigma = calcShrinkCovMatrix(trainingPredictors,pooledS);
        
        for i3 = 1:numClasses
            % Is sphering the data necessary? No. Just zero-centering is
            % enough.
            thisClassIdx = trainingResponse == i3;
            thisClassX = trainingPredictors(:,thisClassIdx);
            
            mu = mean(thisClassX,2);
            
            discriminantScore(i3,:) = calcDiscriminantFunction(testPredictors,sigma,mu,prior);
            
        end
        
        [~,label] = max(discriminantScore,[],1);

        temp = squeeze(label==testResponse);
        validationScore = sum(temp)/length(temp);
    end
    
end
    
    


