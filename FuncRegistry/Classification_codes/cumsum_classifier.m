% STATUS: Inactive(?)
% 
% SYNTAX:
% [performanceLDALedoitHbO] = ...
%   calcCumSum_DiffTLen_AllChn_CompClassifiers_Final_NoSave(sbjNum,...
%       trialsTr,trialsTst,movieListTrain,movieListTest,numClasses,...
%       mlActAuto,timeLen)
% 
% DESCRIPTION:
% Perform cross validation where SS beta coefficients from training fold 
%   is passed to test fold. Perform all-channel classification using rLDA
%   classifier. Test different decision window lengths.
% 
% RESTRICTION:
% None.
% 
% INPUTS:
% sbjNum - string: subject ID
% trialsTr - training dataset. 3D double array: channel x time x trial
% trialsTst - test dataset. 3D double array: channel x time x trial
% movieListTrain - trials info for training dataset. numTrials x 5 double array:
%       col 1: index of target movies in uniqueMovies
%       col 2: index of spatial location
%       col 3: boolean: masker is fixed or random
%       col 4: index of masker movies in fixedMaskerList(?)
%       col 5: boolean: condition is target-alone or target+maskers
% movieListTest - trials info for test dataset. same structure as movieListTrain
% numClasses - int: number of classes for classification.
% mlActAuto - 1x1 cell array containing 1D int array of channel list of 2
%   different wavelengths
% timeLen - 1D double array of window lengths to test, in sec.
%
% RETURNED VARIABLES:
% performanceLDALedoitHbO - decoding performance of rLDA classifier for 
%   each decision window. 1D double array: 1 x time windows
% performanceLDACERNNHbO - decoding performance of CERNN LDA classifier for 
%   each decision window. 1D double array: 1 x time windows
% performanceCosineKNN - decoding performance of KNN classifier for 
%   each decision window. 1D double array: 1 x time windows
% performanceSVMHbO - decoding performance of SVM classifier for 
%   each decision window. 1D double array: 1 x time windows
% performanceBaggingHbO - decoding performance of random forest ensemble 
%   classifier for each decision window. 1D double array: 1 x time windows
% performanceBoostedHbO - decoding performance of boosting ensemble
%   classifier for each decision window. 1D double array: 1 x time windows
% 
% FILES SAVED:
% None.
% 
% PLOTTING:
% None.

function [performanceLDALedoitHbO] = cumsum_classifier(trialsTr,trialsTst,indexMoviesTrain,indexMoviesTest,timeLen, mlList,chOI,fs)

startT = ceil(2*fs);

numChn = size(trialsTr,1);

performanceLDALedoitHbO = zeros(1,length(timeLen));

chOI = chOI(chOI <= length(mlList));
updatedMlList = zeros(size(mlList));

% Set ones at specified indices if original mlList had ones at these indices
updatedMlList(chOI) = mlList(chOI) & 1;

mlList = logical(updatedMlList);

for i2 = 1:length(timeLen)
    
    temp = cumsum(trialsTr(mlList,startT:startT+ceil(timeLen(i2)),:),2);
    cumsumTr = squeeze(temp(:,end,:));

    temp = cumsum(trialsTst(mlList,startT:startT+ceil(timeLen(i2)),:),2);
    cumsumTst = squeeze(temp(:,end,:));

    [performanceLDALedoitHbO(1,i2)] = train_RLDA_Ledoit_TrTst_sudan(...
    cumsumTr,cumsumTst,2,indexMoviesTrain,indexMoviesTest);
       
end

end