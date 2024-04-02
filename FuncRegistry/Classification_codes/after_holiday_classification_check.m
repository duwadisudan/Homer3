%%
close all;
clear all;
clc;
%%
cd 'H:\My Drive\Labs\BOAS Lab\Homer3'
setpaths
%%
% Define the base directory and runs

baseDirectory = 'H:\My Drive\Labs\Sen Lab\Pilot_data_processing_Laura';
runs = {
    struct('name', 'actual_exp_fef', 'path', fullfile(baseDirectory, 'actual_exp_fef', 'ninjaNIRS2022_2023-11-07-16-56-53'), 'indexMoviesTest', [1; 0; 0; 0; 0; 1; 1; 1; 1; 1; 0; 0; 0; 0; 1; 0; 1; 0; 0; 1; 0; 1; 1]),
    struct('name', 'covert', 'path', fullfile(baseDirectory, 'covert', 'ninjaNIRS2022_2023-11-07-16-10-04'), 'indexMoviesTest', [0; 1; 0; 1; 1; 0; 0; 0; 1; 1; 0; 1; 0; 0; 1; 0; 0; 0; 0; 0; 0; 0; 1; 0; 0; 0; 1; 0; 1; 1; 1; 1; 1; 0; 0; 0; 0; 0; 1; 0; 0; 0; 0; 1; 1; 0; 1; 1; 1; 1; 1; 0; 1; 0; 0;]),
    struct('name', 'eye_control_fef', 'path', fullfile(baseDirectory, 'eye_control_fef', 'ninjaNIRS2022_2023-11-07-15-12-01'), 'indexMoviesTest', [1; 0; 0; 0; 1; 1; 1; 0; 0; 1; 0; 0; 1; 1; 1; 0; 1; 1; 0; 0]),
    struct('name', 'eye_control_no_stim_fef', 'path', fullfile(baseDirectory, 'eye_control_no_stim_fef', 'ninjaNIRS2022_2023-11-14-14-04-22'), 'indexMoviesTest', [1; 1; 1; 0; 0; 0; 0; 0; 1; 0; 0; 0; 1; 1; 1; 0; 1; 0; 1; 0]),
    struct('name', 'visual_only', 'path', fullfile(baseDirectory, 'visual_only', 'ninjaNIRS2022_2023-11-07-15-33-23'), 'indexMoviesTest', [0; 0; 0; 1; 1; 0; 0; 1; 0; 0; 0; 0; 1; 0; 1; 1; 1; 1; 1; 1; 0; 1; 1; 0; 0; 1; 1; 1; 0; 0; 1; 1; 0; 1; 0; 1; 1; 0; 0; 1; 1; 0; 1; 0; 1; 1; 0; 1; 1; 1; 0; 1; 0; 1; 1; 0; 0; 0; 1; 0])
    struct('name', 'covert_laura', 'path', fullfile(baseDirectory, 'covert_laura', 'ninjaNIRS2022_2024-01-10-14-36-11'), 'indexMoviesTest', [0;0;0;1;1;0;0;1;1;1;1;0;0;1;1;1;0;1;0;0;0;0;0;0;0;1;1;0;0;0;1;0;1;0;1;0;1;0;1;1;0;1;1;0;0;0;0;0;0;1;1;0;1;1;1;1;0;0;1;0])
    struct('name', 'covert_jonathan', 'path', fullfile(baseDirectory, 'covert_jonathan', 'ninjaNIRS2022_2024-01-11-16-17-24'), 'indexMoviesTest', [0;1;0;1;1;1;1;0;1;1;0;1;0;0;0;1;1;0;0;1;0;1;1;1;1;1;1;1;1;0;1;0;1;0;0;1;1;1;0;1;0;1;0;0;1;1;0;0;1;0;0;0;1;1;0;1;0;1;0;0])
    struct('name', 'covert_jonathan', 'path', fullfile(baseDirectory, 'covert_jonathan', 'ninjaNIRS2022_2024-01-11-16-35-14'), 'indexMoviesTest', [1;0;1;0;1;0;0;1;1;1;0;1;0;1;0;0;1;1;0;0;1;0;0;0;1;1;0;1;0;1;0;0])
    struct('name', 'covert_jonathan', 'path', fullfile(baseDirectory, 'covert_jonathan', 'ninjaNIRS2022_2024-01-11-17-01-03'), 'indexMoviesTest', [1;1;1;1;0;1;0;1;0;1;0;1;1;0;0;0;0;0;1;0;0;0;0;0;0;0])
    struct('name', 'covert_jonathan', 'path', fullfile(baseDirectory, 'covert_jonathan', 'covert_jonathan_concat'), 'indexMoviesTest', [0;1;0;1;1;1;1;0;1;1;0;1;0;0;0;1;1;0;0;1;0;1;1;1;1;1;1;1;1;0;1;0;1;0;0;1;1;1;0;1;0;1;0;0;1;1;0;0;1;0;0;0;1;1;0;1;0;1;0;0])

    };

%%

currentRun = runs{5};


%load("H:\My Drive\Labs\Sen Lab\Pilot_data_processing_Laura\Actual_exp_fef\ninjaNIRS2022_2023-11-07-16-56-53\top_ten_channels\top_ten_channels_actual_exp_fef.mat")
%load('H:\My Drive\Labs\Sen Lab\Pilot_data_processing_Laura\covert\ninjaNIRS2022_2023-11-07-16-10-04\top_ten_channels\top_ten_channels_covert.mat')



num_ch_selected = 5;


[snirf, time, events] = loadDataAndPreprocess(currentRun.path);

% Process events and triggers
[allS, indexMoviesTest] = processEvents(events, time, currentRun.indexMoviesTest);

%% prune channels

fs = 8.988;

mlActAuto = hmrR_PruneChannels(snirf.data, snirf.probe, [], [], [1e-3, 1e7], 10, [0.0, 45.0]);
%% propagation of pruned channels

mlActive = mlActAuto{1};
activeBinary = mlActive(:,3);


firstWavelengthChannel_idx = activeBinary(1:567);
secondWavelengthChannel_idx = activeBinary(568:1134);

secondWavelengthChannel_idx(~firstWavelengthChannel_idx) = 0;

firstWavelengthChannel_idx(~secondWavelengthChannel_idx) = 0;

activeBinary = [firstWavelengthChannel_idx; secondWavelengthChannel_idx];

total_num_channels_check = sum(activeBinary(1:567) == activeBinary(568:end));

mlActive(:,3) = activeBinary;

mlActAuto{1} = mlActive;

%%

% Separate into left and right stimuli based on indexMoviesTest
leftMultiIndex = find(indexMoviesTest == 0);
rightMultiIndex = find(indexMoviesTest == 1);

% Populate stim class for left and right stimuli
snirf.stim(1).data = [allS(leftMultiIndex), ones(length(leftMultiIndex), 1) * 10, ones(length(leftMultiIndex), 1)]; %what is stim duration?
snirf.stim(1).states = [allS(leftMultiIndex), ones(length(leftMultiIndex), 1)];
snirf.stim(1).name = 'left';

snirf.stim(2).data = [allS(rightMultiIndex), ones(length(rightMultiIndex), 1) * 10, ones(length(rightMultiIndex), 1)];
snirf.stim(2).states = [allS(rightMultiIndex), ones(length(rightMultiIndex), 1)];
snirf.stim(2).name = 'right';

%%
timeLen = [0.1 0.2 0.5 1 1.5 2 3 4 5 6 7 8].*fs; %how long the response is
    
nrept = 10;
kFold = 5;

num_of_trials = size(allS,1);

%% conversion into optical dendity
dod = hmrR_Intensity2OD(snirf.data);

%% Initialize perrformance

performanceLDALedoitHbO_All = zeros(nrept*kFold,length(timeLen));
performanceLDALedoitHbR_All = zeros(nrept*kFold,length(timeLen));
performanceLDALedoitHbT_All = zeros(nrept*kFold,length(timeLen));

% Left FEF
% performanceLDALedoitHbO_LFEF = zeros(nrept*kFold,length(timeLen));
% performanceLDALedoitHbR_LFEF = zeros(nrept*kFold,length(timeLen));
% performanceLDALedoitHbT_LFEF = zeros(nrept*kFold,length(timeLen));
% 
% % Right FEF
% performanceLDALedoitHbO_RFEF = zeros(nrept*kFold,length(timeLen));
% performanceLDALedoitHbR_RFEF = zeros(nrept*kFold,length(timeLen));
% performanceLDALedoitHbT_RFEF = zeros(nrept*kFold,length(timeLen));


%% Change Homer and addpaths

cd 'H:\My Drive\Labs\Sen Lab\NN22_processing\ProcessingCodefNIRS-main\Homer3-master'

setpaths %homer function to setpaths for its dependencies

addpath 'H:\My Drive\Labs\Sen Lab\NN22_processing\ProcessingCodeClassifier-main'
addpath 'H:\My Drive\Labs\Sen Lab\NN22_processing\ProcessingCodefNIRS-main'
addpath 'H:\My Drive\Labs\Sen Lab\NN22_processing\GroupProcessingfNIRS-main'
addpath 'H:\My Drive\Labs\Sen Lab\NN22_processing\ProcessingCodefNIRS-main\Classification'
%%
%% CV loop

for iRep = 1:nrept

   cvp = cvpartition(num_of_trials, 'KFold',kFold);

   for iFold = 1:kFold

        foldIdx = (iRep-1)*kFold+iFold;

            leftOnsetMTr = StimClass('leftMulti');
            rightOnsetMTr = StimClass('rightMulti');

            leftOnsetMTst = StimClass('leftMulti');
            rightOnsetMTst = StimClass('rightMulti');

            indexMoviesTrainMultiple = indexMoviesTest(cvp.training(iFold)); %%this indexMoviesMultiple variable is how we index trials only correct ones

            indexMoviesTestMultiple = indexMoviesTest(cvp.test(iFold));

            allSMultipleTr = allS(cvp.training(iFold));
            allSMultipleTst = allS(cvp.test(iFold));

            leftMultiIndexTr = indexMoviesTrainMultiple==0;
            rightMultiIndexTr = indexMoviesTrainMultiple==1;

            leftMultiIndexTst = indexMoviesTestMultiple==0;
            rightMultiIndexTst = indexMoviesTestMultiple==1;

            % Only correct trials
            AddStims(leftOnsetMTr, allSMultipleTr(leftMultiIndexTr));
            AddStims(rightOnsetMTr, allSMultipleTr(rightMultiIndexTr));

            AddStims(leftOnsetMTst, allSMultipleTst(leftMultiIndexTst));
            AddStims(rightOnsetMTst, allSMultipleTst(rightMultiIndexTst));

            updateStates(leftOnsetMTr);
            updateStates(rightOnsetMTr);

            updateStates(leftOnsetMTst);
            updateStates(rightOnsetMTst);

            % I prefer this over SetStim so I can control index
            stimTr(1,1) = leftOnsetMTr;
            stimTr(1,2) = rightOnsetMTr;

            stimTst(1,1) = leftOnsetMTst;
            stimTst(1,2) = rightOnsetMTst;

            %dodBPFilt = hmrR_BandpassFilt(dod,0,4); %sudan come back
            dodBPFilt = hmrR_BandpassFilt(dod, 0.01, 0.5);

            dc = hmrR_OD2Conc(dodBPFilt,snirf.probe,[1  1  1]);

            mlActMan = {};
            tIncMan = {};
            tIncAuto{1, 1} = ones(size(snirf.data.dataTimeSeries,1),1);  
            Aaux = [];
            rcMap = [];


% sudan come back
             [~,~,~,dcNewTr,~,~,~,~,~,~,betaSS] = hmrR_GLM_MN(dc,stimTr,snirf.probe,mlActAuto,Aaux,tIncAuto,rcMap,...
                [-2  15],1,1,[1.0 1.0 0.0 0.0 0.0 0.0],20,2,0,0); %could have done this way earlier


            % format beta var here % sudan come back
            ssBeta = betaSS{1}(end,:,:);

            ssBeta_sq = squeeze(ssBeta); % sudan come back

            % actually this is not needed, can use stimTst on dcNewTr. %
            % sudan come back
            dcNewTst = hmrR_ssBeta_CV(snirf.data,dc, snirf.probe, mlActAuto, tIncAuto, squeeze(ssBeta));

            [trials_HbOM_Tr, trials_HbRM_Tr, trials_HbTM_Tr] ...
                = createSingleTrialHRF_MultiOnly_NoSave_sudan(dcNewTr, allSMultipleTr); % sudan come back

            [trials_HbOM_Tst, trials_HbRM_Tst, trials_HbTM_Tst] ...
                = createSingleTrialHRF_MultiOnly_NoSave_sudan(dcNewTst, allSMultipleTst);

            [idxLFEF,idxRFEF] = channelfinder(dcNewTr);
            idxLRFEF = [idxLFEF,idxRFEF];


            only_active = mlActAuto{1}; %measurement list in both wavelengths
            only_active_binary = only_active(:,3);
            mlActAuto_changed = {only_active_binary};
            mlList = mlActAuto_changed{1};
            numChn = length(mlList)/2; %number of channels
            % mlList_LFEF = mlList([idxLFEF idxLFEF+numChn]); %all present
            % mlList_RFEF = mlList([idxRFEF idxRFEF+numChn]); %imbalanced, some absent
            % mlList_LRFEF = mlList([idxLRFEF idxLRFEF+numChn]); %test


            %% dealing singularity by randomly removing channels for test purpose change mlLis
            % Example initialization of 'a' as a 1134x1 matrix
            a = zeros(1134, 1);

            % Set specified indices to 1 in the first half
            a(idxLFEF) = 1;

            % Calculate corresponding indices in the second half
            second_half_indices = idxLFEF + 567; 

            % Set corresponding indices to 1 in the second half
            a(second_half_indices) = 1;

            c = (mlActAuto_changed{1} == 1 & a == 1);

            mlList = c;
            % %mlList = mlActAuto_changed{1};

%% simulated data test

% % Parameters for training trials
% n_trials_tr = 60;
% duration = 17;  % seconds
% n_samples = 153;
% sampling_rate = n_samples / duration;
% 
% % Creating the canonical HRF
% t = 0:(1/sampling_rate):duration;
% hrf = gampdf(t, 8) - 0.5 * gampdf(t, 6);
% hrf = hrf(1:n_samples);
% 
% % Initialize arrays for trials and binary logic
% trials_tr = zeros(n_trials_tr, length(hrf));
% binary_logic_tr = zeros(n_trials_tr, 1);
% 
% % Set noise level
% noise_level = 0.01;
% 
% for i = 1:n_trials_tr
%     % Randomly decide if the HRF goes up or down
%     direction = randi([0 1]) * 2 - 1;  % Generates -1 or 1
%     noisy_hrf = hrf * direction + noise_level * randn(1, length(hrf));
%     trials_tr(i, :) = noisy_hrf;  % Do not multiply by direction again
% 
%     % Record the direction in the binary logic variable
%     binary_logic_tr(i) = direction > 0;
% end
% 
% 
% 
% % Parameters
% % Parameters for training trials
% n_trials_tst = 30;
% duration = 17;  % seconds
% n_samples = 153;
% sampling_rate = n_samples / duration;
% 
% % Creating the canonical HRF
% t = 0:(1/sampling_rate):duration;
% hrf = gampdf(t, 8) - 0.5 * gampdf(t, 6);
% hrf = hrf(1:n_samples);
% 
% % Initialize arrays for trials and binary logic
% trials_tst = zeros(n_trials_tst, length(hrf));
% binary_logic_tst = zeros(n_trials_tst, 1);
% 
% noise_level = 0.01;
% 
% for i = 1:n_trials_tst
%     % Randomly decide if the HRF goes up or down
%     direction = randi([0 1]) * 2 - 1;  % Generates -1 or 1
%     noisy_hrf = hrf * direction + noise_level * randn(1, length(hrf));
%     trials_tst(i, :) = noisy_hrf;  % Do not multiply by direction again
% 
%     % Record the direction in the binary logic variable
%     binary_logic_tst(i) = direction > 0;
% end
% 
% 
% 
% 
% % Assuming 'trials' is your newly generated trials with size 44 x 153
% 
% % Reshape trials to match the second and third dimensions of trials_HbOM_Tr
% reshaped_trials_tr = permute(trials_tr, [3, 2, 1]);  % reshaped_trials will have size 1 x 153 x 44
% 
% trials_HbOM_Tr = zeros(567, 153, 60);
% 
% % Replace each set of trials in trials_HbOM_Tr
% for i = 1:size(trials_HbOM_Tr, 1)
%     trials_HbOM_Tr(i, :, :) = reshaped_trials_tr;
% end
% 
% reshaped_trials_tst = permute(trials_tst, [3, 2, 1]);  % reshaped_trials will have size 1 x 153 x 44
% 
% trials_HbOM_Tst = zeros(567, 153, 30);
% 
% % Replace each set of trials in trials_HbOM_Tr
% for i = 1:size(trials_HbOM_Tst, 1)
%     trials_HbOM_Tst(i, :, :) = reshaped_trials_tst;
% end
% 
% indexMoviesTrainMultiple = binary_logic_tr;
% indexMoviesTestMultiple = binary_logic_tst;

% figure(1);
% for i =1:60
%     plot(linspace(-2,15,153),trials_HbOM_Tr(1,:,i))
%     hold on;
%     xlim([-2 15])
% end


% %% d prime calculation assuming equal variances
% 
% signal_channel_data = squeeze(trials_HbOM_Tr(1, :, indexMoviesTrainMultiple == 1));
% noise_channel_data = squeeze(trials_HbOM_Tr(1, :, indexMoviesTestMultiple == 0));
% 
% n_time_points = size(trials_HbOM_Tr, 2);
% 
% d_prime_values = zeros(n_time_points, 1);
% 
% for t = 1:n_time_points
%     mu_signal = mean(signal_channel_data,2);
%     mu_noise = mean(noise_channel_data,2);
%     sigma_signal = std(signal_channel_data,[],2);
%     sigma_noise = std(noise_channel_data,[],2);
% 
%     d_prime_values(t) = abs((mu_signal(t) - mu_noise(t)))./ sqrt(0.5.* (sigma_signal(t).^2 + sigma_noise(t).^2));
% end
% 
% 
% figure(2);
% plot(linspace(-2,15,153),d_prime_values);
% title("d' vs sample points")
% xlabel("Time in seconds")
% ylabel("d'")
% xlim([-2 15])
% %% d prime rank channels
% 
% % Number of channels and time points
% n_channels = size(trials_HbOM_Tr, 1);
% n_time_points = size(trials_HbOM_Tr, 2);
% 
% % Initialize arrays
% d_prime_max = zeros(n_channels, 1);
% 
% for ch = 1:n_channels
%     %ch = 301;
%     % Extract signal and noise data for the current channel
%     signal_channel_data = squeeze(trials_HbOM_Tr(ch, :, indexMoviesTrainMultiple == 1));
%     noise_channel_data = squeeze(trials_HbOM_Tr(ch, :, indexMoviesTrainMultiple == 0));
% 
%     % Initialize d-prime array for this channel
%     d_prime_values = zeros(n_time_points, 1);
%     d_prime_mean = [];
% 
%     % Calculate d-prime for each time point
%     for t = 1:n_time_points
%         mu_signal = mean(signal_channel_data,2);
%         mu_noise = mean(noise_channel_data,2);
%         sigma_signal = std(signal_channel_data,[],2);
%         sigma_noise = std(noise_channel_data,[],2);
% 
%         d_prime_values(t) = abs((mu_signal(t) - mu_noise(t)))./ sqrt(0.5.* (sigma_signal(t).^2 + sigma_noise(t).^2));
%         d_prime_mean = [d_prime_mean; cumsum(d_prime_values(18+18:18+27))];
%     end
% 
% 
% 
%     % figure(2);
%     % plot(linspace(-2,15,153),d_prime_values);
%     % title("d' vs sample points")
%     % xlabel("Time in seconds")
%     % ylabel("d'")
%     % xlim([-2 15])
% 
% 
%     % Find the maximum d-prime value for this channel
%     d_prime_max(ch) = max(d_prime_mean);
% end
% 
% 
% 
% % Create an index array corresponding to channel numbers
% channel_indices = (1:length(d_prime_max))';
% 
% % Identify and remove NaN values from d_prime_max and the corresponding indices
% valid_indices = ~isnan(d_prime_max);
% d_prime_valid = d_prime_max(valid_indices);
% channel_indices_valid = channel_indices(valid_indices);
% 
% % Sort the valid d-prime values in descending order
% [d_prime_sorted, sort_idx] = sort(d_prime_valid, 'descend');
% 
% % Get the original channel indices corresponding to the sorted d-prime values
% sorted_channel_indices = channel_indices_valid(sort_idx);
% 
% % Extract the top ten d-prime values and their corresponding original channel indices
% top_ten_d_prime = d_prime_sorted(1:min(10, length(d_prime_sorted)));
% top_ten_channels = sorted_channel_indices(1:min(10, length(sorted_channel_indices)));
% 
% % Create a matrix to display channel number and its rank
% top_ten_rankings = [(1:length(top_ten_d_prime)).' top_ten_channels top_ten_d_prime];



%% post drime sorting channels

            % %% dealing singularity by randomly removing channels for test purpose change mlLis
            % % Example initialization of 'a' as a 1134x1 matrix
            % top_ten_channels = [5];
            % 
            % idxdprime = top_ten_channels(1:num_ch_selected);
            % b = zeros(1134, 1);
            % 
            % % Set specified indices to 1 in the first half
            % b(idxdprime) = 1;
            % 
            % % Calculate corresponding indices in the second half
            % second_half_indices = idxdprime + 567; 
            % 
            % % Set corresponding indices to 1 in the second half
            % b(second_half_indices) = 1;
            % 
            % mlList = b;
            % mlList = mlActAuto_changed{1};



%%

            % Train
            % trials_HbOM_Tr_LFEF = trials_HbOM_Tr(idxLFEF,:,:);
            % trials_HbOM_Tr_RFEF = trials_HbOM_Tr(idxRFEF,:,:);
            % trials_HbOM_Tr_LRFEF = trials_HbOM_Tr(idxLRFEF,:,:);
            % 
            % trials_HbRM_Tr_LFEF = trials_HbRM_Tr(idxLFEF,:,:);
            % trials_HbRM_Tr_RFEF = trials_HbRM_Tr(idxRFEF,:,:);
            % trials_HbRM_Tr_LRFEF = trials_HbRM_Tr(idxLRFEF,:,:);
            % 
            % trials_HbTM_Tr_LFEF = trials_HbTM_Tr(idxLFEF,:,:);
            % trials_HbTM_Tr_RFEF = trials_HbTM_Tr(idxRFEF,:,:);
            % trials_HbTM_Tr_LRFEF = trials_HbTM_Tr(idxLRFEF,:,:);
            % 
            % % Test
            % trials_HbOM_Tst_LFEF = trials_HbOM_Tst(idxLFEF,:,:);
            % trials_HbOM_Tst_RFEF = trials_HbOM_Tst(idxRFEF,:,:);
            % trials_HbOM_Tst_LRFEF = trials_HbOM_Tst(idxLRFEF,:,:);
            % 
            % trials_HbRM_Tst_LFEF = trials_HbRM_Tst(idxLFEF,:,:);
            % trials_HbRM_Tst_RFEF = trials_HbRM_Tst(idxRFEF,:,:);
            % trials_HbRM_Tst_LRFEF = trials_HbRM_Tst(idxLRFEF,:,:);
            % 
            % trials_HbTM_Tst_LFEF = trials_HbTM_Tst(idxLFEF,:,:);
            % trials_HbTM_Tst_RFEF = trials_HbTM_Tst(idxRFEF,:,:);
            % trials_HbTM_Tst_LRFEF = trials_HbTM_Tst(idxLRFEF,:,:);

            sbjNum = '1';
            numClasses = 2;

% %%
%             figure(1);
%             for i =1:18
%                 plot(trials_HbOM_Tr(1,:,i))
%                 hold on;
%             end
%%

            % train classifier
            performanceLDALedoitHbO_All(foldIdx,:) = ...
                calcCumSum_DiffTLen_AllChn_CompClassifiers_Final_NoSave_sudan(sbjNum,trials_HbOM_Tr,trials_HbOM_Tst,...
                indexMoviesTrainMultiple,indexMoviesTestMultiple,numClasses,mlList);

            performanceLDALedoitHbR_All(foldIdx,:) = ...
                calcCumSum_DiffTLen_AllChn_CompClassifiers_Final_NoSave_sudan(sbjNum,trials_HbRM_Tr,trials_HbRM_Tst,...
                indexMoviesTrainMultiple,indexMoviesTestMultiple,numClasses,mlList);

            performanceLDALedoitHbT_All(foldIdx,:) = ...
                calcCumSum_DiffTLen_AllChn_CompClassifiers_Final_NoSave_sudan(sbjNum,trials_HbTM_Tr,trials_HbTM_Tst,...
                indexMoviesTrainMultiple,indexMoviesTestMultiple,numClasses,mlList);

            % % Left FEF
            % performanceLDALedoitHbO_LFEF(foldIdx,:) = ...
            %     calcCumSum_DiffStartT_AllChn_CompClassifiers_NoSave_Subset(sbjNum,trials_HbOM_Tr_LFEF,trials_HbOM_Tst_LFEF,...
            %     indexMoviesTrainMultiple,indexMoviesTestMultiple,numClasses,mlList_LFEF);
            % 
            % performanceLDALedoitHbR_LFEF(foldIdx,:) = ...
            %     calcCumSum_DiffStartT_AllChn_CompClassifiers_NoSave_Subset(sbjNum,trials_HbRM_Tr_LFEF,trials_HbRM_Tst_LFEF,...
            %     indexMoviesTrainMultiple,indexMoviesTestMultiple,numClasses,mlList_LFEF);
            % 
            % performanceLDALedoitHbT_LFEF(foldIdx,:) = ...
            %     calcCumSum_DiffStartT_AllChn_CompClassifiers_NoSave_Subset(sbjNum,trials_HbTM_Tr_LFEF,trials_HbTM_Tst_LFEF,...
            %     indexMoviesTrainMultiple,indexMoviesTestMultiple,numClasses,mlList_LFEF);
            % 
            % % Right FEF
            % performanceLDALedoitHbO_RFEF(foldIdx,:) = ...
            %     calcCumSum_DiffStartT_AllChn_CompClassifiers_NoSave_Subset(sbjNum,trials_HbOM_Tr_RFEF,trials_HbOM_Tst_RFEF,...
            %     indexMoviesTrainMultiple,indexMoviesTestMultiple,numClasses,mlList_RFEF);
            % 
            % performanceLDALedoitHbR_RFEF(foldIdx,:) = ...
            %     calcCumSum_DiffStartT_AllChn_CompClassifiers_NoSave_Subset(sbjNum,trials_HbRM_Tr_RFEF,trials_HbRM_Tst_RFEF,...
            %     indexMoviesTrainMultiple,indexMoviesTestMultiple,numClasses,mlList_RFEF);
            % 
            % performanceLDALedoitHbT_RFEF(foldIdx,:) = ...
            %     calcCumSum_DiffStartT_AllChn_CompClassifiers_NoSave_Subset(sbjNum,trials_HbTM_Tr_RFEF,trials_HbTM_Tst_RFEF,...
            %     indexMoviesTrainMultiple,indexMoviesTestMultiple,numClasses,mlList_RFEF);
            % 
            % %Both LRFEF
            % 
            % performanceLDALedoitHbO_LRFEF(foldIdx,:) = ...
            %     calcCumSum_DiffStartT_AllChn_CompClassifiers_NoSave_Subset(sbjNum,trials_HbOM_Tr_LRFEF,trials_HbOM_Tst_LRFEF,...
            %     indexMoviesTrainMultiple,indexMoviesTestMultiple,numClasses,mlList_LRFEF);
            % 
            % performanceLDALedoitHbR_LRFEF(foldIdx,:) = ...
            %     calcCumSum_DiffStartT_AllChn_CompClassifiers_NoSave_Subset(sbjNum,trials_HbRM_Tr_LRFEF,trials_HbRM_Tst_LRFEF,...
            %     indexMoviesTrainMultiple,indexMoviesTestMultiple,numClasses,mlList_LRFEF);
            % 
            % performanceLDALedoitHbT_LRFEF(foldIdx,:) = ...
            %     calcCumSum_DiffStartT_AllChn_CompClassifiers_NoSave_Subset(sbjNum,trials_HbTM_Tr_LRFEF,trials_HbTM_Tst_LRFEF,...
            %     indexMoviesTrainMultiple,indexMoviesTestMultiple,numClasses,mlList_LRFEF);



   end 
end

%% plotting
yLimAxis = [0 1];
    
    performanceLDALedoitHbOCI = zeros(length(timeLen),2);
    performanceLDALedoitHbRCI = zeros(length(timeLen),2);
    performanceLDALedoitHbTCI = zeros(length(timeLen),2);


    % mean
    performanceLDALedoitHbOMean = zeros(length(timeLen),1);
    performanceLDALedoitHbRMean = zeros(length(timeLen),1);
    performanceLDALedoitHbTMean = zeros(length(timeLen),1);
     
    conf = 0.95;

    for i = 1:length(timeLen)

        temp = performanceLDALedoitHbO_All(:,i);
        se = std(temp(:))/sqrt(nrept*kFold);
        performanceLDALedoitHbOMean(i) = mean(temp(:));
        performanceLDALedoitHbOCI(i,:) = calcCITDist(nrept*kFold,se,conf);

        temp = performanceLDALedoitHbR_All(:,i);
        se = std(temp(:))/sqrt(nrept*kFold);
        performanceLDALedoitHbRMean(i) = mean(temp(:));
        performanceLDALedoitHbRCI(i,:) = calcCITDist(nrept*kFold,se,conf);

        temp = performanceLDALedoitHbT_All(:,i);
        se = std(temp(:))/sqrt(nrept*kFold);
        performanceLDALedoitHbTMean(i) = mean(temp(:));
        performanceLDALedoitHbTCI(i,:) = calcCITDist(nrept*kFold,se,conf);

    end

    numSubSets = 4;
    cmap = jet(numSubSets);

    figure('units','normalized','outerposition',[0 0 1 1]); hold on;
    subplot(1,3,1);hold on;

    errorbar(timeLen./fs,performanceLDALedoitHbOMean,performanceLDALedoitHbOCI(:,2),'Color',cmap(1,:));

    ylim(yLimAxis);
    title(sprintf('Sbj %s: Δ[HbO] Multi',num2str(sbjNum)));
    xlabel('Time [s]');ylabel('Accuracy');
    hold off;
    grid on;

    subplot(1,3,2);hold on;

    errorbar(timeLen./fs,performanceLDALedoitHbRMean,performanceLDALedoitHbRCI(:,2),'Color',cmap(1,:));

    ylim(yLimAxis);
    title(sprintf('Sbj %s: Δ[HbR] Multi',num2str(sbjNum)));
    legend({'subset'},'Location','southwest');
    xlabel('Time [s]');ylabel('Accuracy');
    hold off;
    grid on;

    subplot(1,3,3);hold on;

    errorbar(timeLen./fs,performanceLDALedoitHbTMean,performanceLDALedoitHbTCI(:,2),'Color',cmap(1,:));

    ylim(yLimAxis);
    title(sprintf('Sbj %s: Δ[HbT] Multi',num2str(sbjNum)));
    xlabel('Time [s]');ylabel('Accuracy');
    hold off;
    grid on;




%%
function [snirf, time, events] = loadDataAndPreprocess(runPath)
    % Load SNIRF file and preprocess data
    snirf = SnirfLoad([runPath, '.snirf']);
    snirf.data(1).dataTimeSeries = hmrR_PreprocessIntensity_Negative(snirf.data);
    events = snirf.aux(1).dataTimeSeries;
    time = snirf.aux(1).time;
end

function [allS, indexMoviesTest] = processEvents(events, time, indexMoviesTest)
    % Process events to find onsets and filter triggers
    diff_events = diff(events);
    onsets = diff_events > 0.0001;
    t_onset = time(onsets == 1);

    % Filter triggers
    filtered_triggers = t_onset(1);
    previous_trigger_time = t_onset(1);
    for i = 2:length(t_onset)
        current_trigger_time = t_onset(i);
        if (current_trigger_time - previous_trigger_time) >= 3
            filtered_triggers = [filtered_triggers; current_trigger_time];
            previous_trigger_time = current_trigger_time;
        end
    end

    % Ensure indexMoviesTest matches the number of triggers
    minimum_trials = min([length(indexMoviesTest), length(filtered_triggers)]);
    allS = filtered_triggers(1:minimum_trials);
    indexMoviesTest = indexMoviesTest(1:minimum_trials);
end

function [idxLFEF_after,idxRFEF_after] = channelfinder(dcNewTr)
%LFEF_source = [6 6 6 6 6 6 6 6 29 29 29 29 29 29 29 29 26 26 26 26 26 26 26 26 25 25 25 25 25 25 25 25 53 53 53 53 53 53 53 53 27 27 27 27 27 27 27 27];
%LFEF_detector = [17 33 10 12 8 14 15 19 14 16 20 21 19 17 22 15 17 19 21 22 26 24 34 36 34 24 26 28 29 31 35 41 38 36 34 24 35 71 39 40 22 21 23 25 27 28 26 24];

 %LFEF_source = [27 27 27 27 27 27 27 27 27 27 27 27];
 %LFEF_detector = [22 26 21 28 23 27 24 20 25 30 29 19];

 %LFEF_source = [7 7 7 7 7 7 7 7 7 7 7];
 %LFEF_detector = [10 1 5 8 3 7 12 114 118 113 33];

 %LFEF ch 1 from color plot
 %LFEF_source = [29];
 %LFEF_detector = [18];

 %LFEF_source = [6 6 6 6 6 6 6 6 29 29 29 29 29 29];
 %LFEF_detector = [17 33 10 12 8 14 15 19 14 16 20 21 19 17];

 % LFEF_source = [35 35 35 35 35 35 35];
 % LFEF_detector = [43 49 30 45 47 42 32];

 LFEF_source = [30 30 30 30 30 30 30];
 LFEF_detector = [27 30 25 23 28 32 42];

 foundIndices = cell(length(LFEF_source), 1);

% Loop over each source-detector pair
for j = 1:length(LFEF_source)
    givenSourceIndex = LFEF_source(j);
    givenDetectorIndex = LFEF_detector(j);
    matchIndices = [];

    % Inner loop to check each measurement
    for i = 1:3:567*3
        currentSourceIndex = dcNewTr.measurementList(1, i).sourceIndex;
        currentDetectorIndex = dcNewTr.measurementList(1, i).detectorIndex;

        if currentSourceIndex == givenSourceIndex && currentDetectorIndex == givenDetectorIndex
            channelNumber = ceil(i / 3);
            matchIndices = [matchIndices, channelNumber]; % Append the index to the list of matches
        end
    end

    % Store the indices where matches were found for this pair
    foundIndices{j} = matchIndices;
end

% Display results
for j = 1:length(foundIndices)
    if ~isempty(foundIndices{j})
        fprintf('Matches for Source %d and Detector %d found at indices: %s\n', LFEF_source(j), LFEF_detector(j), mat2str(foundIndices{j}));
    else
        fprintf('No matches found for Source %d and Detector %d in the first 567 elements.\n', LFEF_source(j), LFEF_detector(j));
    end
end
%%

%% find right fef channels

 %RFEF_source = [9 9 9 9 9 9 9 9 12 12 12 12 12 12 12 12 10 10 10 10 10 10 10 10 10 10 52 52 52 52 52 52 52 52 51 51 51 51 51 51 11 11 11 11 11 11 11 11];
 %RFEF_detector = [136 139 143 118 119 117 120 134 136 134 132 130 131 135 137 139 129 130 132 134 119 117 120 88 127 88 127 120 117 115 116 81 83 86 116 115 117 81 83 82 125 127 88 86 87 80 129 120];

 %RFEF_source = [9 9 9 9];
 %RFEF_detector = [134 141 136 139];

 %RFEF_source = [54 54 54 54 54 54 54 54 54 54];
 %RFEF_detector = [119 118 113 117 141 115 116 114 3 120];

 % RFEF_source = [54 54 54 54 54 54 54 54 54 54];
 % RFEF_detector = [119 118 113 117 141 115 116 114 3 120];

 RFEF_source = [34];
 RFEF_detector = [54];

 foundIndicesRight = cell(length(RFEF_source), 1);

% Loop over each source-detector pair
for j = 1:length(RFEF_source)
    givenSourceIndex = RFEF_source(j);
    givenDetectorIndex = RFEF_detector(j);
    matchIndices = [];

    % Inner loop to check each measurement
    for i = 1:3:567*3
        currentSourceIndex = dcNewTr.measurementList(1, i).sourceIndex;
        currentDetectorIndex = dcNewTr.measurementList(1, i).detectorIndex;

        if currentSourceIndex == givenSourceIndex && currentDetectorIndex == givenDetectorIndex
            channelNumber = ceil(i / 3);
            matchIndices = [matchIndices, channelNumber];
        end
    end

    % Store the indices where matches were found for this pair
    foundIndicesRight{j} = matchIndices;
end

% Display results
for j = 1:length(foundIndicesRight)
    if ~isempty(foundIndicesRight{j})
        fprintf('Matches for Source %d and Detector %d found at indices: %s\n', RFEF_source(j), RFEF_detector(j), mat2str(foundIndicesRight{j}));
    else
        fprintf('No matches found for Source %d and Detector %d in the first 567 elements.\n', RFEF_source(j), RFEF_detector(j));
    end
end

%%

idxLFEF_after = [foundIndices{:}];
idxRFEF_after = [foundIndicesRight{:}];
end