% Option 1: No GLM, just use SS coefficients = 1 to get individual trials
    % In this option, no CV is performed here. 
% Option 2: Use GLM only for SS coefficients since it's the same for all
% conditions. However, require building design matrix. Also not 100%
% independent of class conditions, still need to know trials labels.
    % In this option, CV is performed here.
% Option 3: Perform GLM on single trial using Gaussian basis and use HRF as
% input to classification.
% 4/5/2022: updated to work with HbT GLM fitting.

% TODO:
    % Specify beta weights
        % betaSS is nChn x 2
        % default all ones
    % Set beta weights to 0 if noisy SS channels.


% Note: data_raw is raw data. data_y is conc data.
function data_ynew = hmrR_ssBeta_CV_sudan(data_raw,data_y, probe, mlActAuto, tIncAuto,...
    betaSS)

    for iBlk=1:length(data_y)

        y      = data_y(iBlk).GetDataTimeSeries('reshape');
        y_raw  = data_raw(iBlk).GetDataTimeSeries('reshape');
        ml     = data_y(iBlk).GetMeasListSrcDetPairs();
        t      = data_y(iBlk).GetTime();
        data_ynew(iBlk)    = DataClass(data_y(iBlk));
        tInc = tIncAuto{iBlk};
        
        %rhoSD_ssThresh = 60;
        rhoSD_ssThresh = 20;
        
        SrcPos = probe.GetSrcPos();
        DetPos = probe.GetDetPos();
        mlAct = mlActAuto{iBlk};

        mlAct = mlAct(:,3);

        lst = 1:size(ml,1);
        rhoSD = zeros(length(lst),1);
        posM = zeros(length(lst),3);
        for iML = 1:length(lst)
            rhoSD(iML) = sum((SrcPos(ml(lst(iML),1),:) - DetPos(ml(lst(iML),2),:)).^2).^0.5;
            posM(iML,:) = (SrcPos(ml(lst(iML),1),:) + DetPos(ml(lst(iML),2),:)) / 2;
        end

        lstSS = lst(find(rhoSD<=rhoSD_ssThresh & mlAct(lst)==1));
        
        y_mean_ss = mean(y(:,:,lstSS),3); %sudan this is where mean of ss signal is
        y_ss = y(:,:,lstSS);
        
        nT = length(t);
        nCh = size(y,3);
        ynew    = zeros(nT,3,nCh);
        
        if ~exist('betaSS','var')
            betaSS = ones(nCh,3);
        end
        
        % Check noisy channels

        %% Meryem's method
    
        for conc = 1:3

            lstML = logical(mlAct(1:length(mlAct)/2));
        
            ytmp = squeeze(y(:,conc,lstML));

            y_mean_ss_tmp = y_mean_ss(:,conc);

            beta_tmp = betaSS(lstML,conc);

            ynew(:,conc,lstML) = ytmp - (y_mean_ss_tmp*beta_tmp');

        end
        
%% matthew's method
        % for conc = 1:3
        % 
        %     lstML = logical(mlAct(1:length(mlAct)/2));
        % 
        %     ytmp = squeeze(y(:,conc,lstML));
        % 
        %     beta_tmp = betaSS(lstSS,conc);
        % 
        %     y_ss_tmp = squeeze(y_ss(:,conc,:));
        % 
        %     ynew(:,conc,lstML) = ytmp - (y_ss_tmp*beta_tmp);
        % 
        % end
        %ynew(:,3,:) = ynew(:,1,:) + ynew(:,2,:);
        
        data_ynew(iBlk).SetDataTimeSeries(ynew);
        
    end
    
end