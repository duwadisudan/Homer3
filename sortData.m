function ysorted = sortData(SD, data_yavg, nTask)

    mlSD = SD.MeasList;
    
    nMeas = length(mlSD);
    nChan = nMeas/2;
    
    ml_glm = data_yavg.measurementList(1:nChan*3);
    indHbO = 1:3:length(ml_glm);
    indHbR = 2:3:length(ml_glm);
    
    yUnsorted = data_yavg.dataTimeSeries; % HRFs
    for t = 1:nTask
        yUnsorted_HbO(:,:,t) = yUnsorted(:,indHbO + nChan*3*(t-1));
        yUnsorted_HbR(:,:,t) = yUnsorted(:,indHbR + nChan*3*(t-1));
    end
    mlMatrix = zeros(nChan,2);
    i = 1;
    for m = 1:length(ml_glm)
        
        if ismember(m, indHbO)
            mlMatrix(i,1) = ml_glm(m).sourceIndex;
            mlMatrix(i,2) = ml_glm(m).detectorIndex;        
            i = i+1;
        end
    end
   % mlMatrix = mlMatrix(1:length(mlMatrix),:);
    
   %%
    ysorted_HbO = zeros(size(yUnsorted_HbO));
    ysorted_HbR = zeros(size(yUnsorted_HbR));

    for iM=1:nChan
    
        % get source and detector index in SD
        iS = mlSD(iM,1);
        iD = mlSD(iM,2);
    
        % find corresponding row in mlMatrix 
        ind = find(ismember(mlMatrix, [iS,iD], 'rows'));
    
        % if isempty(ind)
        % 
        %     ysorted_HbO(:,iM) = 0;
        %     ysorted_HbR(:,iM) = 0;
        % 
        % else

        % place 
            ysorted_HbO(:,iM,:) = yUnsorted_HbO(:,ind,:);
            ysorted_HbR(:,iM,:) = yUnsorted_HbR(:,ind,:);
        % end
    
    end

    ysorted = [ysorted_HbO, ysorted_HbR];
        
end

