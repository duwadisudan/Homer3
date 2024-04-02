function data_ynew = sortData_ynew(SD, data_ynew, ml)

    mlSD = SD.MeasList;
    mlSD = mlSD(1:567, 1:2);

    % size(ml) = [567, 2] and size(mlSD) = [567, 2]
    % data_ynew.DataTimeSeries is 13097x3x567
    
    % Initialize an array to hold the new order indices
    newOrderIndices = zeros(size(mlSD, 1), 1);
    
    % Find the index in ml for each row in mlSD
    for i = 1:size(mlSD, 1)
        % Find the row in ml that matches the pair in mlSD
        for j = 1:size(ml, 1)
            if isequal(mlSD(i, :), ml(j, :))
                newOrderIndices(i) = j;
                break;
            end
        end
    end
    
    % Reorder the third dimension of data_ynew.DataTimeSeries
    reorderedData = data_ynew.dataTimeSeries(:, :, newOrderIndices);
     
    % Set data vectors for the dc-parallel data
    data_ynew.SetDataTimeSeries(reorderedData);
        
end

