function OD = conc2OD(conc, SD)

    mlSD = SD.MeasList;

    E = GetExtinctions(SD.Lambda);
    E = E(:,1:2) / 10; % convert from /cm to /mm
    nTpts = size(conc, 1); 
    ppf = [1,1];
    lst = find( mlSD(:,4)==1 );
    
    OD = zeros(size(conc));
    for idx = 1:length(lst)
    
        idx1 = lst(idx); % index for first wavelength
        idx2 = find( mlSD(:,4)>1 & mlSD(:,1)==mlSD(idx1,1) & mlSD(:,2)==mlSD(idx1,2) ); % index for second wavelength 
    
        rho = norm(SD.SrcPos3D(mlSD(idx1,1),:)-SD.DetPos3D(mlSD(idx1,2),:)); % SD separation
        OD(:,[idx1 idx2']) = (E * conc(:,[idx1, idx2])')' .* (ones(nTpts,1)*rho*ppf);     
    end

end