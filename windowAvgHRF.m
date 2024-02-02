function yavg_conc = windowAvgHRF(y_conc, tHRF, trange)

    lst = find(tHRF>=trange(1) & tHRF<=trange(2));
    yavg_conc = mean(y_conc(lst,:),1); % take average 

end