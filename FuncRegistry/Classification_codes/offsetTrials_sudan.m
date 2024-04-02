% trials is channels x time x trials
function trials = offsetTrials_sudan(trials,fs)


% Offset
for i1 = 1:size(trials,1)
    for i2 = 1:size(trials,3)
        baseline = mean(trials(i1,1:round(2*fs),i2));
        trials(i1,:,i2) = trials(i1,:,i2) - baseline;
    end
end

end