% Based on student's t distribution
% population mean and SD are unknown
% estimated standard error using n-1 correction factor
% se = std(x)/sqrt(n);
% This is a bad way of calculating confidence interval but is the standard
% among medical journals.
% Can try more accurate method using paper "Cross-validation: what does it
% estimate and how well does it do it?" by Bates, Hastie & Tibshirani 2021.
% Summarized here: https://towardsdatascience.com/your-cross-validation-error-confidence-intervals-are-wrong-heres-how-to-fix-them-abbfe28d390
% conf = 0.95

% n and se can be vectors.

function ci = calcCITDist(n,se,conf)
    %z = norminv((1-conf)/2);
    
    %ebm = abs(z)*SD/sqrt(n);
    nSz = length(n);
    ci = zeros(nSz,2);
    
    for i = 1:nSz
        dof = n(i)-1;
        %pLo = (1-conf)/2.*ones(1,nSz);
        %pUp = 1-(1-conf)/2.*ones(1,nSz);

        pLo = (1-conf)/2;
        pUp = 1-(1-conf)/2;

        t = tinv([pLo pUp], dof);

        ci(i,:) = se(i).*t;
    end
end