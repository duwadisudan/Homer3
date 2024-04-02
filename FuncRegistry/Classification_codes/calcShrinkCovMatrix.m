% Compare against clsutil_shrinkage for validation. Done.
% Pick Bessel's correction factor n-1
% X is p x n matrix.
% Use pooled covariance matrix instead of covariance matrix computed on
% entire dataset.
% Try different target matrices.

function [covHat,gamma] = calcShrinkCovMatrix(X,S)
    [p, n] = size(X);
    
    if p==1
        covHat = var(X);
        gamma = 0;
        return;
    end
    
    % Zero-centering the data first.
    % Assume X is iid
    Xn     = X - repmat(mean(X,2), [1 n]);
    if nargin == 1
        S      = Xn*Xn'/(n-1);
    end

    % Need to use pooled covariance matrix instead.
    
    % asymptotic mu
    mu = mean(diag(S));
    % asymptotic delta
    temp = S-mu*eye(p);
    deltaSquared = mean(diag(temp*temp'));
    % asymptotic beta
    betaSquared = 0;
    for i = 1:n
        temp = X(:,i)*X(:,i)'-S;
        betaSquared = betaSquared + mean(diag(temp*temp'));
    end
    betaSquared = betaSquared/((n-1)^2);
    betaSquaredCeil = max(0,min(betaSquared,deltaSquared));

    gamma = betaSquaredCeil/deltaSquared;

    covHat = gamma*mu*eye(p) + (1-gamma)*S;

end