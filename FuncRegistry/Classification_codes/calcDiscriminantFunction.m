% Calculate linear discriminant function for a given class.
% Required covariance matrix, mean and class prior (can use empirical)
% Is sphering the data necessary if using Ledoit's shrinkage method?
% Also, deal with ties.
% x is p x n vector.
% mu is p x 1 vector.
% sigma is p x p vector. It is symmetric and positive semi-definite.
% discriminantScore is 1 x n vector.
% If X is matrix, more computational efficient to loop.

function discriminantScore = calcDiscriminantFunction(X,sigma,mu,prior)

    % V is p x p matrix, D is diagonal p x p matrix.
    % No need to zero-center the data, just sphere the data.
    % This is computationally more efficient than inv for symmetric matrix.
    [V,D] = eig(sigma);

    invD = 1./D;
    invD(D==0)=0;
    
    [p,n] = size(X);
    
    discriminantScore = zeros(1,n);
    
    for i = 1:n
        % the weight vector normal to the discriminant (separating) hyperplane is 
        % w proportional to (sigma0 + sigma1)^-1 * (mu1-mu0)
        % The decision rule maximizing the posterior probability is the
        % same as the decision rule finding the weight vector w that 
        % maximize the distance between two classes as can be seen in the
        % original discriminant function.
        discriminantScore(1,i) = -1/2*log(det(sigma))-1/2*((V'*(X(:,i)-mu))'*invD*(V'*(X(:,i)-mu)))+log(prior);
    end
    
end