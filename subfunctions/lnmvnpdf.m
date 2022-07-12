% Log-Likelihood function: Normal with given cov matrix & sim data
function lnP = lnmvnpdf(X, MU, COV)
% INPUT:
%   
%   X    = samples, i.e. data vector (length n)
%   MU   = location, i.e. prediction vector (length n)
%   COV  = covariance matrix (n x n)
% OUTPUT:
%   lnP   = ln(likelihood)
%
% 2016-07-07 tsonne: created
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loglikelihood:
D = X(:) - MU(:);
R = chol(COV); % upper triangular matrix by Cholesky factorization of cov.
logdetC = 2*sum(log(diag(R))); % log(det(C))
Ld = (R\(R'\D))'*D; % (obs-sim)' * C^-1 * (obs-sim)
lnP = -0.5*(logdetC + Ld);