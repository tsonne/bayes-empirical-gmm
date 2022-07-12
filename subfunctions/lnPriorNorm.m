% log of Normal Prior probability function
function lnPr = lnPriorNorm(x,Mu,Cov)
%   x      = parameter values
%   Mu     = normal prior parameter mean value column vector
%   Cov    = normal prior parameter covariance matrix
%
% 2016-06-30 tsonne
%  - speedup by assuming only diagonal matrix coaviance?
%
% check if prior covariance is diagonal matrix: (i.e. indep.par.)
% ISDIAG = sum(diag(Cov).^2)==sum(sum(Cov.^2));
% if ISDIAG
% log Prior (Normal): (if prior covariance is diagonal)
lnPr = -0.5*sum((x(:)-Mu(:)).^2 ./diag(Cov));
% else:
% ddd = (x - Mu)';
% Cho = chol(Cov);
% lnPr = -0.5*(Cho\(Cho'\ddd))'*ddd;
% end