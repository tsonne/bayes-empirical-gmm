function lnP = ln_re_mvnpdf(X, MU, n, vSIGMA)
% Log-Likelihood function: Normal pdf, given obs&sim data & [tau sigma].
% Only for random effects model, assuming that the covariance matrix
% is a block-diagonal matrix according to the data vector which is
% sorted by event, such that each block corresponds to one event.
% Only one single frequency for all events is assumed.
% 
% INPUT:  
%   X    = samples, i.e. data vector (length n)
%   MU   = location, i.e. prediction vector (length n)
%   n    = number of stations for each event (vector)
%   vSIGMA = [sigma_1 sigma_2] % first event cov, then diagonal sigma
% OUTPUT:
%   lnP   = ln(likelihood)
%
% 2017-03-09 tsonne: created
% 2019-02-22 tsonne: added isreal() check, fixed issue: sim shall be real.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREPARATION:
if ~isreal(MU) % check if imaginary parts exist in simulated values
    lnP = -inf; % return -inf loglikelihood to get rejected
    return
end
D = X(:) - MU(:); % residuals
N = sum(n);    % number of data points
M = numel(n);  % number of events
T = vSIGMA(1); % tau
S = vSIGMA(2); % sigma
S2 = S*S;
T2 = T*T;
S2nT2 = S2 + n*T2;
% LOG_DETERMINANT:  log(det(COV))
LogDetC = 2*(N-M)*log(S) + sum(log(S2nT2));
% MAHALANOBIS DISTANCE:  (obs-sim)' * C^-1 * (obs-sim)
MD = 0;
D2=D.*D;
p = cumsum([0;n]);
for a = 1:M
    S = 0;
    for i =  p(a)+1 : p(a+1)-1
        for j = i+1 : p(a+1)
            S = S + D(i)*D(j);
        end
    end
    B = - T2/S2/S2nT2(a);
    A = 1/S2 + B;
    MD = MD + ( A*sum(D2(p(a)+1:p(a+1))) + 2*B*S );
end
% LOG_LIKELIHOOD: log(2*pi) as const
lnP = -0.5*(N*1.837877066409345 + LogDetC + MD);