function [WCOV,WMEAN,WSD] = wcov(w,X)
% weighted covariance
% w - weight vector
% X - value matrix, each column is one parameter
%
% weights are expected to be probability values, but maybe not necessary
%
% 2017-03-08 tsonne: created
[N,NPAR] = size(X);
assert(NPAR>1,'Require at least 2 columns in second input term!');
sw = sum(w); % sum of weights
WMEAN = wmean(w,X); % weighted means
WSD = nan(1,NPAR);
WCOV = nan(NPAR,NPAR);
for a = 1:NPAR
    da = X(:,a)-WMEAN(a);
    WSD(a) = sqrt(sum(w.*(da.*da))/(sw-1)); % weighted st.dev.
    for b = 1:a
        WCOV(a,b) = da'*((X(:,b)-WMEAN(b)).*w)./(sw-1);
        if (a ~= b)
            WCOV(b,a) = WCOV(a,b); % weighted cov.
        end
    end
end
end