function WMEAN = wmean(w,X)
% weighted arithmetic mean
% 2017-03-08 tsonne: created
[N,NPAR] = size(X);
sw = sum(w); % sum of weights
WMEAN = nan(1,NPAR); % all weighted means
for a = 1:NPAR
    WMEAN(a) = sum(w.*X(:,a))/sw; % weighted mean
end