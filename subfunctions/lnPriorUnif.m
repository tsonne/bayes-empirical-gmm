% log of Uniform Prior probability function
function lnPr = lnPriorUnif(x,lowB,uppB)
%   x     = parameter values
%   lowB  = lower boundaries for x
%   uppB  = upper boundaries for x
%
% 2016-06-30 tsonne
%
% check whether parameters out of bounds:
if sum(x>=lowB)+sum(x<=uppB)~=2*numel(x)
    lnPr = -inf; % return -inf loglikelihood to get rejected
else
    lnPr = 0; % log(1), probability = 1
end
%