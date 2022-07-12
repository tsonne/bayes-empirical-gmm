% log of Prior probability functions: mix of uniform and normal
function lnPr = lnPriorNU(x,PT,C1,C2)
%   x   = parameter value vector
%   PT  = prior type char string (each parameter has either 'u' or 'n')
%   C1  = prior coefficient 1 vector (lowB if unif, mean if norm)
%   C2  = prior coefficient 2 vector (uppB if unif, sigma if norm)
%
% Example:
%   x  = [0.23 1.947 3.87 15.2]';
%   PT = 'unuu';  C1 = [0 2 -4  0]';  C2 = [3 1 10 50]';
%   --> -0.0014
%   PT = 'unuu';  C1 = [1 2 -4  0]';  C2 = [3 1 10 50]'; % C1(1) changed
%   --> -Inf
% 2019-02-26 Tim: Created
lp = nan(numel(x),1); % each log(prior)
for a=1:numel(x)
    if PT(a)=='u' % uniform
        if x(a)>=C1(a) && x(a)<=C2(a)
            lp(a) = 0;   % p=1, log(p)=0
        else
            lnPr = -inf; % p=0, log(p)=-Inf (rejected)
            return
        end
    else % normal (constant terms dropped)
        lp(a) = -0.5*((x(a)-C1(a))/C2(a))^2;
    end
end
lnPr = sum(lp);