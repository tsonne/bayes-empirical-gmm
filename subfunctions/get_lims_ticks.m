% Give data and get appropriate plot limits and tick positions:
function [pl,ticks] = get_lims_ticks(px,lb,ub,FREE,AXIS_LIM,BURN,NB,perc)
% INPUT:
%  px: data matrix (1 col/par) or struct px(npar).rnd(nt,nc)
% optional:
%  lb: lower boundary vector
%  ub: upper boundary vector
%  FREE: (vector) which variables are free (only used by SBMRVTMCMCvX.m)
%  AXIS_LIM: if 'tight' axis is min,max(data), 
%           if 'robust' clip at 1 to 99 percentile, 
%           else use given lb,ub
%  BURN: 'yes' (use burn-in data) or 'no' (use only post-burn-in data)
%  NB: burn-in length
%  perc: limit percentiles for 'robust' setting, DEF:[0.01 099]
%
% OUTPUT:
%  pl: matrix (npar,2) low and high limit for each parameter|column
%  ticks: cell with suggested tick position vector for each par|col
%
% 2016-03-27 tsonne: Outsourced routine to save space in main script.
%                    Rewrite to make it more useful for general cases.
%
cpar = {'lb','ub','FREE','AXIS_LIM','BURN','NB','perc'};
% check all input variables, assign default if empty
for a=1:size(cpar,1)
    if ~exist(cpar{a},'var') || isempty(eval(cpar{a}))
        eval([cpar{a} '= {};']);
    end
end
if ~isstruct(px) % assume matrix instead! one column per parameter
    [NT,NPAR] = size(px);
    xx = struct('rnd',nan(NT,1));
    xx = repmat(xx,1,NPAR);
    for a = 1:NPAR
        xx(a).rnd = px(:,a);
    end
    px = xx;
else % assume my typical px.rnd struct from MCMC runs
    NPAR = length(px);
end
if ~any( strcmpi(AXIS_LIM,{'tight','robust'}) ) % use low|upp bounds
    if isempty(lb) || isempty(ub) % if not given, switch to min|max
        AXIS_LIM = 'tight';
        pl = nan(NPAR,2);
    else
        pl = [lb' ub']; % parameter plot limits
    end
end
% avoid error if burn wanted but no number given:
if isempty(NB), BURN='yes'; end
% set limits according to setting:
if strcmpi(AXIS_LIM,'tight') % min,max data
    for a = 1:NPAR
        if any(strcmpi(BURN,{'n','no'})) % only data after burn-in
            so = px(a).rnd(NB+1:end,:);
            pl(a,:) = [.97*min(so(:)) 1.03*max(so(:))];
        else % use all data
            pl(a,:) = [.97*min(px(a).rnd(:)) 1.03*max(px(a).rnd(:))];
        end
    end
elseif strcmpi(AXIS_LIM,'robust') % data from 1 to 99 percentile
    if isempty(perc) || numel(perc)~=2
        perc = [0.01 0.99];
    end
    for a = 1:NPAR
        if any(strcmpi(BURN,{'n','no'})) % only data after burn-in
            so = px(a).rnd(NB+1:end,:);
            so = sort(so(:));
            nn = sum(~isnan(so));
            ii = round(nn*perc);
            pl(a,:) = [so(ii(1)) so(ii(2))];
        else
            so = sort(px(a).rnd(:)); % use all data
            ii = round(numel(so)*perc);
            pl(a,:) = [so(ii(1)) so(ii(2))];
        end
    end
end
ticks = cell(NPAR,1);
for a = 1:NPAR
    rx = pl(a,2)-pl(a,1);
    x  = 10^floor(log10(rx));
    if rx/x>7, mt=4*x;
    elseif rx/x>4, mt=2*x;
    elseif rx/x>2, mt=x;
    elseif rx/x>1, mt=x/2;
    else mt=x/4;
    end
    ticks(a) = {floor(pl(a,1)/mt)*mt : mt : ceil(pl(a,2)/mt)*mt};
end