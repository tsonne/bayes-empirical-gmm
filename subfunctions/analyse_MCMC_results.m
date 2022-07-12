% Analysis of MCMC Results (Statistics)
function px = analyse_MCMC_results(px,NPAR,NC,NB,NL,PV,PP,maxlag)
% INPUT:
%  px   = struct px: numel(px) = NPAR,
%                    size(px(a).rnd) = [NB+NL NC],
%                    size(px(a).pdf) = [NB+NL NC].
%  NPAR = number of variables
%  NC   = number of chains
%  NB   = number of burn-in samples per chain
%  NL   = number of after-burn-in samples per chain
%  PV   = input chain cell: numel(PV) = NC, size(PV{a}) = [NB+NL NPAR].
%  PP   = input post.prob. cell: numel(PP) = NC, size(PP{a}) = [NB+NL 1].
%  maxlag = max. lag to calculate autocorrelation for each chain
%
% OUTPUT:
%  px   = struct containing chain samples, probabilities and all statistics
%         for each chain (mean,mode,stdev,percentiles,autocorr,acc.rate)
%         and for all chains combined (only using after-burn-in values)
%
% 2016-11-22 tsonne: created
%                    This function has been copied from script to script,
%                    now it is finally a function for easy insertion.
% 2016-12-05 tsonne: only read probability values from first px element,
%                    instead of going through all parameters for the same
%                    ln_posterior values (values stored only in first par
%                    to save space).
for a =  1:NPAR
    for c = 1:NC % gather chains from cells to matrix
        px(a).rnd(:,c) = PV{c}(:,a);
        if a==1, px(a).pdf(:,c) = PP{c}(:); end % ln_posterior
    end
end
% statistics: (mean, std, percentiles: 2.5, 25, 50, 75, 97.5 %)
% new: logarithm stats (mean, std, mode of log(px)) & autocorr
ppm = px(1).pdf(NB+1:end,:);
[~,im] = max(ppm);
im = sub2ind(size(ppm),im,1:NC); % mode indices per chain
[~,imT] = max(reshape(ppm,NL*NC,1)); % total mode index
for a = 1:NPAR
    % vectors: values for each chain:
    pmat = px(a).rnd(NB+1:end,:);
    px(a).mean = mean(pmat); % 1 x NC
    px(a).std = std(pmat); % 1 x NC
    p_sort = sort(pmat); % NL x NC
    px(a).perc = [  p_sort(ceil(NL*0.025),:);
        p_sort(ceil(NL*0.160),:);
        p_sort(ceil(NL*0.500),:);
        p_sort(ceil(NL*0.840),:);
        p_sort(ceil(NL*0.975),:)]; % 5 x NC
    px(a).mode = pmat(im); % mode
    % scalars: values for all data together:
    pvec = reshape(pmat,NL*NC,1);
    px(a).Tmean = mean(pvec);
    px(a).Tstd = std(pvec);
    p_sort = sort(pvec); % NL*NC x 1
    px(a).Tperc = [p_sort(ceil(NL*NC*0.025)); % -2 sigma
        p_sort(ceil(NL*NC*0.160)); % -1 sigma
        p_sort(ceil(NL*NC*0.500)); % median
        p_sort(ceil(NL*NC*0.840)); % +1 sigma
        p_sort(ceil(NL*NC*0.975))]; % +2 sigma % 5 x 1
    px(a).Tmode = pvec(imT);
    dmat = diff(pmat);
    px(a).acr = sum(sum(abs(dmat)>1e-16))/numel(dmat); % acceptance
    % Gelman's parameters:
    if NC>1
        [px(a).R,px(a).neff,px(a).lag1_corr,px(a).acr] = gpar(pmat);
    end
    % log stats:
    ln_px = log(pmat);
    px(a).ln_mean = mean(ln_px); % 1 x NC
    px(a).ln_std = std(ln_px); % 1 x NC
    px(a).ln_perc = log(px(a).perc); % 5 x NC
    px(a).Tln_mean = mean(reshape(ln_px,NL*NC,1)); % 1x1
    px(a).Tln_std = std(reshape(ln_px,NL*NC,1)); % 1x1
    px(a).Tln_perc = log(px(a).Tperc); % 5x1
    % autocorrelation (lag1 to lag50):
    if NC>1
        for b=1:NC
            pxc = xcorr(px(a).rnd(NB+1:end,b) - px(a).mean(b), ...
                maxlag, 'coeff');
            px(a).autocorr = px(a).autocorr + pxc(maxlag+1:end)/NC;
        end
    end
end