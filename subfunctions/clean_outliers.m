function [px, mask] = clean_outliers(px,NB)
% Partly reassess given MCMC result struct by removing outliers.
% Assume outlier when:  abs(x_i - MEDIAN(X)) > 4*1.4826*MAD(X).
% Update struct: mean, std, perc (vector, scalar, log),
%                set outliers to NaN in px.rnd
% No change to be expected for other values (Gelman-stat, autocorr), but
% recalculated anyway. This way, unanalyzed struct can be filled up.
%
% Note:
% It is possible in this function that a NAN might be returned as mode, if
% the highest probability (px.pdf) is represented by an assumed outlier.
% This should be very unlikely and might indicate deep trouble.
%
% 2016-03-29 tsonne: created
% 2016-10-20 tsonne: added output 'mask' to indicate which chain samples
%                    were not outliers (1=true), and outliers (0=false).
NP = length(px);
[NT,NC] = size(px(1).rnd);
NL = NT - NB;
maxlag = size(px(1).autocorr,1)-1;
mask = true(NL,NC); % [2016-10-20]
% statistics: (mean, std, percentiles: 2.5, 25, 50, 75, 97.5 %)
% new: logarithm stats (mean, std, mode of log(px)) & autocorr
for a = 1:NP
    % vectors: values for each chain:
    M1 = px(a).rnd(NB+1:end,:);
    VL = M1(:);
    mm = median(VL); % robust mean substitute
    ss = 1.4826*mad(VL,1); % robust standard deviation substitute
    MI = abs(M1-mm)>4.0*ss; % find values beyond 4 sigma (0.006334%)
    ML = M1;
    ML(MI) = nan;
    mask(MI) = false; % [2016-10-20]
    % outliers are now NaNs, continue...
    px(a).mean = nanmean(ML); % 1 x NC
    px(a).std  = nanstd(ML); % 1 x NC
    Mso = sort(ML); % NL x NC
    for b = 1:NC
        nln = sum(~isnan(Mso(:,b)));
        px(a).perc(:,b) = [  Mso(ceil(nln*0.025),b);
            Mso(ceil(nln*0.160),b);
            Mso(ceil(nln*0.500),b);
            Mso(ceil(nln*0.840),b);
            Mso(ceil(nln*0.975),b)]; % 5 x NC
    end
    ppm = px(a).pdf(NB+1:end,:);
    [~,im] = max(ppm);
    im = sub2ind(size(ppm),im,1:NC);
    px(a).mode = ML(im);
    % scalars: values for all data together:
    VL = reshape(ML,NL*NC,1);
    px(a).Tmean = nanmean(VL);
    px(a).Tstd  = nanstd(VL);
    Vso = sort(VL); % NL*NC x 1
    nln = sum(~isnan(Vso));
    px(a).Tperc = [Vso(ceil(nln*0.025)); % -2 sigma
        Vso(ceil(nln*0.160)); % -1 sigma
        Vso(ceil(nln*0.500)); % median
        Vso(ceil(nln*0.840)); % +1 sigma
        Vso(ceil(nln*0.975))]; % +2 sigma % 5 x 1
    [~,im] = max(reshape(ppm,NL*NC,1));
    px(a).Tmode = VL(im);
    dmat = diff(M1);
    px(a).acr = sum(sum(abs(dmat)>1e-16))/numel(dmat); % acceptance
    % Gelman's parameters:
    if NC>1
        [px(a).R,px(a).neff,px(a).lag1_corr,px(a).acr] = gpar(M1);
    end
    
    % log stats:
    ln_px = log(ML);
    px(a).ln_mean  = nanmean(ln_px); % 1 x NC
    px(a).ln_std   = nanstd(ln_px); % 1 x NC
    px(a).ln_perc  = log(px(a).perc); % 5 x NC
    px(a).Tln_mean = nanmean(reshape(ln_px,NL*NC,1)); % 1x1
    px(a).Tln_std  = nanstd(reshape(ln_px,NL*NC,1)); % 1x1
    px(a).Tln_perc = log(px(a).Tperc); % 5x1
    % autocorrelation (lag1 to lag50):
    if NC>1
        px(a).autocorr = zeros(maxlag+1,1);
        for b=1:NC
            pxc = xcorr(px(a).rnd(NB+1:end,b) - px(a).mean(b), ...
                maxlag, 'coeff');
            px(a).autocorr = px(a).autocorr + pxc(maxlag+1:end)/NC;
        end
    end
    px(a).rnd(NB+1:end,:) = ML; % update struct: rnd with nan for outliers
end