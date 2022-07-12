% calc random event-effect (Abrahamson&Youngs1992)
function [etai,mag] =  random_effect_etai(yij,mij,vNS,MRS,tau2,sigma2)
% INPUT:   intensity measures in log10 units!)
% yij = log10(PSAobs) (row: station-event, column: frequency)
% mij = log10(PSAsim)
% vNS = vector with number of stations per event
% MRS = matrix [Mw,R,soil] for each event-station
% tau2 = log10-based intra-event covariance
% sigma2 = log10-based inter-event variance
%
% 2016-03-27 tsonne: created
nf = size(yij,2);
M  = length(vNS);
etai  = nan(M,1);
mag   = nan(M,1);
for a = 1:M
    i2 = sum(vNS(1:a));
    i1 = i2 -vNS(a) +1;
    ii = i1:i2;
    ni = vNS(a)*nf;
    mag(a) = MRS(i1,1);
    etai(a) = sum(sum(yij(ii,:)-mij(ii,:)))*tau2/(ni*tau2+sigma2);
end