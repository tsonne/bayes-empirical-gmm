function S = config_Y2()
% GMM FUNCTION
% used model naming system:
%   - char 1: all start with 'Y'
%   - char 2: functional form ID number
%   - dash '-'
%   - char 4+: indicates prior PDFs: 'up' if all uniform priors, otherwise
%              'npXYZ' where XYZ are coeffecients with normal priors
%   - dash '-'
%   - next 5 char: info about units (log-base, unit, component)
%                  'LD'=decadic log [log10()], 'LN' = natural log [log()],
%                  'c'=[cm/s^2], 'm'=[m/s^2], 'g'=[g],
%                  'GM'=geom.mean, 'RIA'=rot.invariant avg., ...
%                      'LV'=larger value, etc...
%
% Example: 'Y1np124-LDcGM' means model form 1, parameters c1,c2,c4 have
%          normal priors, ground motion unit is in log10([cm/s^2]) of the
%          components geometric mean.
%
% 2020-05-12 Tim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% names (keep it short, will be part of file name)
mod_n = {
    'Y2-up-LDmRIA' % all unconstrained
    'Y2-np4-LDmRIA' % only M-R is constrained
    };
% set how data shall be scaled for processing from raw:
% input raw data expected non-log, cm/s^2,
% some typical conventions of GMMs are available:
% LOGT can be 'LD'=decadic log [log10()], 'LN' = natural log [log()],
% UNIT can be 'c'=[cm/s^2], 'm'=[m/s^2], 'g'=[g]
LOGT = 'LD';
UNIT = 'm';
% function definition strings ready for use as Matlab functions
sf_x = {
    'c(1)+c(2)*M+(c(3)+c(4)*M).*log10(sqrt(c(5)^2+R.^2))+c(6)*S';
    'c(1)+c(2)*M+(c(3)+c(4)*M).*log10(sqrt(c(5)^2+R.^2))+c(6)*S';
    };
% original paper coeff values:
PCS = coeff_Am05(); % Periods, Coefficients, Sigma
% chosen periods:
per = [0 0.2 0.3 1];
%per = 'all';
%per = [0.00 0.05 0.10 0.20 0.30 0.50 0.75 1.00 2.00]; 
% global optimization starting values for each parameter incl. sigmas
X0 = [
       0,   0.30,   -1, 0.314,     5,  0.35, 0.05,  0.2];
% optimization/prior lower boundaries (keep as tight as possible)
lowB = [
     -10,     -2,   -5,    -1,     0,    -1,    0,    0];
% optimization/prior upper boundaries (keep as tight as possible)
uppB = [
      10,      2,    2,     1,    12,     2,  0.5,  0.7];
% Prior density function types for each parameter (u=unif,n=norm)
PriorType = {
    'uuuuuuuu'
    'uuuNuuuu'
    };
RatioMean2SD = 10; % ratio of mean to standard dev. of normal prior
%
S = struct();
S.mod_n = mod_n;
S.sf_x = sf_x;
S.PCS = PCS;
S.per = per;
S.X0 = X0;
S.lowB = lowB;
S.uppB = uppB;
S.PriorType = PriorType;
S.RatioMean2SD = RatioMean2SD;
S.LOGT = LOGT;
S.UNIT = UNIT;
