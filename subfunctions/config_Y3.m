function S = config_Y3()
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
    'Y3-up-LNgRIA' % all unconstrained
    'Y3-np45-LNgRIA' % here, M-R1 and M-R2 are constrained (M-saturation)
    };
% set how data shall be scaled for processing from raw:
% input raw data expected non-log, cm/s^2,
% some typical conventions of GMMs are available:
% LOGT can be 'LD'=decadic log [log10()], 'LN' = natural log [log()],
% UNIT can be 'c'=[cm/s^2], 'm'=[m/s^2], 'g'=[g]
LOGT = 'LN';
UNIT = 'g';
% function definition strings ready for use as Matlab functions
sf_x = {
    'c(1)+c(2)*M+c(3)*log(R+c(4)*exp(c(5)*M))+c(6)*H+c(7)*S';
    'c(1)+c(2)*M+c(3)*log(R+c(4)*exp(c(5)*M))+c(6)*H+c(7)*S';
    };
% original paper coeff values:
PCS = coeff_LL08(); % Periods, Coefficients, Sigma
% chosen periods:
per = [0 0.2 0.3 1];
%per = 'all';
%{
per = [ % removed 100 Hz
    0.00 0.02 0.03 0.04 0.05 0.06 0.09 0.10 0.12 0.15 0.17 0.20 ...
    0.24 0.30 0.36 0.40 0.46 0.50 0.60 0.75 0.85 1.00 1.50 2.00 ...
    3.00 4.00 5.00];
%}
% global optimization starting values for each parameter incl. sigmas
X0 = [
    -3.757 1.293 -2.138 2.092 0.237 0.019 0.88 0.05 0.378];
     % -4,  1.50, -2.3, 0.515, 0.633,     0,  0.35, 0.05, 0.2 ];
% optimization/prior lower boundaries (keep as tight as possible)
lowB = [
     -25,    -1,   -5,     0,     0,    -1,    -1,    0,    0];
% optimization/prior upper boundaries (keep as tight as possible)
uppB = [
       2,     5,    1,    12,     2,     1,     2,  2.2,  0.9];
% Prior density function types for each parameter (u=unif,n=norm)
PriorType = {
    'uuuuuuuuu'
    'uuuNNuuuu'
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