function Y = convert_GMM_units(Y,fromLOGT,fromUNIT,toLOGT,toUNIT)
% PURPOSE: CONVERT INPUT GMM DATA BASED ON UNIT INFO.
% INPUT:
%  Y = (matrix) acceleration values in fromLOGT, fromUNIT
%  fromLOGT = input logarithm-type: 'LD'=log10(), 'LN'=log(), ''=no log
%  fromUNIT = input motion unit: 'c'=[cm/s^2], 'm'=[m/s^2], 'g'=[g]
%  toLOGT = output log type
%  toUNIT = output motion unit
% OUTPUT
%  Y = (matrix) acceleration values in toLOGT, toUNIT
%
% 2020-05-12 Tim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% calc antilog using correct base
if strcmpi(fromLOGT,'ld')
    Y = 10.^(Y);
elseif strcmpi(fromLOGT,'ln')
    Y = exp(Y);
end
% convert to [m/s^2] if not already
if strcmpi(fromUNIT,'c') % from [cm/s^2]
    Y = Y/100;
elseif strcmpi(fromUNIT,'g') % from [g]
    Y = Y*9.81;
end
% convert to final unit
if strcmpi(toUNIT,'c') % to [cm/s^2]
    Y = Y*100;
elseif strcmpi(toUNIT,'g') % to [g]
    Y = Y/9.81;
end
% convert to final log
if strcmpi(toLOGT,'ld')
    Y = log10(Y);
elseif strcmpi(toLOGT,'ln')
    Y = log(Y);
end