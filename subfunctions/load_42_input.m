function SD = load_42_input()
% Read input data files
%
% OUTPUT
%   SD = struct of data
%
% 2020-05-12 Tim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input file paths
input_dir = fullfile('.','input_data');
f_expl = fullfile(input_dir, 'Explanatory.csv');
f_pgsa = fullfile(input_dir, 'Data_PSA_values.csv');
%
% Read files
fid = fopen(f_expl);
X = textscan(fid,'%s%f%f%f%f%f','delimiter',',');
fclose(fid);
Y = load(f_pgsa,'-ascii');
%
% Deal data to separate variables
SD = struct();
SD.name = 'is83ria';
SD.LOGT = ''; % data are scaled by log ('LD'=log10, 'LN'=log, ''=no log)
SD.UNIT = 'c'; % data in motion unit ('c'=[cm/s^2], 'm'=[m/s^2], 'g'=[g])
SD.SID = X{1,1}; % station ID codes
SD.EvID = X{1,2}; % Event IDs
SD.M = X{1,3}; % moment magnitude
SD.R = X{1,4}; % source-receiver distance (Joyner-Boore type) [km]
SD.S = X{1,5}; % site type (0=rock, 1=stiff soil)
SD.H = X{1,6}; % hypocentral depth [km]
SD.per = Y(1,:); % pseudospectral oscillation period [s]
SD.SA = Y(2:end,:); % spectral acceleration [cm/s^2]
