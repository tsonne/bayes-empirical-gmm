function [SMOD,SDAT,SSET] = setup_BayesInf(mod_conf,mcmc_conf,S_rawdata)
% PURPOSE: wrap functions to get config settings for a model and MCMC run
% and convert given dataset to required units and structure.
%
% 2020-05-12 Tim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Model and coefficient configuration:
Sconf = eval(mod_conf);
%
% Setup ready-to-process input from model and coefficient config:
SMOD = extract_coeff(Sconf);
%
% INPUT DATASETS:
% for each dataset, get required period data
Ndsets = numel(S_rawdata);
SDAT = repmat(struct('name','','logY',[],'IMRSH',[],'per',[]),Ndsets,1);
for b=1:Ndsets
    X = S_rawdata(b);
    [A,id]=ismember(Sconf.per, X.per);
    assert(all(A),'Some periods not available in input data!!!')
    SDAT(b).name  = X.name; % dataset name
    SDAT(b).per   = X.per(id); % selected periods
    SDAT(b).IMRSH = [X.EvID, X.M, X.R, X.S, X.H]; % explanatory var.
    SDAT(b).logY  = convert_GMM_units(...  % motion data
        X.SA(:,id), X.LOGT, X.UNIT, Sconf.LOGT, Sconf.UNIT);
end
%
% SETTINGS FOR MCMC:
SSET = eval(mcmc_conf);