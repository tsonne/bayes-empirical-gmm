% Demonstration script for Chapter 4.2.
%
% Load dataset,
% plot dataset,
% set Bayesian inference parameters,
% run Bayesian MCMC,
% evaluate results,
% plot results
%
% 2020-05-12 Tim: created
% 2020-12-08 Tim: refactoring and more documentation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
clear all
close all
clc
addpath('subfunctions'); % add subfunctions to Matlab path
%
% Load dataset (could append multiple different datasets here)
SD = load_42_input();
[Nd,Np] = size(SD.SA);
%
% Name of model configuration function (Matlab m-file without extension),
% the config function can define one simulation model and multiple prior
% settings, resulting in multiple statistical models, e.g. first all priors
% are uniform, the second model has normal priors for coefficient 2, etc.
% It also sets which periods are to be inferred ("per").
mod_conf = 'config_Y2';
%
% Name of MCMC run config function (Matlab m-file without extension),
% you can choose to run the inference or use a saved state from the past,
% further important settings for MCMC and result figures can be set there.
mcmc_conf = 'config_M1';
%
% Combine model and mcmc configs with converted input data:
[SMOD,SDAT,SSET] = setup_BayesInf(mod_conf,mcmc_conf,SD);
%
%
% run function for each GMM for each Dataset:
Nm = numel(SMOD);
Nds = numel(SDAT);
for a=1:Nm % for each prior setting
    for b=1:Nds % for each input dataset
        time0 = tic;
        fprintf('############\n######\n###\n');
        fprintf('##  Starting job: Model %d, Dataset %d\n',a,b);
        GMM_BayesInf(SMOD(a),SDAT(b),SSET);
        fprintf('##  Finished job: Model %d, Dataset %d\n',a,b);
        runtime0 = toc(time0);
        fprintf('##  Job Runtime: %02.0f:%02.0f:%05.2f (%.2f sec)',...
            s2hms(runtime0),runtime0);
        fprintf('###\n######\n############\n');
    end
end
