function SSET = config_M1()
% Bayesian MCMC and output figure/table configuration file.
%
% To run the same model again only to change output figures and tables,
% set RUN_MODE to [0 0] to use the saved inference data.
%
% 2020-05-12 Tim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SETTINGS FOR MCMC:
SSET = struct(...
    'version','1',... % inference version, can give it new names
    'MCMC_ALG','DRAM',... % MCMC algorithm type
    'RUN_MODE',[1 1],... % run the inversion [1 1] or use saved
    ... % data from past run [0 0], first number is for global
    ... % optimization, second number is MCMC inference.
    'F_NL'   , 400,... % NL=round(F_NL*NPAR^1.5) chain length after burn-in
    'NC'     ,  12,... % number of chains to compute
    'NBratio', 0.4,... % Burn-in ratio: NB = round(NBratio*NL)
    'WRIT_MODE',1,... % write mode from local search
    'WRIT_HESS',0,... % write Hessian matrix (on local search mode)
    'WRIT_COVM',0,... % write proposal cov matrix (based on Hessian)
    'PLOT_TRAC',1,... % plot Markov chains (trace plots)
    'PLOT_AUTC',1,... % plot autocorrelation next to chains
    'REMV_OUTL',0,... % remove extreme chain values outside 4 sigma
    'PLOT_HIST',1,... % plot histograms
    'PLOT_ENVL',0,... % plot max.likelihood envelopes over each par.
    'WRIT_INF1',1,... % write MCMC info line (number samples etc...)
    'WRIT_TAB1',1,... % write table: posterior mode & 95% interval values
    'PLOT_ATTC',1,... % plot data and curves of inferred best model
    'PLOT_COR1',1,... % plot large correlation figure (2D hists)
    'PLOT_COR2',1); % plot correlation matrix AND bias (double figure)
