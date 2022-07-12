% Bayesian inference on GMPE parameter marginal posterior densities,
%   using global and local optimization for mode finding, then MCMC.
%   Can define different GMM functions (GMPEs) and prior assumptions.
%   Different MCMC algorithms can be used:
%   + Roberts et al. (1997) based multivar.normal proposal in MCMC.
%     (problematic when very close to value boundaries,
%      inappropriate in case of multimodal densities...)
%   + Delayed Rejection Adaptive Metropolis (DRAM) [through full chain]
%     (Ref: Haario, Laine, Mira & Saksman 2006 Stat Comput)
%     (most reliable, most efficient, results in densest/smoothest PDFs)
%     (using it after Burn-in is non-Markovian, needs justification...)
%   + Burn-in DRAM (BiDRAM): DRAM only for Burn-in of chain, afterwards
%     uses cov.-based multivar.normal proposal (as Roberts et al. (1997))
%     (reliable, much less efficient sampling after burn-in, but Markovian)
%   + Multi-Stage Adaptive Metropolis (MSAM)
%     experimental 'robust' adaptation scheme through multiple stages with
%     simple proposal tuning during Burn-in. Difficulties with non-singular
%     covariance matrices make this rather useless for now. 
%     Attempt to combine multiple burn-in stages with tuning, AM and
%     RobertsEtAl(1997). Needs >1 chains (for median after each stage).
%   MATLAB TOOLBOXES REQUIRED:
%   The MCMC algorithms use Matlab's parfor loop, ParallelComputing Toolbox
%   needed. Statistics Toolbox needed as well for some functions.
%
%
%  2016-06-08 Tim (Version 1)
%  2016-06-30: tsonne:
%    Choose different posterior functions (and priors) with strings.
%    Remove Unused MCMC algorithms (for now) to clean up script.
%  2016-07-01: tsonne:
%    Clear separation of functions: Prior, Likelihood, Posterior, Simul..
%    Use global search for mode estimate to make manual start value less
%    significant. Still have to define parameter boundaries.
%    Now with loop over available periods. <-- getting problem at one
%      period: Warning: Matrix is singular to working precision. (Hessian)!
%      Reason: mode(tau) is assumed 3e-06, numHessian h=0.0001 ...
%      Quickfix by setting lower limit for sigma parameters: 0.002
%  2016-07-04: tsonne:
%    Quickfix (mode) can cause Hessian diagonal to be positive, then cov
%    has negative value on diag.--> not pos. def. matrix anymore!
%    Fix: find neg. var and change sign OR apply waterlevel on diagonal:
%    works!
%    Save mode values and statistics of all periods as tables and mat.
%  2016-07-05: tsonne:
%    Made all figures and output tests optional (PLOT_TRAC etc.).
%    New Output: 2 ASCII files (input data table, sim data table with ETAi)
%  2016-07-07: tsonne (VERSION 6)
%    New construction of probability functions: Have one function for
%    Prior, Sim., Cov.matrix and log-multivar.normal pdf ==> Posterior.
%  2016-11-11: tsonne (VERSION 7)
%    Added DRAM MCMC option. 
%  2016-11-30: tsonne (VERSION 8)
%  2016-12-05: tsonne:
%    +migrated folder structure to 'project' directory.
%    +implemented LOAD_SETTINGS option, can store and load specific GMPE-
%      settings using premade functions like 'get_settings_Am2005'.
%    +output coefficient MLE ascii file in /dat/V8_REF_coeff.txt.
%    +zip/unzip all SET-specific files now, keep output file number small. 
%  2017-02-06 tsonne (VERSION 9)
%  2017-02-11 tsonne:
%    + strong motion data source now EQS file database (optional)
%    + removed EXCEL reading capability, but can use *.mat files instead.
%      --> much less problematic than Excel files on all systems/versions.
%  2017-03-07 tsonne: (VERSION 10)
%    + added Multi-Stage Adaptive Metropolis (MSAM) routine.
%    + added Burn-in DRAM (BiDRAM), then continue without adaptations.
%    + optimized loglikelihood-fct for random-(event)-effect model, now
%      MCMC is about twice as fast (using 65 records --> COV: [65 x 65]).
%    + improved automatic handling of parallel processing: number of cores
%      in matlabpool will be limited according to system info, so usage on
%      different machines should always be optimal and without errors.
%    + now collecting and zipping up all figure files in one ZIP-archive
%      file (speeds up sync, less clutter).
%  2017-03-13 tsonne: (VERSION 10MK)
%    + some restructuring to allow standalone version with code subdir.
%    + supposed to run on any OS now, in self-contained folder.
%  2017-03-27 tsonne: (VERSION 11MK)
%    ???
%  2019-02-24 tsonne: VERSION 12
%    + fixed two plots, Matlab R2015a is too different! WHY?!
%      + histograms (labels fixed, added more prencetile lines)
%      + correlation matrix plot (labels fixed)
%    + parallel core determination fails on hyperthead cores,
%      (get half of threads), removed auto-determination, fixed value to 6
%  2019-02-26 tsonne: VERSION 13
%    + different Bayesian prior for each variable,
%      give char string for prior type and the usual 2 vectors for
%      2 prior parameters (unif: low & upp, norm: mean & std)
%  2019-03-12 tsonne: VERSION 14 = GMM_BayesInf VERSION 3
%    + can give X0, lowB, uppB as matrices, one row per period.
%    + if only row vector given, same values for all periods (like before)
%  2019-03-20 tsonne: GMM_BayesInf VERSION 4
%    + extra summary figures...
%  2020-12-09 tsonne: lots of minor changes, refactoring for standalone.
%
function GMM_BayesInf(SMOD,SDAT,SSET)
%==========================================================================
%% SETUP
FS   = filesep; % OS-dep. file seperator
%=========================================================================%
%% OUTPUT DATA LOCATIONS
m_name = SMOD.name; % model name
d_name = SDAT.name; % dataset name
VERSION = SSET.version; % inference script/setting version
%
DOCREF = [m_name '-' d_name '-V' VERSION]; % project reference string
DMAT = ['dat' FS]; % local data dir (for output files)
DMAT2 = [DMAT DOCREF FS ]; % specific output folder
DMATV = [DMAT2 DOCREF '_']; % output file path base
MEMONAME = ['MCMC-' DOCREF]; % file name of LaTeX PDF
Date = date; % LaTeX PDF date under title
Subject = ['MCMC ' DOCREF]; % LaTeX PDF subject/title line
%=========================================================================%
%% MODEL SETTINGS
%
% GMM sim function string, e.g. 'c(1)*M+c(2)*log(R)',
% GMM_BayesInf() assigns variables M,R,S,H from data to be used in
% evaluated simulation function string sf.SIM
sf.SIM = SMOD.sf;
RND_EFFECT = SMOD.RND_EFFECT; % use random effects model?
% cell of 'u'|'n', OR char:'uniform'|'normal':
PRIOR_TYPE = lower(SMOD.Prior);
% Periods to loop through:
Per = SDAT.per; % periods [s] if fixed 'SET_PERIODS = 1'
%=========================================================================%
%% MODE FINDING AND MCMC SETTINGS
%   RUN_MODE(1) = 1 : do Mode & Hessian calculation,
%   RUN_MODE(2) = 1 : do MCMC calculation,
%   Otherwise load from previous save files (if avail.).
RUN_MODE = SSET.RUN_MODE;
Nreanneal =  200; % number of reannealing intervals to run (global search)
fmNmaxiter= 2000; % max. number of fminsearch iterations allowed
MCMC_ALGR =SSET.MCMC_ALG;% 'Roberts1997','DRAM','BiDRAM' or 'MSAM'
NSTAGES   =    1; % number of adaptive stages, last stage produces finals
% factor F_NL of chain length after burn-in: NL = round(F_NL * NPAR^1.5)
F_NL = getpar(SSET,'F_NL',400);
MaxWorkers=    6; % maximum parallel workers to allow,
                  % if MaxWorkers>NC, will be set to number of chains.
% number of chains to compute (parallel computing):
NC = getpar(SSET,'NC',12);
% Burn-in ratio: NB = round(NBratio*NL):
NBratio = getpar(SSET,'NBratio',0.4);
% after first set, use cov matrix of previous set to guide next set?
USE_LAST_COV = getpar(SSET,'USE_LAST_COV',true);
% numerical stability factor in DRAM MCMC
pfact = getpar(SSET,'pfact',[]);
%=========================================================================%
%% OUTPUT SETTINGS
CORR_TYPE = 'Pearson'; % 'Pearson' or 'Spearman' [if relationship is 
%                        linear, Pearson is fine. Spearman takes more time]
 % write mode from local search
WRIT_MODE = getpar(SSET,'WRIT_MODE',1);
% write Hessian matrix (on local search mode)
WRIT_HESS = getpar(SSET,'WRIT_HESS',0);
% write proposal cov matrix (based on Hessian)
WRIT_COVM = getpar(SSET,'WRIT_COVM',0); 
% plot Markov chains (trace plots)
PLOT_TRAC = getpar(SSET,'PLOT_TRAC',0); 
% plot autocorrelation next to chains
PLOT_AUTC = getpar(SSET,'PLOT_AUTC',1);
% remove extreme chain values outside 4 sigma
REMV_OUTL = getpar(SSET,'REMV_OUTL',0); 
% plot histograms
PLOT_HIST = getpar(SSET,'PLOT_HIST',1);
% plot max.likelihood envelopes over each par.
PLOT_ENVL = getpar(SSET,'PLOT_ENVL',0);  
% write MCMC info line (number samples etc...)
WRIT_INF1 = getpar(SSET,'WRIT_INF1',1);
% write table: posterior mode & 95% interval values
WRIT_TAB1 = getpar(SSET,'WRIT_TAB1',1);
% plot data and curves of inferred best model
PLOT_ATTC = getpar(SSET,'PLOT_ATTC',0);
% plot large correlation figure (2D hists)
PLOT_COR1 = getpar(SSET,'PLOT_COR1',0); 
% plot correlation matrix AND bias (double figure)
PLOT_COR2 = getpar(SSET,'PLOT_COR2',1);
%=========================================================================%
%% GRAPHICS SETTINGS
FIGTYPE = 'png'; % output file type of figures
VIS   = 'off';
AXIS_LIM = 'tight'; % if 'tight' axis is min,max(data), else use lb,ub
UN  = {'Units','Normalized'};
HAlr  = {'HorizontalAlignment','right'};
pLW   = {'LineWidth',1.5}; % linewidth
gray  = .5*[1 1 1];
LIMr = [0.1  200]; % plot axis limits: distance
LIMm = [4.5    7]; % plot axis limits: magnitude
LIMb = [-1.1 1.1]; % plot axis limits: bias
LIMe = [-0.3 0.3]; % plot axis limits: event term
% correlation plot colormap: [2D density plots]
c = [255, 0  , 255        % violet
     25 , 25 , 225        % blue [mod]
     35 , 225, 225        % cyan [mod]
     0  , 255, 0  ]./255; % green
x = [0  0.333  0.667  1]; % rel. position
CMap = make_cmap(c,x,256); % make this colormap happen
CMap = cmap_whiten(CMap,2); % fade colormap to white for low values
% correlation plot colormap: [correlation matrix]
c = [0  , 0  , 255        % blue
     255, 255, 255        % white
     255, 0  , 0  ]./255; % red
x = [0  0.5  1]; % rel. position
CMap2 = make_cmap(c,x,256); % make this colormap happen
% optional settings
option_pdflatex = '-q';
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ASSIGN INPUT DATA AND EXPLANATORY VARIABLES
%
LogY = SDAT.logY; % expect some log(amp), e.g. log10([m/s^2])
EvID = SDAT.IMRSH(:,1); % earthquake ID in table
M = SDAT.IMRSH(:,2); % Mw list
R = SDAT.IMRSH(:,3); % Rjb list
S = SDAT.IMRSH(:,4); % Site number list (typically 0 or 1)
H = SDAT.IMRSH(:,5); % hypocentral depth list
Nm  = numel(unique(EvID)); % number of events
Ny  = size(SDAT.IMRSH,1); % number of all data points (each period)
vNS = nan(Nm,1); % number of records for each event
Mw = nan(Nm,1); % unique Mw unsorted (as input)
for a=1:Nm
    vNS(a) = sum(EvID==a);
    Mw(a) = unique(M(EvID==a));
end
MRS = [M R S]; % set matrix as required for GMM function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ASSIGN LIKELIHOOD FUNCTION
%
% The data covariance matrix C is input to the likelihood function.
% If all data are treated as independent, C is a unity matrix multiplied
% by sigma^2 (variance). However, if random effects model used, i.e.
% data of same event assumed to be dependent (event effect), then C is
% a block-diagonal matrix as produced by cov_1effect(vNS,x) with two
% variance parameters (tau and sigma) and structure based on number of
% data per event (vNS).
%
% Covariance matrix definition:
switch lower(RND_EFFECT)
    case 'none'
        sf.COV = 'eye(Ny)*x^2'; % simple
        NSIG = 1;
    case 'event'
        %sf.COV = 'cov_1effect(vNS,x)'; % not optimized function
        sf.COV = '[optimized_in_llikf]'; % implemented in ln_re_mvnpdf()
        NSIG = 2;
    otherwise
        error('Unknown RND_EFFECT option!');
end
% Parameter Indices:
NPAR = numel(SMOD.X0(1,:)); % number of parameters
iX = 1:NPAR-NSIG; % model parameters
iZ = NPAR-NSIG+1:NPAR; % statistical deviance parameters (sigma_x)
NiX = numel(iX); NiZ = numel(iZ);
%
% Likelihood function definition:
switch lower(RND_EFFECT)
    case 'none'
        sf.LIK = 'lnmvnpdf(y,fh.SIM(x(iX)),fh.COV(x(iZ)))';
    case 'event'
        %sf.LIK = 'lnmvnpdf(y,fh.SIM(x(iX)),fh.COV(x(iZ)))';
        sf.LIK = 'ln_re_mvnpdf(y,fh.SIM(x(iX)),vNS,x(iZ))';
    otherwise
        error('Unknown RND_EFFECT option!');
end
%
% MCMC settings:
NL = round(F_NL * NPAR^1.5); % sample size L to be drawn (after burn-in)
NB = round(NBratio*NL); % burn-in samples at beginning
NT = NL + NB; % total number of iterations
MPoolN = min(NC,MaxWorkers); % as many parallel as chains, or less
if NC<MaxWorkers
    fprintf('Using %d workers due to number of requested chains: %d\n',...
        MPoolN,NC);
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ASSIGN LABEL STRINGS
Par_labels = cell(1,NPAR); % Matlab plots
Par_names  = cell(1,NPAR); % ASCII
header_out = cell(1,NPAR); % LaTeX (Tables etc.)
format_out = cell(1,NPAR); % LaTeX tables number format
for a = 1:NiX % for all model parameters
    Par_labels{a} = ['c_{' num2str(a) '}'];
    Par_names{a}  = ['c_{' num2str(a) '}'];
    header_out{a} = ['$c_{' num2str(a) '}$'];
    format_out{a} = '%6.3f';
end
for a = 1:NiZ % for all covariance parameters
    Par_labels{NiX+a} = ['\sigma_' num2str(a)];
    Par_names{NiX+a}  = ['sigma_' num2str(a)];
    header_out{NiX+a} = ['$\sigma_' num2str(a) '$'];
    format_out{NiX+a} = '%6.3f';
end
sum_tab_fmt = cell(1,NPAR);
sum_tab_fmt(:) = {'%8.5f'};
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SET PERIODS
% (SetA, SetB from older script version, always fixed here)
NPER = numel(Per); % number of processed periods (output)
SetA = 1; % first period to use
SetB = NPER; % last one
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SET OUTPUT FILE PATHS
% make data dir:
[~,~,~] = mkdir(DMAT);
[~,~,~] = mkdir(DMAT2);
% cell for all par.values at all periods:
% {per,mode,2.5%,50%,97.5%,mean,stdev,corr,etai,Mw}
OUTstats = cell(NPER,10);
sOUTstats= [DMATV 'PostStats.mat']; % output file name
%
sOUTtab3 = [DMATV 'coeff.txt']; % MLE coeff. output table file name
sOUTtab4 = ... % MLE coeff. output table file name
    [DMATV 'rec_coe_s' num2str(NiZ) '.mat'];
sARCHIVE = [DMATV 'ARCHIVE.zip']; % zip archive name for all SET files
% unzip archive if specified:
if any(RUN_MODE==0)
    assert(exist(sARCHIVE,'file')>0,'Cannot find ARCHIVE file!');
    unzip(sARCHIVE); % unzips all the SET mode,hess,info,px files in ./dat/
end
% keep all output data file names to zip them into archive at end:
cSAVED = cell(NPER,4); % {mode, hess, info, px} for each period
%
% Input data table file
xOUT1 = [M S H R LogY(:,SetA:SetB)];
sOUTmat1 = [DMATV 'input_table.mat'];
if RUN_MODE(2) == 1
    save(sOUTmat1,'xOUT1');
end
%
% save input settings/configs and values
if any(RUN_MODE==1)
    s_InStructs = [DMATV 'SModDatSet.mat']; % used input values/settings
    save(s_InStructs,'SMOD','SDAT','SSET');
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SET FIGURE ARCHIVE AND LATEX FILE
% LaTeX output file for tables, figures:
[~,p_fig,~,p_aux,fbase,fhead,fbody,fpdf,fpdf2] = ...
        ts_prep_tex(MEMONAME,Date,Subject);
fid = fopen(fbody,'w');
%
% ZIP-Archive to save all figures:
sFIGARCH = [p_fig FS fbase '.zip'];
Nfig = (3*NPAR+5)*NPER+2; % (max) number of figure-files produced
cFIG = cell(Nfig,1); % list of all figure file names
cntF = 0; % figure counter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% START LOOP OVER PERIODS
%
for SET  = SetA:SetB; % for each chosen period...
    %
    % START OF LOOP FOR EACH PERIOD
    %
cnt1 = SET - SetA + 1; % loop counter from 1 corresponding to SetA
sSET = sprintf('_%d',SET); % current output filename suffix
xPer = Per(SET); % current period or frequency
yy = LogY(:,SET); % data of current period
%
%
% PRIOR TYPES:
% Assign log(prior) function,
% set lower/upper bounds for global search lowBx and uppBx, which also
% determines simulated annealing initial temperature for each parameter.
if size(SMOD.X0,1)>1 % check if each period has its own values
    X0 = SMOD.X0(SET,:);
else
    X0 = SMOD.X0;
end
if size(SMOD.lowB,1)>1 % check if each period has its own values
    lowB = SMOD.lowB(SET,:);
else
    lowB = SMOD.lowB;
end
if size(SMOD.uppB,1)>1 % check if each period has its own values
    uppB = SMOD.uppB(SET,:);
else
    uppB = SMOD.uppB;
end
PriC1 = lowB; % prior pdf coefficients,
PriC2 = uppB; % go into function string to be evaluated
%
if any(strcmp(PRIOR_TYPE,'uniform')) || all(PRIOR_TYPE=='u')
    % if all parameters have uniform priors:
    PRIOR_TYPE = repmat('u',1,NPAR);
    sf.PRI = 'lnPriorUnif(x,PriC1,PriC2)';
    lowBx = lowB; % global optimization lower boundaries
    uppBx = uppB; % global optimization upper boundaries
elseif any(strcmp(PRIOR_TYPE,'normal')) || all(PRIOR_TYPE=='n')
    % if all parameters have normal priors:
    PRIOR_TYPE = repmat('n',1,NPAR);
    sf.PRI = 'lnPriorNorm(x,PriC1,PriC2)';
    lowBx = X0 - 4*uppB; % +- 4 sigma
    uppBx = X0 + 4*uppB;
else
    % parameters have both uniform and normal priors,
    % assume char string of 'u' and 'n' mixed
    assert(all(PRIOR_TYPE=='n'|PRIOR_TYPE=='u'),...
        ['PRIOR_TYPE must be "uniform", "normal" or combination of ',...
        '"u" and "n": %s'],PRIOR_TYPE);
    assert(numel(PRIOR_TYPE)==NPAR,...
        ['PRIOR_TYPE "un"-string must have same number as ',...
        'parameters: %s = %d, NPAR = %d'],...
        PRIOR_TYPE,numel(PRIOR_TYPE),NPAR);
    sf.PRI = 'lnPriorNU(x,PRIOR_TYPE,PriC1,PriC2)';
    lowBx = lowB;
    uppBx = uppB;
    iPN = PRIOR_TYPE=='n';
    lowBx(iPN) = lowB(iPN)-4*uppB(iPN); % +- 4 sigma
    uppBx(iPN) = lowB(iPN)+4*uppB(iPN);
end
%
% Function strings and handles:
fh.PRI = eval(['@(x)   ' sf.PRI]); % log(prior) handle
fh.SIM = eval(['@(c)   ' sf.SIM]); % simulation handle
fh.COV = eval(['@(x)   ' sf.COV]); % cov.matrix fct handle
fh.LIK = eval(['@(x,y) ' sf.LIK]); % log(likelihood) handle
sf.POS = sprintf('%s + %s',sf.LIK,sf.PRI); % log(posterior) string
fh.POS = eval(['@(x,y) ' sf.POS]); % log(posterior) handle
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% WRITE SOME INFO TO LATEX FILE ON SET 1
if SET==SetA % first set
    put(fid,'\chapter{INFO}','');
    % in case this is loaded from existing calculations
    if RUN_MODE(2) == 0
        s_INFO = [DMATV 'INFO_SET_1.mat'];
        assert(exist(s_INFO,'file')>0,'Cannot find saved INFO struct.');
        load(s_INFO,'INFO');
        sf = INFO.sf;
    end
    %
    % Table: Used statistical functions
    header_1 = 'Purpose & Function handle';
    cBOD1    = { 'Simulation',sf.SIM; 'Covariance',sf.COV;...
        'Prior',sf.PRI; 'Likelihood',sf.LIK; 'Posterior',sf.POS;};
    cap = 'Used functions in this script.';
    tsLatexTable(fid,header_1, cBOD1, [], 'll', cap,[],'literal'); % TAB
    %
    % Figure: Mw vs. distance
    put(fid,'\chapter{Magnitude distance distribution}');
    s = fullfile(p_fig,'MR.png');
    plot_mag_dist(M,R,S,'off',s);
    % SAVE IMAGE AND UPDATE TEX FILE:
    cap = ['Magnitude-distance distribution of dataset used in this ',...
        'study. Site types are indicated for rock (black diamonds) ',...
        'and soil (grey squares).'];
    ts_WriteLaTeXstringsNxM(fid,s,[1 1 1 1],cap,0.5);
    close
    cntF = cntF+1;
    cFIG{cntF} = s;
    %
    put(fid,'','\newpage','');
end
%
fprintf(fid,'\\chapter{Period %4.2f s}\n',xPer);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODE and HESSIAN
%
% simulannealbnd() requires 'Global Optimization Toolbox'
% fminsearch() might require 'Optimization Toolbox'
%
s_mode    = [DMATV 'PostMode' sSET '.mat'];
s_hess    = [DMATV 'PostHess' sSET '.mat'];
cSAVED{cnt1,1} = s_mode;
cSAVED{cnt1,2} = s_hess;
if RUN_MODE(1) == 1
    % MODE:
    if SET==SetA
        % Global Search: Simulated Annealing
        ReannItv = 100;
        Nmaxiter = ReannItv*Nreanneal;
        MaxFeval = round(Nmaxiter*1.1);
        TolFun   = 1e-8;
        StallIL  = 50000;
        Tinit    = 2*(uppBx-lowBx);
        saoptions= saoptimset('StallIterLimit',StallIL,...
            'MaxIter',Nmaxiter, 'MaxFunEvals',MaxFeval, 'TolFun',TolFun,...
            'ReannealInterval',ReannItv, 'InitialTemperature',Tinit...
            );
        X1 = simulannealbnd(@(x) -fh.POS(x,yy),X0,lowBx,uppBx,saoptions);
        pause(0.1);
        close
        if fh.POS(X0,yy) > fh.POS(X1,yy) % simulannealbnd isn't always good
            X1 = X0;
        end
    else
        X1 = MODE;
    end
    % Local Search: Simplex
    fmoptions = optimset('MaxIter',fmNmaxiter);
    mode_th = fminsearch(@(x) -fh.POS(x,yy), X1, fmoptions);
    % Force sigma terms to be significantly nonzero:
    mode_th(iZ) = max(mode_th(iZ), 0.002);
    % HESSIAN:
    h = ones(1,numel(X0))*0.0001; % step size numerical Hessian
    He = my_hessian(fh.POS, mode_th, h, yy);
    % save it:
    save(s_mode,'mode_th');
    save(s_hess,'He');
    pause(0.1);
    close
    if SET~=SetA
        mode_th = LAST_MEDIAN;
    end
else
    load(s_mode,'mode_th');
    Hx = load(s_hess);
    if isfield(Hx,'H');
        He = Hx.H;
    else
        He = Hx.He;
    end
end
% COVARIANCE:
scale = 2.38^2/NPAR; % scaling for cov Gelman etal 1996, Roberts etal 1997
Pcov = -scale*(He\eye(NPAR)); % covariance of proposal pdf
% waterlevel: set diagonal cov to min 0.0001
dCz = diag(Pcov)<=0;
if any(dCz)
    ii = find(dCz);
    Pcov(ii,ii) = max(0.001,-Pcov(ii,ii));
end
% WRITE TO TEX FILE:
if WRIT_MODE
put(fid,'Mode:');
WriteLaTeXmatrix(fid,'\pmb{\hat{\theta}}',mode_th,'%5.3f');
end
if WRIT_HESS
put(fid,'Hessian:');
WriteLaTeXmatrix(fid,'\pmb{H}(\pmb{\hat{\theta}})',He,'%.1e');
end
if WRIT_COVM
put(fid,'Proposal covariance:');
WriteLaTeXmatrix(fid,'\pmb{C}(\pmb{\hat{\theta}})',Pcov,'%.1e');
end
% Check before doing potentially stupid things:
%assert(~any(isnan(Pcov(:))),'Covariance is NaN! Abort.');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MCMC
%
s_INFO = [DMATV 'INFO_SET' sSET '.mat']; % save file names
s_PX   = [DMATV 'PX_SET'   sSET '.mat']; % save file names
cSAVED{cnt1,3} = s_INFO;
cSAVED{cnt1,4} = s_PX;
%
% sampling parameters have been set above in script...
%
% NSTAGES set at top of script
NBs = round(NB/NSTAGES);
NB0 = round(NBs*0.1);
In0 = nan(NSTAGES-1,1);
In1 = In0; In2 = In0;
assert(NSTAGES>0,'This script requires NSTAGES>0.');
if strcmp(MCMC_ALGR,'MSAM') % burn-in stages + final stage, so N>1
    assert(NSTAGES>1,'Method MSAM requires NSTAGES>1.');
end
if NSTAGES == 1
    NB0 = 0; In0 = []; In1=[]; In2=[];
else
    for stg = 1:NSTAGES-1
        In0(stg) = 1 + (stg-1)*NBs; % first index
        In1(stg) = In0(stg) + NB0; % first after burn-in
        In2(stg) = stg*NBs; % last index
    end
end
%
maxlag = 50; % maximum lag for later autocorrelation analysis
%
% memory init:
FIX   = zeros(1,NPAR);
FREE  = 1:NPAR;
% first proposal, will force values to be within given bounds:
if any(strcmp(PRIOR_TYPE,'uniform')) || all(PRIOR_TYPE=='u')
    lb = lowB;
    ub = uppB;
elseif any(strcmp(PRIOR_TYPE,'normal')) || all(PRIOR_TYPE=='n')
    lb = mode_th - 2*uppB; % limit to a few stdev +- mode
    ub = mode_th + 2*uppB;
else
    lb = lowB;
    ub = uppB;
    iPN = PRIOR_TYPE=='n';
    lb(iPN) = mode_th(iPN) - 2*uppB(iPN); % +- 2 prior sigma
    ub(iPN) = mode_th(iPN) + 2*uppB(iPN);
end
F.sp = @(x,sd) proposal_norm_lim(x,sd,NPAR,lb,ub);
%
ddd = {NC, NL, NB, NT, NPAR, FREE, FIX, X0, lb, ub, Par_names, ...
    Par_labels, maxlag, NSTAGES, NBs, NB0, In0, In1, In2, sf, fh};
if strcmp(MCMC_ALGR,'DRAM')
    [px,INFO,PV,PP] = get_PxInfoStruct(ddd{:});
else
    [px,INFO,PV,PP,PR] = get_PxInfoStruct(ddd{:});
end
%
% Probability functions (handles in struct F):
F.po = @(x) fh.POS(x,yy); % log(Posterior)
TuneItv = 100; % tune interval to assess acceptance rate (also for output)
drscale = 3; % DR proposal variance scaling factor
Pwid  = (ub-lb); % unif.PDF width [proposal]
if RUN_MODE(2) == 1 % create initial proposal covariance
    if SET==SetA || ~USE_LAST_COV || ~strcmpi(MCMC_ALGR,'dram')
        inicov = diag((Pwid*0.1).^2); % initial covariance
    else
        inicov = NEW_INICOV; % previous DRAM run result covariance
    end
end
%
%%
time1 = tic;
if RUN_MODE(2) == 1 % do the MCMC run
    %%%%%%%%%%%%%%%%%
    %%% MCMC run: %%%
    %%%%%%%%%%%%%%%%%
    % (attention: chains are parallel, using parfor)
    if isempty(gcp('nocreate')) %matlabpool('size') == 0
        %matlabpool('open',MPoolN)
        poolobj = parpool(MPoolN);
        addAttachedFiles(poolobj, 'subfunctions/') 
    end
    switch MCMC_ALGR
        case 'Roberts1997' % SIMPLE MCMC WITH SET COV.MATRIX
            SD0  = sqrt(diag(Pcov)')*0.5; % start value deviation
            % NormProposal by Roberts etal 1997:
            F.pp = @(x,r) x + (chol(Pcov)'*r')';
            % MCMC by Roberts et al 1997:
            parfor c = 1:NC
                a = 1;
                acn  = 1; % acceptance counter
                acr  = 0; % acceptance rate
                x_t  = F.sp(mode_th,SD0); % propose start values
                p_t  = F.po(x_t); % logposterior
                PV{c}(a,:) = x_t; % save
                PP{c}(a)   = p_t; % save
                for a = 2:NT
                    r  = PR{c}(a,:); % get saved random numbers
                    x  = F.pp(x_t,r); % propose new values
                    p  = F.po(x);     % logposterior
                    if p-p_t > log(rand) % Metropolis step: accept
                        p_t = p;
                        x_t = x;
                        acn = acn+1;
                    end
                    PV{c}(a,:) = x_t; % save
                    PP{c}(a)   = p_t; % save
                    if a<NB && mod(a,TuneItv)==0 % assess acc.rate and tune
                        acr  = acn/TuneItv; % rate
                        acn  = 0; % reset counter
                        fprintf('W:%d C:%6d AR:%5.3f\n',c,a,acr);
                    end
                end
            end
            %
        case 'DRAM' % DELAYED REJECTION ADAPTIVE METROPOLIS
            model.ssfun    = F.po; % log posterior function
            params.par0    = mode_th; % initial value
            options.nsimu    = NT;
            options.adaptint = TuneItv;
            options.drscale  = drscale;
            options.qcov     = inicov; % initial covariance
            options.pfact    = pfact; % num.stability factor (20190324TS)
            results = cell(NC,1);
            parfor c = 1:NC
                [results{c},PV{c},PP{c}] = tsdramrun(...
                    model,params,options);
            end
            % store best acceptance cov matrix, median has been not
            % positive definite on occasions:
            [~,Ixx]=max(cellfun(@(x)x.accepted,results));
            NEW_INICOV = results{Ixx}.cov;
            %
        case 'BiDRAM' % DELAYED REJECTION ADAPTIVE METROPOLIS BURN-IN
            model.ssfun    = F.po; % log posterior function
            params.par0    = mode_th; % initial value
            options.nsimu    = NB;
            options.adaptint = TuneItv;
            options.drscale  = drscale;
            options.qcov     = inicov; % initial covariance
            PV1   = cell(NC,1); % Burn-in model parameter value cell
            PP1   = cell(NC,1); % Burn-in log-likelihood value cell
            PV1(:)= {nan(NB,NPAR)};
            PP1(:)= {nan(NB,1)};
            results = cell(NC,1);
            F.pp = @(x,r,c) x + r*c;  % final proposal function
            parfor c = 1:NC
                [results{c},PV1{c},PP1{c}] = tsdramrun(...
                    model,params,options);
                PV{c}(1:NB,1:NPAR) = PV1{c};
                PP{c}(1:NB) = PP1{c};
                % Run final stage (no tuning):
                cR = results{c}.R; % R = scale*chol(COV+eps*eye(NPAR));
                x_t  = PV{c}(NB,1:NPAR); % last burn-in value
                p_t  = PP{c}(NB);        % loglik
                acn  = 0; % acceptance counter
                for a = NB+1 : NT
                    r  = PR{c}(a,:);  % get saved random numbers
                    x  = F.pp(x_t,r,cR); % propose new values
                    p  = F.po(x);     % loglik
                    if p-p_t > log(rand) % Metropolis step: accept
                        p_t = p;
                        x_t = x;
                        acn = acn+1;
                    end
                    PV{c}(a,:) = x_t; % save
                    PP{c}(a)   = p_t; % save
                    if mod(a,TuneItv)==0 % just print iter#,acr
                        acr = acn/(a-NB); % acc.rate of entire stage
                        fprintf('W:%d C:%6d AR:%5.3f\n',c,a,acr);
                    end
                end
            end
            %
        case 'MSAM' % MULTI-STAGE ADAPTIVE METROPOLIS
            % adaptation by intra-stage tuning and inter-stage var.analysis
            model.ssfun    = F.po; % log posterior function
            model.stfun    = F.sp; % start value generating function
            params.par0    = mode_th; % initial value
            options.adaptint = TuneItv;
            options.qcov     = inicov; % initial covariance
            options.SD0     = SD0; % initial stdev vector
            [results,PV,PP] = ts_MSAMv2(...
                model,params,options,INFO,PV,PP,PR);
            % This performs rather badly, the adaptation is not good enough
        otherwise
            error('Available MCMC methods: Roberts1997, DRAM, MSAM.');
    end % choice of MCMC Method
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Analysis of Results
    % gather chains to struct, get various statistical point values:
    px = analyse_MCMC_results(px,NPAR,NC,NB,NL,PV,PP,maxlag);
    INFO.fh = [];
    save(s_INFO,'INFO');
    save(s_PX,'px');
else % DO NOT SIMULATE, load saved data instead:
    assert(exist(s_INFO,'file')>0,'Cannot find saved INFO struct.');
    assert(exist(s_PX,'file')>0,'Cannot find saved data struct PX.');
    load(s_INFO,'INFO');
    load(s_PX,'px');
    NC = INFO.NC; NL = INFO.NL; NB = INFO.NB; NT = INFO.NT;
    NPAR = INFO.NPAR; FREE = INFO.FREE; FIX = INFO.FIX;
    X0 = INFO.X0; lb = INFO.lb; ub = INFO.ub;
    Par_names = INFO.Par_names; Par_labels = INFO.Par_labels;
    NSTAGES = INFO.NSTAGES; NBs = INFO.NBs; NB0 = INFO.NB0;
    In0 = INFO.In0; In1 = INFO.In1; In2 = INFO.In2;
    sf = INFO.sf; fh = INFO.fh;
    if isempty(fh)
        fh.SIM = eval(['@(c)   ' sf.SIM]); % prediction handle
    end
end
tEl1 = toc(time1);
fprintf('MCMC set %d: %10.3f sec, %8.6f s/exec\n',SET,tEl1,tEl1/NT/NC);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODE DECISION
% Which is better: local search or MCMC?
if F.po(mode_th) > F.po([px(:).Tmode])
    MODE = mode_th;
else
    MODE = [px(:).Tmode];
end
tperc = [px(:).Tperc];
LAST_MEDIAN = tperc(3,:); % 50.0 percentile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT CHAINS and HISTOGRAMS
%
if PLOT_TRAC % if plot wanted
% get plot value limits and reasonably spaced ticks:
[pl,ticks] = get_lims_ticks(px,lb,ub,FREE,AXIS_LIM,'yes',NB);
% Traceplots:
for a = 1:NPAR
    ts_fig(8,1.2,VIS);
    if PLOT_AUTC % chain position
        subplot('Position',[0.09 0.19 0.70 0.75])
    else
        set(gca,'Position',[0.09 0.19 0.89 0.75])
    end
    plot(px(a).rnd) % chains
    hold on
    % indicate burn-in length:
    for b = 1:NSTAGES-1
        plot(In1(b)*[1 1],pl(a,:),'k--','LineWidth',2)
        plot(In2(b)*[1 1],pl(a,:),'k','LineWidth',2)
    end
    plot(NB*[1 1],pl(a,:),'k','LineWidth',2)
    ylabel(Par_labels{a},UN{:},'Position',[-0.11,0.4, 0],'rotation',0)
    axis([1 NT pl(a,:)])
    set(gca,'ytick',ticks{a})
    xtnn= get(gca,'XTick');
    xtss= get(gca,'XTickLabel');
    xtlf = round(log10(xtnn(end) / str2double(xtss(end,:))));
    if PLOT_AUTC
        % Autocorrelation:
        subplot('Position',[0.84 0.19 0.14 0.75])
        bar((0:maxlag)',px(a).autocorr,'FaceColor',gray,'EdgeColor',gray)
        axis([0 maxlag 0 1])
        set(gca,'xtick',[0 25 50])
    end
    % SAVE IMAGE AND UPDATE TEX FILE:
    cap = sprintf(['Markov chains and autocorrelations of free ',...
        'parameters. Acceptance rate: %5.3f. TickLabelFactor $10^%d$.'],...
        px(a).acr, xtlf);
    s = sprintf('%s%sSET_%d_Chain_%d',p_fig,FS,SET,a);
    s = ts_print(s,FIGTYPE,300);
    ts_WriteLaTeXstringsNxM(fid,s,[NPAR NPAR 1 a],cap,0.99);
    close
    cntF = cntF+1;
    cFIG{cntF} = s;
end
end % if PLOT_TRAC
%
% Remove Outliers (rogue chain parts to low prob. regions):
if REMV_OUTL
    % Out: abs(x_i - MEDIAN(X)) > 4*1.4826*MAD(X)
    PXn = clean_outliers(px,NB);
else
    PXn = px;
end
%
% get new limits and ticks for next plots:
[pl,ticks] = get_lims_ticks(PXn,lb,ub,FREE,'robust','no',NB,[0.01 0.99]);
%
% Histograms:
if PLOT_HIST % if plot wanted
for a = 1:NPAR
    MM = PXn(a).rnd(NB+1:end,:);
    MM = MM(:);
    MM(isnan(MM)) = [];
    ts_fig(1.5,1.5,VIS);
    [ybin,xval] = hist(MM, round(sqrt(numel(MM))));
    bar(xval,ybin,'FaceColor',gray,'EdgeColor',gray)
    hold on
    ybmx = max(ybin);
    ymax = 1.05*ybmx;
    if PRIOR_TYPE(a)=='u' % plot prior dirstibution (uniform)
        plot([lowB(a) lowB(a) uppB(a) uppB(a)],[0 ybmx ybmx 0],'r--')
    else % normal prior
        xPR = linspace(pl(a,1),pl(a,2),50)';
        yPR = normpdf(xPR,lowB(a),uppB(a));
        yPR = yPR/max(yPR)*ybmx;
        plot(xPR,yPR,'r--')
    end
    plot(MODE(a)*[1 1],[0 ymax],'k','LineWidth',1.5) % best model value
    plot(px(a).Tmean*[1 1],[0 ymax],'k--','LineWidth',1.5) % mean
    plot(px(a).Tperc(1)*[1 1],[0 ymax],'k:','LineWidth',.5) % 2.5perc.
    plot(px(a).Tperc(2)*[1 1],[0 ymax],'k--','LineWidth',.5) % 16perc.
    plot(px(a).Tperc(3)*[1 1],[0 ymax],'k','LineWidth',.5) % 50perc.
    plot(px(a).Tperc(4)*[1 1],[0 ymax],'k--','LineWidth',.5) % 84perc.
    plot(px(a).Tperc(5)*[1 1],[0 ymax],'k:','LineWidth',.5) % 97.5perc.
    set(gca,'Position',[0.08 0.32 0.84 0.675])
    xlabel(Par_labels{a},UN{:},'Position',[0.5,-0.19, 0]);
    axis([pl(a,:) 0 ymax])
    set(gca,'ytick',[])
    %set(gca,'xtick',ticks{a})
    set(gca,'TickDir','out','ticklength',2.5*get(gca,'ticklength'))
    
    % SAVE IMAGE AND UPDATE TEX FILE:
    cap =['Model parameter chain histograms after burn-in. ',...
        'Percentiles (2.5, 16, 50, 84, 97.5) are shown as thin ',...
        '(dotted, dashed, solid, dashed, dotted) lines, respectively. ',...
        'Mode and mean are shown as thick solid and dashed lines, ',...
        'respectively. Red dashed lines indicate prior densities.'];
    s = sprintf('%s%sSET_%d_Hist_%d',p_fig,FS,SET,a);
    s = ts_print(s,FIGTYPE,300);
    ts_WriteLaTeXstringsNxM(fid,s,[NPAR 3 5 a],cap,0.95);
    close
    cntF = cntF+1;
    cFIG{cntF} = s;
end
end % if PLOT_HIST
%
% PLOT: maximum likelihood envelope over parameters
%
if PLOT_ENVL % if plot wanted
N1 = NB;%round(NT/2);
N2 = NT-N1;
N3 = N2*NC;
pdf_v = (reshape(PXn(1).pdf(N1+1:end,:),N3,1));
pdf_v = exp(pdf_v);%-min(pdf_v));
pdf_lims = [min(pdf_v) max(pdf_v)];
pdf_d = diff(pdf_lims);
yminmax = [pdf_lims(1)-0.02*pdf_d pdf_lims(2)+0.02*pdf_d];
for a = 1:NPAR % dots: probability
    ts_fig(1.5,1.5,VIS);
    
    vv   = reshape(PXn(a).rnd(N1+1:end,:),N3,1); % par.values
    IN   = isnan(vv);
    vv(IN) = []; % remove nan from parameter values
    pdf_vn = pdf_v;
    pdf_vn(IN)=[]; % remove nan from prob.values
    vvmm = [min(vv) max(vv)]; % min,max par.values
    vvv  = linspace(vvmm(1),vvmm(2),101); % 100 segments over par.range
    vvmp = nan(1,100); 
    for b = 1:100 % get max.probability of each segment
        ii = find(vv>vvv(b)&vv<=vvv(b+1));
        vvsm = max(pdf_vn(ii));
        if ~isempty(vvsm), vvmp(b) = vvsm; end
    end
    vvpx = vvv(1:100)+0.5*(vvv(2)-vvv(1)); % plot prob.values here
    vvmp = smooth(vvmp,9,'loess'); 
    % sample density:
    % openGL set by ndhist, but cannot print properly, so use zbuf
    set(gcf, 'renderer', 'zbuffer');
    ndhist(vv,pdf_vn,'axis',[pl(a,:) yminmax]);
    %plot(vv, pdf_v, 'k.'); % black dots for all pdf-value pairs
    hold on
    plot(vvpx,vvmp,'color',gray,'LineWidth',1) % max.prob. envelope!
    %plot(vvx,f) % kernel density estimate plot
    plot(MODE(a)*[1 1],yminmax,'r','LineWidth',1) % best model value
    plot(px(a).Tperc(2)*[1 1],yminmax,'r--','LineWidth',1) % 16perc.
    plot(px(a).Tperc(4)*[1 1],yminmax,'r--','LineWidth',1) % 84perc.
    set(gca,'Position',[0.05 0.30 0.84 0.695])
    xlabel(Par_labels{a})%,UN{:},'Position',[0.5,-0.19, 0]);
    axis([pl(a,:) yminmax])
    set(gca,'ytick',[])
    set(gca,'xtick',ticks{a})
    set(gca,'TickDir','out')
    
    % SAVE IMAGE AND UPDATE TEX FILE:
    cap =['Sampling density of unscaled posterior probability ',...
        'versus the corresponding parameter values (after ',...
        'burn-in). The gray line represents a ',...
        'maximum probability envelope over the parameter ranges. ',...
        'For comparison, the red lines are identical to the black ',...
        'lines of the previous figure.'];
    s = sprintf('%s%sSET_%d_PDFs_%d',p_fig,FS,SET,a);
    s = ts_print(s,FIGTYPE,300);
    ts_WriteLaTeXstringsNxM(fid,s,[NPAR 3 5 a],cap,0.95);
    close
    cntF = cntF+1;
    cFIG{cntF} = s;
end
end % if PLOT_ENVL
% infodump:
if WRIT_INF1 % if text wanted
fprintf(fid,['\n\nParameters: %d, Chains: %d, ',...
    'each chain: Burn-In: [%d x %d]=%d, ',...
    'Length after burn-in: %d\n\n'],...
    NPAR,NC,NSTAGES,NBs,NB,NL);
end % if WRIT_INF1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TABLE: Posterior statistics
inpar = X0;
inpar(~FIX) = MODE;
Nh   = size(inpar,1);
%
if WRIT_TAB1 % if text wanted
Ncol = size(inpar,2);
cALIGN = cell(Ncol,1);
for a=1:Ncol % alignment string (different colors for FIX and ~FIX)
    if FIX(a), sAL='>{\color{red}}c'; else sAL='c'; end
    cALIGN{a,1} = sAL;
end
align = sprintf('%s',cALIGN{:});
% prepare special sub/super-script table strings:
cBOD = cell(1,Ncol);
for a = 1:Ncol
    fm = format_out{a};
    if ~FREE(a) % indicate fixed with asterisk
        cBOD{a} = sprintf([fm '*'], inpar(a));
    else % write mode,97.5 & 2.5 perc as sub/super-script
        vv = [MODE(a); PXn(FREE(a)).Tperc([5 1])];
        cBOD{a} = sprintf(['$' fm '^{' fm '}_{' fm '}$'],vv);
    end
end
%put(fid,'\subsection{Best-fitting Model}');
cap = ['Optimal model parameter values given with 2.5th percentile ',...
    '(subscript) and 97.5th percentile (superscript). ',...
    'Red color and * indicate fixed values.'];
put(fid,'\tabcolsep=0.11cm'); % make tables narrower (global in doc)
tsLatexTable(fid,header_out, cBOD, [], align, cap); % TAB
end % if WRIT_TAB1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT DATA AND GMM ATTENUATION CURVES
if PLOT_ATTC
% Plot both for rock (left) and stiff soil (right)
Styp = {'Rock', 'SSoil'}; % fig name part
StypAnnot = {'rock', 'stiff soil'}; % plot annotation
% convert data logtype/units to m/s^2, no log:
Ymss = convert_GMM_units(yy,SMOD.LOGT,SMOD.UNIT,'','m');
XM = SDAT.IMRSH(:,2); % Mw list
XR = SDAT.IMRSH(:,3); % Rjb list
XS = SDAT.IMRSH(:,4); % Site number list (typically 0 or 1)
XH = SDAT.IMRSH(:,5); % hypocentral depth list
% prepare GMM simulation for some meaningful MW and distance
Mmax = max(M); % MW 1 to simulate
Mmin = min(M); % MW 2 to simulate
Mmv = [Mmin Mmax];
Mmid = (Mmax+Mmin)/2;
vRS = [0 1]; % values meaning rock and soil
Ndsim = 25; % points to simulate over distance range
Rmin = max(min(R)/2, 0.3); % min distance to simulate
Rmax = max(R)*2; % max distance to simulate
l10R = [log10(Rmin), log10(Rmax)];
R = logspace(l10R(1), l10R(2), Ndsim)'; % distances to simulate
H = mean(H)*ones(Ndsim,1); % depth to simulate
% get total model standard deviance (max likelihood):
if NSIG == 2 % if two sigmas, combine 
    sigT = sqrt(inpar(iZ(1))^2 + inpar(iZ(2))^2);
else
    sigT = inpar(iZ(1));
end
sigT = convert_GMM_units(sigT,SMOD.LOGT,'','',''); % remove log from sigma
% set plot limits
Rrng = l10R(2)-l10R(1);
Rrng = max(fix(Rrng*2+0.99)/2, 3);
Rlim = 10.^((l10R(2)+l10R(1))/2 + Rrng/2*[-1 1]);
Ymax = max(Ymss);
% colors blue-green-red (all light) from min(Mw) to max(Mw)
Mmm = [Mmin; Mmid; Mmax];
Ccc = [0.5 0.5 1; 0.5 1 0.5; 1 0.5 0.5];
for iS = 1:2 % soil types (Fig1 and optional Fig2)
    soilType = vRS(iS);
    S = soilType*ones(Ndsim,1); % soil types to simulate
    ts_fig(5.4,5.15,VIS);
    set(gca,'position',[0.105 0.105 0.875 0.875])
    % Data points:
    for b=1:Nm % each earthquake gets color for Mw
        iI = EvID==b & XS==soilType; % rock/soil
        if any(iI)
            cmw = unique(XM(iI)); % current Mw
            mcol = interp1(Mmm,Ccc,cmw);
            lh=loglog(XR(iI),Ymss(iI),'d',...
                'markeredgecolor',mcol,...
                'markerfacecolor',mcol,'markersize',7);
            ll1 = lh;
            hold on
        end
    end
    % Plot calibrated GMM:
    for b=1:2 % at two magnitudes
        M = Mmv(b)*ones(Ndsim,1); % MW to sim
        gmm_fun = eval(['@(c)   ' sf.SIM]); % updated function handle
        Ysim = gmm_fun(inpar);
        Ysim = convert_GMM_units(Ysim,SMOD.LOGT,SMOD.UNIT,'','m');
        lhA = loglog(R,Ysim,'-','color','k','linewidth',1.5);
        lhB = loglog(R,Ysim*sigT,'--','color','k','linewidth',1.0);
        loglog(R,Ysim/sigT,'--','color','k','linewidth',1.0);
        hold on
    end
    LEGSTR = {'Data','GMM','\pm 1 \sigma'};
    ll = [ll1; lhA; lhB];
    tsSquareLogAxes(gca, Ymax*2, Rlim)  % fix plot extent and shape
    lgd = legend(ll,LEGSTR,'location','sw');
    set(lgd,'color','none')
    text(0.05,0.2,StypAnnot{iS},UN{:})
    xlabel('$\mathrm{R_{JB} [km]}$','interpreter','latex')
    ylabel('$\mathrm{SA [m/s^2]}$','interpreter','latex')
    set(gca,'xgrid','on','ygrid','on')
    % add COLORBAR for MW of data:
    x = linspace(Mmm(1),Mmm(3),64)';
    CMap3 = interp1(Mmm,Ccc,x);
    axes('position',[0.85 0.76 0.02 0.2]);
    imagesc([1 1],[Mmm(1) Mmm(3)],x)
    colormap(CMap3)
    set(gca,'YDir','normal')
    set(gca,'Yaxislocation','right')
    set(gca,'xtick',[])
    text(-1.5,-0.08,'Data Mw','units','normalized')
    %
    s = sprintf('%s%sSET_%d_Attenuation_%s', ...
        p_fig,FS,SET,Styp{iS});
    s = ts_print(s,'png',300);
    cap = sprintf(['Attenuation curves for the lowest (%3.1f) and ',...
        'highest (%3.1f) magnitudes ',...
        'of the inferred best (ML) solution. Input ',...
        'data (diamonds) colorcoded by magnitude. ',...
        'Data and model curves for rock (left) and stiff soil ',...
        '(right).'],Mmv);
    ts_WriteLaTeXstringsNxM(fid,s,[2 1 2 iS],cap,0.99);
    close
    cntF = cntF+1; cFIG{cntF} = s;
end
M = SDAT.IMRSH(:,2); % Mw list
R = SDAT.IMRSH(:,3); % Rjb list
S = SDAT.IMRSH(:,4); % Site number list (typically 0 or 1)
H = SDAT.IMRSH(:,5); % hypocentral depth list
end % if PLOT_ATTC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SCATTER PLOTs
Mcv = nan(NC*NL,NPAR);
labels2 = cell(2,NPAR);
for a = 1:NPAR
    Mcv(:,a) = reshape(px(a).rnd(NB+1:end,:),NL*NC,1);
    p=regexp(Par_labels{a},' \['); % if given, put unit on second line
    if isempty(p)
        labels2(:,a) = {Par_labels{a};''};
    else
        labels2(:,a) = {Par_labels{a}(1:p-1);Par_labels{a}(p+1:end)};
    end
end
% get indices for best models:
llik = reshape(px(1).pdf(NB+1:end,:),NL*NC,1);
[~,Ix]=sort(llik);
Ib1   = Ix(end); % best model
% correlation
Rh = corr(Mcv,'type',CORR_TYPE);
%
% figure, plot:
if PLOT_COR1 % if plot wanted
ts_fig(7,6.5,VIS);
% openGL set by ndhist, but cannot print properly, so use zbuf
set(gcf, 'renderer', 'zbuffer');
Nsub  = NPAR-1;
yLxs  = 0.02-0.10*Nsub; % ylabel position shift in x-y direction,...
yLys  = 0.5 -0.06*Nsub; %   functions based on trial-error.
mar_h = [0.075 0.020]; % bottom top
mar_w = [0.100 0.025]; % left right
gap   = [0.012 0.015]; % hor. vert.
hy = (1 - sum(mar_h) - (Nsub-1)*gap(2))/Nsub;
ytop = 1-mar_h(2);
wx = (1 - sum(mar_w) - (Nsub-1)*gap(1))/Nsub;
xs = wx+gap(1);

for a = 1:Nsub
    y0 = ytop - a*hy - (a-1)*gap(2);
    for b = 1:a
        x0 = mar_w(1) + (b-1)*xs;
        pos = [x0 y0 wx hy];
        subplot('Position',pos)
        [~,density,X,Y] = kde2d([Mcv(:,b),Mcv(:,a+1)]); % kde2d
        h = pcolor(X,Y,density); hold on
        set(h,'edgecolor','none');
        colormap(CMap)
        axis([pl(b,:) pl(a+1,:)])
        plot(Mcv(Ib1,b),Mcv(Ib1,a+1),'rd','MarkerSize',6,pLW{:}) % best mod
        text(0.97,0.03,sprintf('%5.2f',Rh(b,a+1)),...
            UN{:},HAlr{:},'VerticalAlignment','Bottom')
        set(gca,'ytick',[]);
        if b == 1 % first column: show ytick and ylabel
            set(gca,'ytick',ticks{a+1})
            ylabel(labels2(:,a+1), 'rotation',0,...
                UN{:},'Position',[yLxs, yLys, 0])
        end
        set(gca,'xtick',[]);
        if a == Nsub % last row: show xtick and xlabel
            set(gca,'xtick',ticks{b})
            xlabel(Par_labels{b}, UN{:}, 'Position', [0.5, -0.22, 0])
        end
    end
end
% SAVE IMAGE AND UPDATE TEX FILE:
cap =sprintf(['Posterior probability 2D histograms indicating ',...
    'model parameter correlations as ' CORR_TYPE '''s correlation ',...
    'coefficient in the lower right of each plot. The red diamond ',...
    'represents the best fitting model.']);
s = sprintf('%s%sSET_%d_Scatter',p_fig,FS,SET);
s = ts_print(s,FIGTYPE,300);
ts_WriteLaTeXstringsNxM(fid,s,[1 1 1 1],cap,0.85);
close
cntF = cntF+1;
cFIG{cntF} = s;
end % if PLOT_COR1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulate (MLE), get BIAS and Random Effect ETAi
mBMS = nan(Nh,2);
SAobs = LogY(:,SET);
SAsim = fh.SIM(inpar);
% Get Bias, Slopes(R,M), Random Effect eta_i
bias = SAobs-SAsim;%log10(SAobs./SAsim);
BMS2 = [mean(bias(:)) std(bias(:))]; % total bias over entire freq. range
if NSIG == 2 % if random effects model, calculate ETAi:
    [etai,mag] =  random_effect_etai(...
        SAobs,SAsim,vNS,MRS,inpar(iZ(1))^2,inpar(iZ(2))^2);
else
    etai = nan(Nm,1);
    mag  = Mw;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CORRELATION MATRIX
if PLOT_COR2 % if plot wanted
ts_fig(6,5.5,VIS);
imagesc(Rh,[-1 1])
axis image
format_ticks(gca,Par_labels,Par_labels,[],[],[],[],[0.06 0.01]);
colormap(CMap2);colorbar
set(gca,'Position',[0.08 0.00 0.77 0.99])
% SAVE IMAGE AND UPDATE TEX FILE:
cap =sprintf(['Left: ' CORR_TYPE '''s correlation coefficient matrix ',...
    'for posterior model parameters. Right: Bias versus distance ',...
    '(top) and Mw (middle), event term $\\eta_i$ (bottom).']);
s = sprintf('%s%sSET_%d_CorrMat',p_fig,FS,SET);
s = ts_print(s,FIGTYPE,300);
ts_WriteLaTeXstringsNxM(fid,s,[2 1 2 1],cap,0.98);
close
cntF = cntF+1;
cFIG{cntF} = s;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT BIAS AND EVENT TERMS
Irock = S==0; % logical index of rock sites
%
CC = cell(3,5);
% Slope: Bias - log10(Rjb)
[abc,slo,~,~,s_ab,XYarr,~]= ts_lsq(log10(MRS(:,2)),bias,0.05,'-q');
CC{1,1} = abc; CC{1,2} = slo; CC{1,3} = s_ab;
CC{1,4} = squeeze(XYarr(:,5,:)); % conf. interval y-values (:,[1 2])
CC{1,5} = 10.^(XYarr(:,1,1)); % conf. interval x-values (:)
% Slope: Bias - Mw
[abc,slo,~,~,s_ab,XYarr,~]= ts_lsq(MRS(:,1),bias,0.05,'-q');
CC{2,1} = abc; CC{2,2} = slo; CC{2,3} = s_ab;
CC{2,4} = squeeze(XYarr(:,5,:)); % conf. interval y-values (:,[1 2])
CC{2,5} = XYarr(:,1,1); % conf. interval x-values (:)
if NSIG == 2 % if random effects model,
    % Slope: RandomEff. Eta - Mw
    [abc,slo,~,~,s_ab,XYarr,~]= ts_lsq(mag,etai,0.05,'-q');
    CC{3,1} = abc; CC{3,2} = slo; CC{3,3} = s_ab;
    CC{3,4} = squeeze(XYarr(:,5,:)); % conf. interval y-values (:,[1 2])
    CC{3,5} = XYarr(:,1,1); % conf. interval x-values (:)
end
%
ts_fig(6,5.5,VIS);
subplot('Position',[0.12 0.75 0.86 0.24])
semilogx(LIMr,[0 0],'k:'); hold on % zero line
semilogx(MRS(Irock,2),bias(Irock),'k.','markersize',8)
semilogx(MRS(~Irock,2),bias(~Irock),'.','color',gray,'markersize',8)
semilogx(CC{1,5},CC{1,4}(:,1),'k--',CC{1,5},CC{1,4}(:,2),'k--')
axis([LIMr LIMb])
set(gca,'ytick',-1:0.5:1)
xlabel('R_{jb} [km]',UN{:},'Position',[0.5,-0.17, 0])
ylabel('log_{10}(Obs./Pred.)')

subplot('Position',[0.12 0.41 0.86 0.24])
plot(LIMm,[0 0],'k:'); hold on % zero line
plot(MRS(Irock,1),bias(Irock),'k.','markersize',8)
plot(MRS(~Irock,1),bias(~Irock),'.','color',gray,'markersize',8)
plot(CC{2,5},CC{2,4}(:,1),'k--',CC{2,5},CC{2,4}(:,2),'k--')
axis([LIMm LIMb])
set(gca,'ytick',-1:0.5:1)
xlabel('M_W',UN{:},'Position',[0.5,-0.15, 0])
ylabel('log_{10}(Obs./Pred.)')

if NSIG == 2 % if random effects model,
    subplot('Position',[0.12 0.08 0.86 0.24])
    plot(LIMm,[0 0],'k:'); hold on % zero line
    plot(mag,etai,'kd','markersize',5,'markerfacecolor','k')
    plot(CC{3,5},CC{3,4}(:,1),'k--',CC{3,5},CC{3,4}(:,2),'k--')
    axis([LIMm LIMe])
    xlabel('M_W',UN{:},'Position',[0.5,-0.13, 0])
    ylabel('\eta_i','rot',0,'fontsize',14)
end

% SAVE IMAGE AND UPDATE TEX FILE:
s = sprintf('%s%sSET_%d_Bias',p_fig,FS,SET);
s = ts_print(s,FIGTYPE,300);
ts_WriteLaTeXstringsNxM(fid,s,[2 1 2 2],cap,0.98);
close
cntF = cntF+1;
cFIG{cntF} = s;
end % if PLOT_COR2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVE POSTERIOR STATISTICS
% {per,mode,2.5%,50.0%,97.5%,mean,stdev,corr,etai}
tperc = [PXn(:).Tperc];
OUTstats{SET,1} = xPer; % period
OUTstats{SET,2} = MODE; % mode
OUTstats{SET,3} = tperc(1,:); %  2.5 percentile
OUTstats{SET,4} = tperc(3,:); % 50.0 percentile
OUTstats{SET,5} = tperc(5,:); % 97.5 percentile
OUTstats{SET,6} = [PXn(:).Tmean]; % mean
OUTstats{SET,7} = [PXn(:).Tstd]; % standard deviation
OUTstats{SET,8} = Rh; % correlation matrix
OUTstats{SET,9} = etai'; % event term with magnitude,
%                                 (same order as input data)
OUTstats{SET,10} = mag'; % event magnitudes (corresp. to ETAi)
end % for SET = 1:NSETs; % for each chosen period...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% END OF LOOP OVER PERIODS %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVE ALL STATISTICS TO FILE
if RUN_MODE(2) == 1
    save(sOUTstats,'OUTstats'); % save full statistics cell
    % zip data archive:
    zip(sARCHIVE,cSAVED(:)); % add all SET files to zip archive
end
mass_safe_delete(cSAVED); % remove all those files
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% WRITE STATISTICS TABLES
put(fid,'\chapter{APPENDIX}','');
% collect values of all periods in matrices:
per = [OUTstats{:,1}]';
moz = cell2mat(OUTstats(:,2)); % Mode (MLE)
p02 = cell2mat(OUTstats(:,3)); %  2.5%
p50 = cell2mat(OUTstats(:,4)); % 50.0%
p97 = cell2mat(OUTstats(:,5)); % 97.5%
mez = cell2mat(OUTstats(:,6)); % Mean
sdz = cell2mat(OUTstats(:,7)); % SD
etz = cell2mat(OUTstats(:,9)); % ETAi
ccc = {moz, p02, p50, p97, mez, sdz};
mag = OUTstats{1,10};
% Output Table Header and Number Format: (LaTeX)
hdr = [{'Period'} header_out];
fmt = [{'%5.3f'}  sum_tab_fmt];
% table captions:
cap = {'Maximum likelihood model parameter values at all periods.';
    '2.5 percentile model parameter values at all periods.';
    '50.0 percentile (median) model parameter values at all periods.';
    '97.5 percentile model parameter values at all periods.';
    'Mean model parameter values at all periods.';
    'Standard deviations of model parameters at all periods.'};
for a = 1:numel(cap)
    tsLatexTable(fid, hdr, [per ccc{a}], fmt, [], cap{a}); % TAB
end
% Event term table: (goes off page if too many events!)
hdr1 = cell(1,Nm);
for a = 1:Nm, hdr1{a}=['$M_' num2str(a) '$']; end
hdr = [ {'Period'} hdr1 ; {'$\text{[sec]}$'} num2cell(mag)];
cap = 'Event terms for input earthquakes at all periods';
tsLatexTable(fid, hdr, [per etz], '%5.3f', [], cap); % TAB
%
% save mat file: periods and recalibrated coeff.:
xOUT2 = [per moz];
if RUN_MODE(2) == 1
    save(sOUTtab4,'xOUT2'); % save onle periods and MLE coeff.
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT SUMMARY RESIDUALS
% calc and store all periods' predicted amplitudes:
cSAsim = {nan(Ny,NPER), nan(Ny,NPER)};
for a=1:NPER
    cSAsim{1}(:,a) = fh.SIM(moz(a,:)); % maxlik pred
    cSAsim{2}(:,a) = fh.SIM(p50(a,:)); % median pred
end
% make obs and pred units non-logarithmic:
if strcmpi(SMOD.LOGT,'LN')
    SAobs = exp(LogY);
    cSAsim{1}= exp(cSAsim{1});
    cSAsim{2}= exp(cSAsim{2});
else
    SAobs = 10.^(LogY);
    cSAsim{1}= 10.^(cSAsim{1});
    cSAsim{2}= 10.^(cSAsim{2});
end
% calculate residuals (log10)
cRES = cell(2,1);
for a=1:NPER
    cRES{1} = log10(SAobs./cSAsim{1}); % maxlik pred res
    cRES{2} = log10(SAobs./cSAsim{2}); % median pred res
end
%
% PLOT BIAS OVER ALL PERIODS, R, M
ResYLIM = {0.5*[-1 1]; 1.1*[-1 1]}; % BIAS|residual Y axis limits
fXLIM = [0.3 30]; % Frequency plot limits
vFreq = 1./Per;
Rstri = 'JB'; % TO DO: distance type should be metadata input!
BFRM_Ssettings = struct();
BFRM_Ssettings.VIS = VIS;
BFRM_Ssettings.R_LIM = ResYLIM{1};
BFRM_Ssettings.X_LIM = fXLIM;
BFRM_Ssettings.myRenderer = 'painters';
PtEst = {'MaxLik','Median'}; % point estimate type for file names
EstTypTxt = {'maximum likelihood','median'};
if any(isinf(vFreq))  % if have PGA, i.e. period = 0 s, avoid "Inf" freq
    MaxF = max(vFreq(~isinf(vFreq)));
    sFreqUsed = sprintf('[%4.2f to %5.2f Hz and PGA]',vFreq(end),MaxF);
else
    MaxF = max(vFreq(~isinf(vFreq)));
    sFreqUsed = sprintf('[%4.2f to %5.2f Hz]',vFreq(end),MaxF);
end
for a=1:2  % plot each for max.likelihood and median
    [vFH,~,BMS2] = plot_bias_vs_FRM(...
        MRS,cRES{a},Ny,vFreq,Nm,BFRM_Ssettings);
    cTyp = {'BIAS','RSLO','MSLO'};
    cap = sprintf([...
        '(a) Bias of residuals $\\pm$ one standard error ',...
        'versus frequency %s using \\emph{%s} estimates ',...
        'of model parameters to the %d earthquakes using ',...
        '%d records. The gray shaded area represents the 90 \\%% ',...
        'confidence interval of the bias. ',...
        'Total bias: %.3f, standard deviation $\\sigma$: %.3f. ',...
        '(b) and (c), slopes of lines fitted to residuals versus ',...
        '$\\log_{10}(R_{' Rstri '})$ and $M_W$, respectively,  $\\pm$ ',...
        'one standard error for each frequency. '],...
        sFreqUsed,EstTypTxt{a},Nm,Ny,BMS2);
    % save figures, update tex-file:
    for b = 1:3
        s = sprintf('%s%sBias_%s_%svsF', p_fig,FS,PtEst{a},cTyp{b});
        s = ts_print(vFH(b),s,FIGTYPE,300);
        ts_WriteLaTeXstringsNxM(fid,s,[3 1 3 b],cap,0.97);
        cntF = cntF+1;
        cFIG{cntF} = s;
    end
    close all
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FULL RESIDUAL PLOTs
RMR_Ssettings = struct();
RMR_Ssettings.VIS = VIS;
RMR_Ssettings.R_LIM = ResYLIM{2};
RMR_Ssettings.myRenderer = 'painters';
cTyp = {'M','R'};
for a=1:2  % plot each for max.likelihood and median
    [vFH,~] = plot_residuals_vs_MR(MRS,cRES{a},vFreq,RMR_Ssettings);
    cap = sprintf([
        'PSA residuals using \\emph{%s} model parameter estimates '...
        'with 95\\%% confidence intervals for a least-squares line fit '...
        'versus $M_W$ (top) and $\\log_{10}(R_{JB})$ (bottom).'],...
        EstTypTxt{a});
    % save figures, update tex-file:
    for c=1:2
        s = sprintf('%s%sResiduals_%s_Bvs%s', p_fig,FS,PtEst{a},cTyp{c});
        s = ts_print(vFH(c),s,FIGTYPE,300);
        ts_WriteLaTeXstringsNxM(fid,s,[2 2 1 c],cap,0.99);
        cntF = cntF+1;
        cFIG{cntF} = s;
    end
    close all
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% WRITE OUTPUT ASCII FILES
%
% Write ASCII table file: [periods and MLE coefficient values]
fid2 = fopen(sOUTtab3,'w'); % writing calibrated coeff. to ascii table file
fprintf(fid2,[sprintf('%s ',fmt{:}) '\n'],[per ccc{1}]');
fclose(fid2);
% (BH's ASCII files removed)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FINISH AND TIDY UP
%
% closing the file
fclose(fid);
close all
%
% Create a PDF from result figures and tables with a LaTeX file:
try
    % compile and process latex files to pdf file:
    pdflatex(fhead,option_pdflatex);
    % move pdf to pdf-folder and rename, move aux files away:
    movefile(fpdf,fpdf2);
catch
    warning(['It looks like the pdflatex command does not work ',...
        'properly and the result PDF file cannot be created. ',...
        'If you do have have LaTeX and want to use it, ',...
        'it is free open source software and you can search for it, ',...
        'for example this site might help: ',...
        'https://www.latex-project.org/. If you do have a LaTeX ',...
        'compiler and something is still wrong, you might need more ',...
        'packages. Move up the head tex file from the aux directory, ',...
        'and run "pdflatex head_*tex" in your system command line. ',...
        'Missing packages will be indicated. After installing all ',...
        'needed LaTeX packages, the Matlab script should work as ',...
        'well. If it still does not work, remove this "try" frame ',...
        'and run again without catching errors to investigate.'])
end
movefile(['head_' fbase '.*'],[p_aux '/'],'f');
%
% zip figure archive:
cFIG=cFIG(~cellfun(@isempty,cFIG)); % remove empty cells
zip(sFIGARCH,cFIG(:));  % add all figure files to zip archive
mass_safe_delete(cFIG); % remove all those files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%