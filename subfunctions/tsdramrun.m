function [results,chain,PP] = tsdramrun(model,params,options)
% DRAMRUN  Metropolis-Hastings MCMC run with adaptive delayed rejection (DRAM)
%
% This function generates MCMC chain using DRAM adaptation for a model defined
% by user supplied log-posterior function.
% It is a modified version of Marko Laine's original code, such that likelihood
% and prior function definitions must be combined before use, and therefore
% arbitrary covariances can be used, like a random effects model. The posterior
% log probabilities are output as well. The data must be included in the log-
% posterior function definition.
%
% 
% INPUT:
%
% model.ssfun =  % log likelihood + log prior function string,
%                % that returns  log(p(y|par)p(par)),
%                % example definition:
%                %   data = [some data...]
%                %   x = parameters (here: 1:6=sim model, 7:8=tau,sigma)
%                %   fll = function: get log likelihood
%                %   fsim = function: simulate data using parameters
%                %   fcov = function: create cov matrix from tau,sigma
%                %   fpri = function: get log prior p using fixed hyperparams
%                %   '@(x) fll(data,fsim(x(1:6)),fcov(x(7:8))) + fpri(x,hyper)'
%
% params.par0   =  ; % initial parameter vector (a row vector)
% params.bounds = ;  % 2*npar matrix of parameter bounds
%                    % default: [-Inf,Inf]
%
% options.nsimu  = 2000;   % length of the chain
% options.qcov   = ;       % proposal covariance matrix
%
% parameters for DRAM
% options.adaptint = 10;  % how often to adapt, if zero, no adaptation
% options.drscale  = 3;   % scale for the second proposal, if zero, no DR
%
% OUTPUT:
%
% results  structure that contains some info about the run
% chain    nsimu*npar MCMC chain
% PP       log posterior prob values


% calls covupd.m for covariance update and (optionally) gammar_mt.m for
% gamma variates

% this is a 'simple' version for demonstration and educational purposes

% Original from:
% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.0 $  $Date: $
%
% Modified version from:
% Tim Sonnemann, 2016-11-11
% + Changed input: removed params.sigma2 and model.priorfun,
% + use model.ssfun instead as defining logprior + loglikelihood.
% + Added input options.pfact to allow normalizing parameter values,
% + mostly experimanetal and ideally shouldn't be needed.
% + Changed output: instead of s2chain (was sigma square chain), give PP,
% + PP is log posterior probability value chain.
% + some refactoring.

%% get values from the input structs
nsimu  = getpar(options,'nsimu',1000);
% initial parameter vector
par0   = getpar(params,'par0');
par0 = par0(:)'; % row vector
% number of parameters
npar   = length(par0);
% 2*npar matrix of parameter bounds
bounds = getpar(params,'bounds',(ones(npar,2)*diag([-Inf,Inf]))');
% log posterior probability function
ssfun  = getpar(model,'ssfun');

%%% parameters for DRAM
% how often to adapt, if zero, no adaptation
adaptint = getpar(options,'adaptint',100);
% scale for the second proposal, if zero, no DR
drscale  = getpar(options,'drscale',3);
% scale for adapting the propsal
adascale = getpar(options,'adascale',2.4/sqrt(npar));
% blow factor for covariace update
qcovadj  = getpar(options,'qcovadj',1e-5);
% proposal covariance
qcov = getpar(options,'qcov');

% temporary Parameter scaling Factors for numerical stability
PF = getpar(options,'pfact',[]);

% to DR or not to DR
if drscale<=0, dodr=0; else dodr=1;end

printint  = getpar(options,'printint',500);
verbosity = getpar(options,'verbosity',0);

% IF NOT USING A NUMERICAL STABILITY FACTOR, RUN NORMAL CODE:
if isempty(PF)
    %
    % Cholesky factor of proposal covariance
    R       = chol(qcov + eye(npar)*qcovadj*5); % *adascale;
    if dodr
        R2      = R./drscale; % second proposal for DR try
        iR      = inv(R);
    end
    chain   = zeros(nsimu,npar);  % we store the chain here
    PP      = zeros(nsimu,1);  % store the log posterior prob values here
    
    oldpar       = par0;                % first row of the chain
    oldss        = feval(ssfun,oldpar);% first log posterior
    acce         = 1;                       %  how many accepted moves
    chain(1,:)   = oldpar;
    PP(1)        = oldss;
    %%%
    % covariance update uses these to store previous values
    chaincov = []; chainmean = []; wsum = []; lasti = 0;
    %%% the simulation loop
    for isimu=2:nsimu
        
        % info on every printint iteration
        if isimu/printint == fix(isimu/printint)
            fprintf('isimu=%d, %d%% done, accepted: %d%%\n',...
                isimu,fix(isimu/nsimu*100),fix((acce/isimu)*100));
        end
        
        newpar = oldpar+randn(1,npar)*R;     % a new proposal
        
        accept = 0;
        % check bounds
        if any(newpar<bounds(1,:)) || any(newpar>bounds(2,:))
            newss = Inf;
            newprior = 0;
            alpha12 = 0;
        else % inside bounds, check if accepted
            newss  = feval(ssfun,newpar);   % log posterior
            alpha12 = min(1,exp(newss-oldss));
            if rand < alpha12 % we accept
                accept   = 1;
                acce     = acce+1;
                oldpar   = newpar;
                oldss    = newss;
            end
        end
        if accept == 0 && dodr % we reject, but make a new try (DR)
            newpar2 = oldpar+randn(1,npar)*R2;  % a new try
            
            if any(newpar2<bounds(1,:)) || any(newpar2>bounds(2,:))
                newss2 = Inf;
                newprior2 = 0;
            else % inside bounds
                newss2    = feval(ssfun,newpar2);
                alpha32 = min(1,exp(newss-newss2));
                l2 = exp(newss2-oldss);
                q1 = exp(-0.5*(norm((newpar2-newpar)*iR)^2 - ...
                    norm((oldpar -newpar)*iR)^2));
                alpha13 = l2*q1*(1-alpha32)/(1-alpha12);
                if rand < alpha13 % we accept
                    accept = 1;
                    acce     = acce+1;
                    oldpar   = newpar2;
                    oldss    = newss2;
                end
            end
        end
        chain(isimu,:) = oldpar;
        PP(isimu)      = oldss;
        
        if adaptint>0 && fix(isimu/adaptint) == isimu/adaptint
            % adapt the proposal covariances
            if verbosity, fprintf('adapting\n'); end
            % update covariance and mean of the chain
            [chaincov,chainmean,wsum] = covupd(...
                chain((lasti+1):isimu,:),1, ...
                chaincov,chainmean,wsum);
            lasti = isimu;
            [Ra,is] = chol(chaincov + eye(npar)*qcovadj);
            if is % singular cmat
                fprintf('Warning cmat singular, not adapting\n');
            else
                R = Ra*adascale;
                if dodr
                    R2 = R./drscale;     % second proposal for DR try
                    iR = inv(R);
                end
            end
        end
        
    end
    
    % calculate covariance and mean of the chain
    [chaincov,chainmean,~] = covupd(chain((lasti+1):isimu,:),1, ...
        chaincov,chainmean,wsum);
    
    
    results.class = 'MCMC results';
    results.accepted=acce./nsimu;              % acceptance ratio
    results.mean = chainmean;
    results.cov  = chaincov;
    results.qcov = R'*R;
    results.R = R;
    results.nsimu = nsimu;
    results.drscale = drscale;
    results.adascale = adascale;
    results.adaptint = adaptint;
    
else % IF USING A NUMERICAL STABILITY FACTOR, RUN ADJUSTED CODE:
    PF = PF(:)'; % row vector
    FP = 1./PF; % inverse factor
    bounds(1,:) = bounds(1,:).*PF;
    bounds(2,:) = bounds(2,:).*PF;
    %
    % Cholesky factor of proposal covariance (ADJ)
    R       = chol(qcov.*(PF'*PF) + eye(npar)*qcovadj); % *adascale;
    if dodr
        R2      = R./drscale; % second proposal for DR try
        iR      = inv(R);
    end
    chain   = zeros(nsimu,npar);  % we store the chain here
    PP      = zeros(nsimu,1);  % store the log posterior prob values here
    
    oldpar       = par0.*PF;               % first row of the chain (ADJ)
    oldss        = feval(ssfun,oldpar.*FP);% first log posterior (-ADJ)
    acce         = 1;                      %  how many accepted moves
    chain(1,:)   = oldpar;
    PP(1)        = oldss;
    %%%
    % covariance update uses these to store previous values
    chaincov = []; chainmean = []; wsum = []; lasti = 0;
    %%% the simulation loop
    for isimu=2:nsimu
        
        % info on every printint iteration
        if isimu/printint == fix(isimu/printint)
            fprintf('isimu=%d, %d%% done, accepted: %d%%\n',...
                isimu,fix(isimu/nsimu*100),fix((acce/isimu)*100));
        end
        
        newpar = oldpar+randn(1,npar)*R;     % a new proposal
        
        accept = 0;
        % check bounds
        if any(newpar<bounds(1,:)) || any(newpar>bounds(2,:))
            newss = Inf;
            newprior = 0;
            alpha12 = 0;
        else % inside bounds, check if accepted
            newss  = feval(ssfun,newpar.*FP);   % log posterior (-ADJ)
            alpha12 = min(1,exp(newss-oldss));
            if rand < alpha12 % we accept
                accept   = 1;
                acce     = acce+1;
                oldpar   = newpar;
                oldss    = newss;
            end
        end
        if accept == 0 && dodr % we reject, but make a new try (DR)
            newpar2 = oldpar+randn(1,npar)*R2;  % a new try
            
            if any(newpar2<bounds(1,:)) || any(newpar2>bounds(2,:))
                newss2 = Inf;
                newprior2 = 0;
            else % inside bounds
                newss2    = feval(ssfun,newpar2.*FP); % log post. (-ADJ)
                alpha32 = min(1,exp(newss-newss2));
                l2 = exp(newss2-oldss);
                q1 = exp(-0.5*(norm((newpar2-newpar)*iR)^2 - ...
                    norm((oldpar -newpar)*iR)^2));
                alpha13 = l2*q1*(1-alpha32)/(1-alpha12);
                if rand < alpha13 % we accept
                    accept = 1;
                    acce     = acce+1;
                    oldpar   = newpar2;
                    oldss    = newss2;
                end
            end
        end
        chain(isimu,:) = oldpar;
        PP(isimu)      = oldss;
        
        if adaptint>0 && fix(isimu/adaptint) == isimu/adaptint
            % adapt the proposal covariances
            if verbosity, fprintf('adapting\n'); end
            % update covariance and mean of the chain
            [chaincov,chainmean,wsum] = covupd(...
                chain((lasti+1):isimu,:),1, ...
                chaincov,chainmean,wsum);
            lasti = isimu;
            [Ra,is] = chol(chaincov + eye(npar)*qcovadj);
            if is % singular cmat
                fprintf('Warning cmat singular, not adapting\n');
            else
                R = Ra*adascale;
                if dodr
                    R2 = R./drscale;     % second proposal for DR try
                    iR = inv(R);
                end
            end
        end
        
    end
    
    % calculate covariance and mean of the chain
    [chaincov,chainmean,~] = covupd(chain((lasti+1):isimu,:),1, ...
        chaincov,chainmean,wsum);
    
    
    results.class = 'MCMC results';
    results.accepted=acce./nsimu;              % acceptance ratio
    results.mean = chainmean;
    results.cov  = chaincov;
    results.qcov = R'*R;
    results.R = R;
    results.nsimu = nsimu;
    results.drscale = drscale;
    results.adascale = adascale;
    results.adaptint = adaptint;
    
    % ADJUST PARAMETER VALUES BACK
    chain = bsxfun(@times,chain,FP);
    
end
%%%%%%%%
function y=getpar(options,par,default)
%GETPAR get parameter value from a struct
% options   options struct
% par       parameter value to extract from the struct
% default   default value if par is not a member of the options struct

if isfield(options,par)
    y = options.(par);
elseif nargin>2
    y = default;
else
    error('Need value for option: %s',par);
end
