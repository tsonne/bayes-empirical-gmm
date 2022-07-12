function [results,PV,PP] = ts_MSAMv2(model,params,options,INFO,PV,PP,PR);
% MULTI-STAGE ADAPTIVE METROPOLIS
% Adaptation by intra-stage tuning and inter-stage var.analysis
% ATTENTION: USING PARFOR!
%
% 2017-03-08 tsonne: created
%   --> not very good: cov. too often not positive definite!
%       DRAM is much better...
%
results = []; % This is supposed to hold info about the final chain cov...
NC = INFO.NC; NL = INFO.NL; NB = INFO.NB; NT = INFO.NT;
NPAR = INFO.NPAR;
NSTAGES = INFO.NSTAGES; NBs = INFO.NBs; NB0 = INFO.NB0;
In0 = INFO.In0; In1 = INFO.In1; In2 = INFO.In2;
assert(NSTAGES>3,'This scheme needs at least 4 stages!');
% adaptation by intra-stage tuning and inter-stage var.analysis
F.pp = @(x,t,s,r) x+t*s.*r; % Tuned normal proposal function
F.po = model.ssfun; % log posterior
F.sp = model.stfun; % start proposal
s_d  = (2.38)^2/NPAR; % scaling for cov Gelman etal 1996,
%                                       Roberts etal 1997
X0 = params.par0; % initial value;
TuneItv = options.adaptint; % tune interval
iCOV = options.qcov; % initial covariance
SD1 = options.SD0; % initial stdev vector
SD0 = SD1*0.5;
eco = 0.0001*eye(NPAR); % small value to add to proposal cov.matrix
% Go through the burn-in stages:
for stg = 1:2 %NSTAGES-1
    J0 = In0(stg);
    J1 = In1(stg);
    J2 = In2(stg);
    if stg==1
        StX0 = X0;
        StSD = SD0;
        rSD  = SD1;
    else
        StX0 = X0new;
        StSD = SDnew;
        rSD  = SDnew;
    end
    if stg == 3
        % NormProposal by Roberts etal 1997:
        F.pp = @(x,t,r) x + (chol(Pcov)'*r')';
    end
    parfor c = 1:NC
        % Init., Start value (acception guaranteed):
        a = J0;
        acn  = 1; % acceptance counter
        acr  = 0; % acceptance rate
        tune = 1; % tuning factor (adaptive Metropolis)
        x_t = F.sp(StX0,StSD); % propose start values
        p_t = F.po(x_t); % loglik
        PV{c}(a,:) = x_t; % save
        PP{c}(a)   = p_t; % save
        % Chain (2:end):
        for a = J0+1:J2
            r  = PR{c}(a,:); % get saved random numbers
            x  = F.pp(x_t,tune,rSD,r); % propose new values
            p  = F.po(x);       % loglik
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
                tune = adapt_tuning(tune,acr);
                fprintf('W:%d C:%6d AR:%5.3f T:%7.3f\n',...
                    c,a,acr,tune);
            end
        end
    end
    %
    % Stage Analysis:
    NL1 = J2-J1+1;
    w = nan(NL1,NC); % weights of each chain
    for c = 1:NC
        w(:,c) = PP{c}(J1:J2);
    end
    w = exp(w - max(w(:))); % weights shifted, loglik->lik
    aws = nan(NC,NPAR); % all weighted stdevs
    awm = nan(NC,NPAR); % all weighted means
    awc = nan(NPAR,NPAR,NC); % all weighted cov. matrices
    for c = 1:NC
        [WCOV,WMEAN,WSD] = wcov(w(:,c),PV{c}(J1:J2,:));
        awm(c,:) = WMEAN;
        aws(c,:) = WSD;
        awc(:,:,c) = WCOV;
    end
    X0new = median(awm); % new starting point
    SDnew = sqrt(s_d)*median(aws); % est.SD by median of deviations
    WCnew = median(awc,3)+eco; % median of weighted covariance
end
%
%
% NormProposal by Roberts etal 1997, with tuning:
F.pp = @(x,t,c,r) x + t*r*c;%t*(c*r')';
for stg = 3:NSTAGES-1
    J0 = In0(stg);
    J1 = In1(stg);
    J2 = In2(stg);
    StX0 = X0new;
    StSD = SDnew;
    Rnew = sqrt(s_d)*chol(WCnew);%';
    parfor c = 1:NC
        % Init., Start value (acception guaranteed):
        a = J0;
        acn  = 1; % acceptance counter
        acr  = 0; % acceptance rate
        tune = 1; % tuning factor (adaptive Metropolis)
        x_t = F.sp(StX0,StSD); % propose start values
        p_t = F.po(x_t); % loglik
        PV{c}(a,:) = x_t; % save
        PP{c}(a)   = p_t; % save
        % Chain (2:end):
        for a = J0+1:J2
            r  = PR{c}(a,:); % get saved random numbers
            x  = F.pp(x_t,tune,Rnew,r); % propose new values
            p  = F.po(x);       % loglik
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
                tune = adapt_tuning(tune,acr);
                fprintf('W:%d C:%6d AR:%5.3f T:%7.3f\n',...
                    c,a,acr,tune);
            end
        end
    end
    %
    % Stage Analysis:
    NL1 = J2-J1+1;
    w = nan(NL1,NC); % weights of each chain
    for c = 1:NC
        w(:,c) = PP{c}(J1:J2);
    end
    w = exp(w - max(w(:))); % weights shifted, loglik->lik
    aws = nan(NC,NPAR); % all weighted stdevs
    awm = nan(NC,NPAR); % all weighted means
    awc = nan(NPAR,NPAR,NC); % all weighted cov. matrices
    for c = 1:NC
        [WCOV,WMEAN,WSD] = wcov(w(:,c),PV{c}(J1:J2,:));
        awm(c,:) = WMEAN;
        aws(c,:) = WSD;
        awc(:,:,c) = WCOV;
    end
    X0new = median(awm); % new starting point
    SDnew = sqrt(s_d)*median(aws); % est.SD by median of deviations
    WCnew = median(awc,3)+eco; % median of weighted covariance
end
%
%
% Run final stage (no tuning):
%F.pp = @(x,r) x + (chol(WCnew)'*r')';
F.pp = @(x,r) x + r*sqrt(s_d)*chol(WCnew);
parfor c = 1:NC
    a    = J2+1;
    acn  = 1;
    r    = PR{c}(a,:); % get saved random numbers
    x_t  = F.pp(X0new,r); % propose start values
    p_t  = F.po(x_t);   % loglik
    PV{c}(a,:) = x_t; % save
    PP{c}(a)   = p_t; % save
    for a = J2+2 : NT
        r  = PR{c}(a,:); % get saved random numbers
        x  = F.pp(x_t,r); % propose new values
        p  = F.po(x);     % loglik
        if p-p_t > log(rand) % Metropolis step: accept
            p_t = p;
            x_t = x;
            acn = acn+1;
        end
        PV{c}(a,:) = x_t; % save
        PP{c}(a)   = p_t; % save
        if mod(a,TuneItv)==0 % just print iter#,acr
            acr = acn/(a-J2+1); % acc.rate of entire stage
            fprintf('W:%d C:%6d AR:%5.3f\n',c,a,acr);
        end
    end
end
end