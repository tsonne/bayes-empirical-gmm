function SMOD = extract_coeff(S)
% PURPOSE: For list of model definitions, set prior type, lower/upper
% boundaries, start value, normal prior coefficients based on table for
% all chosen periods and return struct with correct settings.
%
% 2020-05-12 Tim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
per = S.per; % (vector) Chosen periods to process
PCS = S.PCS; % (matrix) Coefficient table with at least as many periods,
             % used as input for normal prior assumptions
sf_x = S.sf_x; % (cell) function definition strings
mod_n = S.mod_n; % (cell) model name strings
PRI = S.PriorType; % (cell) prior PDF types for each coefficient and model
M2SD = S.RatioMean2SD; % (float) to set coefficients with normal prior:
                       % mean is taken from PCS, stdev is abs(mean/M2SD)
X0 = S.X0; % (vector) start values
lowB = S.lowB; % (vector) lower boundaries
uppB = S.uppB; % (vector) upper boundaries
LOGT = S.LOGT; % (string) data expected in this log (log10 or log)
UNIT = S.UNIT; % (string) data expected in this unit (m/s^2, cm/s^2, g)
%
% find desired periods 'per' in coefficient table 'PCS':
if isempty(per) || (ischar(per) && strcmpi(per,'all'))
    per = PCS(:,1); % take all periods
    id = find(per>=0); % only take positive, -1 exists as PGV
    per = PCS(id,1);
else
    [A,id]=ismember(per,PCS(:,1));
    assert(all(A),'Check period list, some are not available!!!')
end
NPAR2 = numel(X0); % number of model parameters to use for Bayesian Inf
NPAR1 = NPAR2 - 2; % N of par from prior table (no standard deviations)
Cv = PCS(id,2:1+NPAR1); % coefficient values from table
Nper = numel(per);
% Fill in normal prior values based on coeff. table for each set and period
Ngmm = numel(sf_x);
cX0 = cell(Ngmm,1);
cLB = cell(Ngmm,1);
cUB = cell(Ngmm,1);
for a=1:Ngmm
    % assign values for current model 'a'
    xX0 = repmat(X0,Nper,1); % all periods' start values same
    xLB = nan(Nper,NPAR2);
    xUB = nan(Nper,NPAR2);
    pr = lower(PRI{a});
    iu = pr=='u';
    in = find(~(iu(1:NPAR1)));
    for b=1:Nper % for each period of model 'a'
        xLB(b,iu) = lowB(1,iu);
        xUB(b,iu) = uppB(1,iu);
        if ~isempty(in) % if any normal priors
            xLB(b,in) = Cv(b,in); % mean
            xUB(b,in) = abs(Cv(b,in)/M2SD); % SD
        end
    end
    % store in cells (each element = one model)
    cX0{a}=xX0;
    cLB{a}=xLB;
    cUB{a}=xUB;
end
% collect models into a struct:
SMOD = repmat(...
    struct('name','','sf','','X0',[],'lowB',[],'uppB',[]),...
    Ngmm,1);
for a=1:Ngmm
    SMOD(a).name = mod_n{a};
    SMOD(a).sf   = sf_x{a};
    SMOD(a).X0   = cX0{a};
    SMOD(a).lowB = cLB{a};
    SMOD(a).uppB = cUB{a};
    SMOD(a).RND_EFFECT='event'; % use random effects model
    SMOD(a).Prior= PRI{a};
    SMOD(a).per = per;
    SMOD(a).LOGT = LOGT;
    SMOD(a).UNIT = UNIT;
end
%