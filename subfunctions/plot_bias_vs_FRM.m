function [vFH,cap,BMS2] = plot_bias_vs_FRM(...
    MRS,resd,NsTotal,vFreq,NEV,Ssettings)
% Plot three figures:
%   Bias with 90% conf.intv. and residual stdev vs freq
%   slope of residuals by log10(distance) with stdev vs freq
%   slope of residuals by magnitude with stdev vs freq
%
% INPUT:
%   MRS = matrix: [Mw(ns,1) R(ns,1) S(ns,1)];   "S" not needed
%   resd = matrix of residuals: (ns,nf)
%   NsTotal = total number of used data (receiver-event pairs)
%   NEV  = number of events used
%   vFreq = vector of frequencies for resd: (nf,1)
%   Ssettings = struct of plot settings, see below in code
% OUTPUT:
%   vFH = vector figure handle for all figures
%   cap = suggested figure caption when combining them vertically as one
%   BMS2 = total mean bias and one sigma, (1,2)
%
% 2019-01-29 Tim: created to make scripts shorter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<3, Ssettings=struct(); end
% figure width (inch)
FigW = getpar(Ssettings,'FigW',4);
% figure height (inch)
FigH = getpar(Ssettings,'FigH',2);
% figure visibility, 'on'|'off'
VIS = getpar(Ssettings,'VIS','on');
% plot limits for residual value y-axis
R_LIM = getpar(Ssettings,'R_LIM',[-1.1 1.1]);
% plot limits for frequency value x-axis
fXLIM = getpar(Ssettings,'X_LIM',[0.3 30]);
% distance-type string, e.g. 'JB', 'rup' ...
Rstri = getpar(Ssettings,'Rstri','JB');
% figure renderer: 'opengl'|'painters'
myRenderer = getpar(Ssettings,'myRenderer','opengl');
% line width
LW = getpar(Ssettings,'LW',1.5);
% plot position in each figure
PposXvsF = getpar(Ssettings,'PposXvsF',[0.17 0.23 0.8 0.75]);
% frequency to plot PGA at (vFreq=Inf)
f_PGA = getpar(Ssettings,'f_PGA',30);
% Plot type: errorbar or just lines:
ERROR_BAR = getpar(Ssettings,'ERROR_BAR',true);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT all bias and R-M-slope values:
%
AtoZ = 'abcdefghijklmnopqrstuvwxyz';
UN = {'Units','Normalized'};
% calculate bias and slopes from residuals:
BMS2 = [mean(resd(:)) std(resd(:))]; % total bias over entire freq. range
[BiMS,RsMS,MsMS] = calc_bias_slopes(resd,MRS);
% 90% confidence interval of mean bias at each frequency:
meanB = BiMS(:,1); % mean bias at each f
DOF   = NsTotal - 1; % degrees of freedom of t-distribution
std_B = BiMS(:,2); % standard deviance of residuals at each f
ci90  = tinv(0.95,DOF)*std_B/sqrt(DOF); % one-sided interval (total 90%)
%
mBRM = [BiMS RsMS MsMS]; % all bias and slope values in one matrix
cYL3 = {'log_{10}(Obs./Pred.)',...
    ['log_{10}(R_{' Rstri '})-slope'],'M_W-slope'};
vFH = nan(3,1);
vFreq = vFreq(:);
HAVE_PGA = false;
if any(isinf(vFreq))
    HAVE_PGA = true;
    MaxF = max(vFreq(~isinf(vFreq)));
    vFreq(isinf(vFreq)) = f_PGA; % make sure to plot PGA not at infinity
end
if ERROR_BAR
    fEC = 'k';
else
    fEC = 'none';
end
for a = 1:3
    vFH(a) = ts_fig(FigW,FigH,VIS);
    set(gcf,'Renderer', myRenderer);
    if a==1
        X = [vFreq  ; flipud(vFreq)];
        Y = [meanB+ci90; flipud(meanB-ci90)];
        fill(X,Y,0.8*[1 1 1],'edgecolor',fEC); hold on
    end
    if ERROR_BAR
        plot(fXLIM,[0 0],'k--'); hold on % zero line
        errorbar(vFreq, mBRM(:,2*a-1), mBRM(:,2*a),'k','LineWidth',LW)
    else
        plot(fXLIM,[0 0],'k:'); hold on % zero line
        m = mBRM(:,2*a-1); % mean
        s = mBRM(:,2*a); % 1 sigma
        plot(vFreq, m,'k','linewidth',1) % mean
        plot(vFreq, m+s, 'k--', vFreq, m-s, 'k--') % +- 1 sigma
    end
    set(gca,'xscale','log','yscale','linear')
    set(gca,'xtick',[0.5 1 3 10 20],...
        'ytick',-0.4:0.2:0.4,...
        'XMinorTick','off', 'YMinorTick' , 'off')
    set(gca,'ticklength',2*get(gca,'ticklength'))
    axis([fXLIM R_LIM])
    xlabel('SDOF oscillator eigenfrequency [Hz]',...
        UN{:},'Position',[0.5,-0.17, 0],'fontsize',12)
    ylabel(cYL3{a},UN{:},'Position',[-0.125,0.5, 0],'fontsize',12)
    text(-0.2,-0.15,['(' AtoZ(a) ')'],UN{:},'fontsize',14)
    set(gca,'Position',PposXvsF)
    % Create another tick mark axis at same position
    cax2 = axes('Position', get(gca, 'Position'));
    % Set small tick marks, and remove the X / Y TickLabel
    set(cax2 ...
        ,'YAxisLocation' ,'left' ...
        ,'XAxisLocation' ,'Bottom' ...
        ,'Box'           ,'off' ...
        ,'HitTest'       ,'off' ...
        ,'Box'           ,'off' ...
        ,'Color'         ,'none' ...
        ,'XMinorTick'    ,'off' ...
        ,'YMinorTick'    ,'off' ...
        ,'xscale'        ,'log' ...
        ,'XTick'         ,[0.1:0.1:1 2:1:10 20:10:50] ...
        ,'XTickLabel'    ,[] ...
        ,'XLim'          ,fXLIM ...
        ,'YTick'         ,-1:0.1:1 ...
        ,'YTickLabel'    ,[] ...
        ,'YLim'          ,R_LIM ...
        );
end
if HAVE_PGA
    cap = sprintf([...
    '(a) Bias $(B_{best})$ of residuals $\\pm$ one standard error ',...
    '$(\\sigma_{best})$ versus frequency [%4.2f to %5.2f Hz and PGA] ',...
    'of the best-fit model to the %d earthquakes using ',...
    '%d records. The gray shaded area represents the 90 \\%% ',...
    'confidence interval of the bias. ',...
    'Total bias: %.3f, standard deviation $\\sigma$: %.3f. ',...
    '(b) and (c), slopes of lines fitted to residuals versus ',...
    '$\\log_{10}(R_{' Rstri '})$ and $M_W$, respectively,  $\\pm$ ',...
    'one standard error for each frequency. '],...
    vFreq(end),MaxF,NEV,NsTotal,BMS2);
else
    cap = sprintf([...
    '(a) Bias $(B_{best})$ of residuals $\\pm$ one standard error ',...
    '$(\\sigma_{best})$ versus frequency [%4.2f to %5.2f Hz] ',...
    'of the best-fit model to the %d earthquakes using ',...
    '%d records. The gray shaded area represents the 90 \\%% ',...
    'confidence interval of the bias. ',...
    'Total bias: %.3f, standard deviation $\\sigma$: %.3f. ',...
    '(b) and (c), slopes of lines fitted to residuals versus ',...
    '$\\log_{10}(R_{' Rstri '})$ and $M_W$, respectively,  $\\pm$ ',...
    'one standard error for each frequency. '],...
    vFreq(end),vFreq(1),NEV,NsTotal,BMS2);
end

%
end % function end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function y = getpar(S,field,default)
%GETPAR get parameter value from a struct
% S         struct
% field     field to be extracted from the S struct
% default   default value if field is not a member of the S struct
%
if isfield(S,field)
    y = S.(field);
elseif nargin>2
    y = default;
else
    error('Need value for struct field: %s',field);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%