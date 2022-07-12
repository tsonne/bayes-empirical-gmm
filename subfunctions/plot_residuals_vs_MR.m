function [vFH,cap] = plot_residuals_vs_MR(MRS,resd,vFreq,Ssettings)
% Plot two figures: Residuals vs Magnitude, and Residuals vs log10(R)
%
% INPUT:
%   MRS = matrix: [Mw(ns,1) R(ns,1) S(ns,1)];   "S" not needed
%   resd = matrix of residuals: (ns,nf)
%   vFreq = vector of frequencies for resd: (nf,1)
%   Ssettings = struct of plot settings, see below in code
% OUTPUT:
%   vFH = vector figure handle for both figures
%   cap = suggested figure caption when combining them vertically as one
%
% 2019-01-29 Tim: created to make scripts shorter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<3, Ssettings=struct(); end
% figure visibility, 'on'|'off'
VIS = getpar(Ssettings,'VIS','on');
% plot limits for residual value y-axis
R_LIM = getpar(Ssettings,'R_LIM',[-1.1 1.1]);
% distance-type string, e.g. 'JB', 'rup' ...
Rstri = getpar(Ssettings,'Rstri','JB');
% y-axis label string, e.g. 'log_{10}(Obs./Pred.)'
Ylbl = getpar(Ssettings,'Ylbl','Residual');
% figure renderer: 'opengl'|'painters'
myRenderer = getpar(Ssettings,'myRenderer','opengl');
% plot marker option cell
Mk = getpar(Ssettings,'Mk',...
    {'kd','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize'});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UN={'Units','Normalized'};
MRS(MRS(:,2)<0.1,2) = 0.1; % R_JB = 0 must not cause issues with log()
M_rng = max(MRS(:,1))-min(MRS(:,1));
MR_LIM = [min(MRS(:,1))-M_rng*0.1 max(MRS(:,1))+M_rng*0.1; ...
    0.1    1.1*max(MRS(:,2))];
MR_TIC = {3:10;  [0.3 1 3 10 30 100 300 1000]};
if M_rng<2.2
    MR_TIC{1} = 3:0.5:10;
end
cXL = {'M_W',['R_{' Rstri '} [km]']};
YlabelSet = {Ylbl,UN{1},UN{2},'Position',[-0.14, 0.5, 0],'fontsize',10};
col_CI = [.5 .5 .5]; % color of confidence interval lines
% PLOT bias:
Nper = numel(vFreq);
nco = min(Nper,7);
nro = ceil(Nper/nco);
ii = 1:Nper;
iYL = ii(mod(ii-1,nco)==0);
vFH = nan(2,1); % figure handles
MargLow = 0.405-0.03*nro; % bottom position of lowest subplot
MargLeft = 0.05;
HGap = 0.0265 -0.0022*nro;
for c = 1:2 % 1: Mw, 2: log10(R)
    
    vFH(c) = ts_fig(9,1.2*nro,VIS);
    set(gcf,'Renderer', myRenderer);
    ha = tight_subplot(nro,nco,[HGap .005],[MargLow .01],[MargLeft .005]);
    for a = 1:Nper
        
        set(gcf,'CurrentAxes',ha(a))
        plot( MR_LIM(c,:),[0 0],'k-'); hold on % zero line
        plot( MRS(:,c),resd(:,a),Mk{:},3)
        % Slope: residuals vs Mw|log10(rJB)
        if c==1, MRx=MRS(:,c); else MRx=log10(MRS(:,c)); end
        [~,~,~,~,~,XYarr,~]= ts_lsq(MRx,resd(:,a),0.05,'-q');
        conyy = squeeze(XYarr(:,5,:)); % conf. interval y-values
        conx = XYarr(:,1,1); % conf. interval x-values
        if c==2, conx=10.^conx; end
        plot(conx,conyy(:,1),'color',col_CI)
        plot(conx,conyy(:,2),'color',col_CI)
        if c==2
            set(ha(a),'xscale','log','yscale','linear')
        end
        axis(ha(a), [MR_LIM(c,:) R_LIM])
        if isinf(vFreq(a))
            f_annot = 'PGA';
        else
            f_annot = sprintf('%5.2f Hz',vFreq(a));
        end
        text('Parent',ha(a), UN{:}, 'Position',[0.025 0.92],...
            'String',f_annot,'fontsize',12)
        set(gca,'XTick',MR_TIC{c});
    end
    if nro>1
        set(ha(1:(nro-1)*nco),'XTickLabel','');
    end
    set(ha(ii(mod(ii-1,nco)>0)),'YTickLabel','')
    for b=1:numel(iYL)
        ylabel(ha(iYL(b)), YlabelSet{:})
    end
    for a = (nro-1)*nco+1:Nper
        xlabel(ha(a),cXL{c},UN{:},'Position',[0.5,-0.23, 0],'fontsize',12);
    end
    if nro*nco>Nper % hide superfluous axes
        for a = 1:(nro*nco - Nper)
            set(ha(Nper+a), 'Visible','off')
        end
    end
end
% figure caption (if shown combined)
cap = ['PSA residuals at each frequency for all records with 95\% '...
    'confidence intervals for a least-squares line fit versus $M_W$ '...
    '(top) and $\log_{10}(R_{' Rstri '})$ (bottom).'];
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