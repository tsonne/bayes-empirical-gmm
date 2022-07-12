function h = plot_mag_dist(M,R,S,show,save)
% PLOT magnitude to distance distribution of data.
% Indicate whether rock or soil by two different symbols.
%
% INPUT:
%  M  = (vector) magnitude
%  R  = (vector) distance
%  S  = (vector) site type, expect 0 and 1 values only
%  show = (str) 'on' or 'off' to show figure or hide it and close figure.
%         Default is 'on'.
%  save = (str) output figure file path. Default is no save.
%
% OUTPUT:
%  h  = figure handle (empty if show=='off')
%
% 2020-05-12 Tim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
rr = [0 max(R)]; % min max distances
mm = [min(M) max(M)]; % min max magnitudes
SiteTypes = [0, 1];
gr1 = .5*[1 1 1]; % dark grey
Mk = {'d','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6;
      's','MarkerEdgeColor','k','MarkerFaceColor',gr1,'MarkerSize',6};
if nargin<4
    show = 'on';
end
h = ts_fig(3,3,show); hold on
for a=1:2
    b=S==SiteTypes(a);
    if ~any(b), continue; end
    plot(R(b),M(b),Mk{a,:})
end
XL2 = 10*ceil(rr(2)/10); % right x-axis limit
xlim([0   XL2 ])
% ylim: include next lower Mw if min(Mw) close to floor:
YLM = [floor(mm(1)) ceil(mm(2))]; 
if mm(1)-YLM(1) < 0.5
    YLM(1) = YLM(1)-1;
end
ylim(YLM)
%
Tf = {'interpreter','latex','FontSize',14,'units','normalized'};
xlabel('$\mathrm{R_{JB}}$ [km]',Tf{:},'position',[0.5 -0.09 0])
ylabel('$\mathrm{M_{W}}$','rotation',0,Tf{:},'position',[-.14 0.54 0])
grid on; box on
%set(gca,'Position',[0.18 0.145 0.81 0.835]) % use full space for plot
set(gca,'Position',[0.18 0.145 0.79 0.835]) % more space for last xtick
%
if nargin>4 && any(save)
    x_print(save,300);
end
if strcmpi(show,'off')
    close(h)
    h=[];
end