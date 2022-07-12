function h = ts_fig(W,H,visible,scrW,scrH)
% Input W and H in inches, scrW and scrH in pixels (screen size!)
% If the axis is supposed to be square, use
% the ratio W/H = 5/4.75 due to problems with legends.
% For example, figr(5,4.75), or figr(4,3.8), or figr(4.21,4);
%
% original: B.Halldorsson.
% 2015-08-10 tsonne: visibility option added, clean up
% 2015-12-31 tsonne: papersize and position adjusted for pdf print
% 2016-01-30 tsonne: can give figure handle as 'visible'-input, then no new
%                    figure is created, but given figure is refitted.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SCR = [900 600];
vis = 'on';
% check function input values:
if nargin<1
    W=5;H=4;
elseif nargin<2 && length(W)==2
    H=W(2);W=W(1);
elseif nargin < 3
    % vis
elseif nargin < 4
    vis = visible;
else
    SCR = [scrW scrH];
    vis = visible;
end
if ischar(vis)
    switch lower(vis)
        case {'on','yes','y'}
            vis = 'on';
        case {'off','no','n'}
            vis = 'off';
        otherwise
            vis = 'on';
    end
    h = figure('Visible',vis); % create figure
else
    h = vis; % use given figure handle
end

% set reasonable(?) screen position
hA = findall(0,'type','figure');
Nf = numel(hA)-1; % number of open figures minus the current one
asprat = W/H;
if asprat>1
    SCRpos = [50+Nf*20 50+Nf*20 SCR(1) SCR(1)/asprat];
else
    SCRpos = [50+Nf*20 50+Nf*20 SCR(2)*asprat SCR(2)];
end
% alter figure to fitting papersize
set(h,'PaperUnits','inches',...
    'PaperSize',[W H],...
    'position',SCRpos,...
    'PaperPosition',[0 0 W H]);