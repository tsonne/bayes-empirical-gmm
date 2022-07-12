function tsSquareLogAxes(thisax,FAS,XLIM,YMAG)
%
% This function ensures that log-log axes appear square so that functions
% of f^n appear as straight lines with slope n
%
% Syntax: tsSquareLogAxes(thisax,FAS,XLIM,YMAG)
% Input: thisax = axis handle
%        FAS    = Fourier amplitude spectrum
%        XLIM   = (optional) set x-axis limits (Default: 10.^[-1.5 2.5])
%        YMAG   = (optional) y-axis magnitude range, e.g. 4 to show exactly
%                 4 orders of magnitude as in [0.01 100].
%                 OR: YLIM vector, e.g. [0.01 100] to set. 
%                 (Default: same as x-axis, square plot) 
%                 If not same as x-axis, plot area will not
%                 be square, but x- to y-axis unit lengths will be same, 
%                 such that f^n has slope n on plot.
%
% tsonne 2016-09-17
%        2016-11-02 just corrected description above.
%        2017-02-24 update to allow YMAG setting.
%        2018-04-13 allow YLIM vector instead of YMAG
if nargin<3, XLIM = 10.^[-1.5 2.5]; end
m = log10(max(FAS));  % max Y value exp
% upper Y plot limit exp
% if ceil(m/0.5)-m/0.5 < 0.01
%     % if very close to a 0.5 log10 step, go 0.5 up
%     y = ceil(m/0.5)*0.5+0.5;
% else
%     y = ceil(m/0.5)*0.5;
% end
y = ceil(m/0.5)*0.5;
x = log10(XLIM);      % x-axis exp
rx = x(2)-x(1);        % x value range (to be same on y-axis)
if nargin<4
    axis(thisax, [XLIM 10.^[(y-rx) y]])
    set(thisax, 'xtick',10.^(ceil(x(1)):1:floor(x(2))))
    set(thisax, 'ytick',10.^ceil((y-rx-1):1:y))
    axis(thisax, 'square')
else
    if numel(YMAG)==1
        ry = YMAG;
    else
        y = log10(YMAG(2));
        ry = y-log10(YMAG(1));
    end
    axis(thisax, [XLIM 10.^[(y-ry) y]])
    set(thisax, 'xtick',10.^(floor(x(1)):1:ceil(x(2))))
    set(thisax, 'ytick',10.^ceil((y-ry-1):1:y))
    set(thisax, 'PlotBoxAspectRatio', [rx ry ry])
end