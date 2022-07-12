function finalpath = ts_print(varargin)
% Save figure as eps, pdf (default), both, or png.
%
% SYNTAX:  finalpath = ts_print(figpath,filetype,resolution)
%     OR:  finalpath = ts_print(fhandle,figpath,filetype,resolution)
%
% INPUT: 
%   fhandle    - figure handle; [optional, DEFAULT: current figure]
%   figpath    - filepath without filename extention (added automatic.);
%   filetype   - image file format (string) 'eps','pdf','both','png';
%                  [optional, DEFAULT: 'pdf']
%   resolution - dots-per-inch output resolution;
%                  [optional, DEFAULT: 300]
% OUTPUT:
%   finalpath  - output filepath with extension
% 
% EXAMPLES:
%   finalpath = ts_print('fig/test1','eps',200);
%   finalpath = ts_print(h,'fig/test2','pdf',300);
%
%
% ts 2015-10-07; mod 2015-12-31: direct to pdf with adjusted ts_fig()
% 2016-01-30 tsonne: added PNG option
% 2017-02-20 tsonne: Added optional figure handle input (like the original
%   print-function), such that a new syntax alternative exists (see top).
%   Also updated description, input checks and PNG handling. 
%   Unsolved Linux problem could be completely erased, as its only purpose
%   was automatic trimming of white image borders when first printing to
%   EPS, then transmogrifying to different formats (minor convenience).
% 2019-02-23 Tim: added another optional argument to varargin:
%   can set PROBLEM to 0 (Default:1) to force ghostscript eps trim.
%
%PROBLEM = 1; % if 0, use linux epstopdf to get trimmed pdf file
% Check input:
nv = numel(varargin);
assert(nv>0,'Need at least filepath input!');
if ischar(varargin{1}) % no handle provided
    h = gcf;
    figpath    = varargin{1};
    filetype   = getpar(varargin,2,'pdf');
    resolution = getpar(varargin,3,300);
    PROBLEM = getpar(varargin,4,1);
else
    h = varargin{1};
    assert( ishandle(h) && strcmp(get(h,'type'),'figure'), ...
        'Invalid figure handle (or non-char filename)!');
    assert(nv>1,'Need at least filepath input!');
    figpath = varargin{2};
    assert(ischar(figpath),'Filepath input must be string!');
    filetype   = getpar(varargin,3,'pdf');
    resolution = getpar(varargin,4,300);
    PROBLEM = getpar(varargin,5,1);
end
if ~any(strcmpi(filetype,{'eps','pdf','both','epspdf','png'}))
    error('Filetype must be in {pdf,eps,both,png}!');
end
%
% PRINT TO .EPS
if any(strcmpi(filetype,{'eps','both','epspdf'}))
    finalpath = strcat(figpath, '.eps');
    print(h, '-depsc', sprintf('-r%d', resolution), finalpath);
    %fix_lines(finalpath); % improves dashed and dotted lines in EPS
end
%
% PRINT TO .PDF
% convert the .eps encapsulating the fonts
if any(strcmpi(filetype,{'pdf','both','epspdf'}))
    if PROBLEM == 0 % no problems with system call (see below)
        % if pdf chosen, will produce eps and convert, then remove eps.
        if ~exist(strcat(figpath, '.eps'),'file')
            finalpath = strcat(figpath, '.eps');
            print(h, '-depsc', sprintf('-r%d', resolution), finalpath);
        end
        finalpath = strcat(figpath, '.pdf');
        [~,~] = system( sprintf( ...
            ['epstopdf --gsopt=-dPDFSETTINGS=/prepress ',...
            '--outfile=%s.pdf %s.eps'], ...
            figpath,figpath));
        %
        % if .eps unwanted, remove
        if ~any(strcmpi(filetype,{'both','epspdf'}))
            delete(sprintf('%s.eps', figpath));
        end
    elseif PROBLEM == 1 % problem with epstopdf, but print('dpdf'... ok
        finalpath = strcat(figpath, '.pdf');
        print(h, '-dpdf', sprintf('-r%d', resolution), ...
            finalpath);
        
    else
        error(['Value of ''PROBLEM'' (4th input w/o figure handle) ',...
            'must be either 0 or 1'])
        % problem on Debian 9 with ghostscript 9.16
        % epstopdf returns error:
        % Error: /undefined in --definefont--
        % ...
        % no solution found, maybe bug in gs 9.16,
        % but works perfectly when run in terminal (not through matlab).
        % For now, just produce eps and let pdflatex do the conversion
        % to pdf, which works even though pdflatex is also called by matlab
    end
end
%
% PRINT TO .PNG
% done this way to have automatic trim:
%
% ghostscript error:
% Error: /execstackoverflow in /findfont
% Operand stack:
%    50   Helvetica   ISOLatin1Encoding   Helvetica   Helvetica
% ...
% gs 9.16
if strcmpi(filetype,'png')
    if PROBLEM == 0
        finalpath = strcat(figpath, '.png');
        print(h, '-depsc', sprintf('-r%d', resolution), [figpath '.eps']);
        system( sprintf(...
            'convert -trim -density %d %s.eps %s.png',...
            resolution,figpath,figpath))
        delete(sprintf('%s.eps', figpath));
    else
        finalpath = strcat(figpath, '.png');
        print(h, '-dpng', sprintf('-r%d', resolution), [figpath '.png']);
    end
end
end
%%%%%%%%%%
% Get the value from input cell, or use default value:
function par   = getpar(var,k,defval)
if numel(var)>k-1
    par = var{k};
else
    par = defval;
end
end