function x_print(filename,resolution,h)
% Save figure as eps, pdf (default), or png.
%
% INPUT:
%   filename   - filepath, filetype deduced from extention. Filetypes can
%                be 'eps','pdf','png'. DEFAULT: 'pdf' if no extension.
%   resolution - dots-per-inch output resolution;
%                  [optional, DEFAULT: 300]
%   h          - figure handle; [optional, DEFAULT: current figure]
% OUTPUT:
%   None
% 
% EXAMPLES:
%   x_print('fig/test1.eps',200);
%   x_print('fig/test2.pdf',300,h);
%
%
% 2020-05-12 Tim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Check input:
assert(nargin>0,'Need at least filepath input!');
assert(ischar(filename),'Filepath input must be string!');
if nargin<2, resolution=300; end
if nargin<3, h=gcf; end
[~,~,ext] = fileparts(filename);
if ~any(ext)
    ext='.pdf';
    filename = strcat(filename, ext);
end
%
switch ext
    case '.eps'  % PRINT TO .EPS
        print(h, '-depsc', sprintf('-r%d', resolution), filename);
        %fix_lines(finalpath); % improves dashed and dotted lines in EPS,
                               % must get function from MatlabCentral by
                               % Oliver Woodford
    case '.pdf'  % PRINT TO .PDF
        print(h, '-dpdf', sprintf('-r%d', resolution), filename);
    case '.png'  % PRINT TO .PNG
        print(h, '-dpng', sprintf('-r%d', resolution), filename);
    otherwise
        error('Filetype must be in {pdf,eps,png}!');
end