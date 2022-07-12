function pdflatex(tex,option)
% create pdf from latex file
% runs pdflatex twice to build table of contents.
%
% created: 2014/2015 ts
% 2015-12-18 tsonne: added function to find libgcc on current machine
% 2016-01-30 tsonne: added option for quiet run
% 2017-02-28 tsonne: try to make system independent
if nargin<2, option='';end
%
if isunix
    % must set system's LD_LIBRARY_PATH to this (Matlab doesn't know):
    %lib = '/usr/lib/gcc/x86_64-linux-gnu/4.9';
    lib = get_linux_libgcc; % only linux
    if any(strcmpi(option,{'q','-q','quiet','silent'}))
        cmd = ['export LD_LIBRARY_PATH=' lib ...
            '; pdflatex -interaction=batchmode ' tex];
        [~,~]=system(cmd);
        [~,~]=system(cmd);
    else
        system(['export LD_LIBRARY_PATH=' lib '; pdflatex ' tex]);
        system(['export LD_LIBRARY_PATH=' lib '; pdflatex ' tex]);
    end
else
    if any(strcmpi(option,{'q','-q','quiet','silent'}))
        cmd = ['pdflatex -interaction=batchmode ' tex];
        [~,~]=system(cmd);
        [~,~]=system(cmd);
    else
        system(['pdflatex ' tex]);
        system(['pdflatex ' tex]);
    end
end