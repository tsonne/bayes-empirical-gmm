function ts_WriteLaTeXstringsNxM(...
    fid,epsname,figconfig,caption,width,height)
%
% SYNTAX
%     WriteLaTeXstringsNxM(fid,epsname,figconfig,caption);
%
% This scripts writes LaTeX code to a text file that organizes plots
% into LaTex figures with N rows and M columns of plots in each figure.
% Thus each figure contains NxM plots, and if NS is larger than NxM, it
% writes out multiple figures.
%
% fid       = file identifier
% epsname   = file name of plot
% figconfig = [NS N M n] where
%             NS is the total number of plots.
%             N is the total number of rows in the figure,
%             M is the number of columns,
%             n is this figure's number
% caption   = the caption for the figure (only printed when n=NxM)
% width     = relative width of figure field to page (use 0.9 mostly)
%
% 2016-10-22 tsonne: Added 'height' option, can now set maximal height on
%                    page with 1 being total text height.
%                    Almost complete rewrite: No sprintf() in fprintf()
%                    anymore!
%

I=find(epsname=='\');
if ~isempty(I), epsname(I)='/'; end
if ~exist('width','var'), width = 0.99; end
if ~exist('height','var'), height = 0.99; end

% WRITING THE FILE PARTS
NS= figconfig(1);
N = figconfig(2);
M = figconfig(3);
n = figconfig(4);

FigW=(1/M)*width;
FigH=(1/N)*height;
z = sprintf(['width=%4.2f\\textwidth,'...
    'height=%4.2f\\textheight,keepaspectratio'],FigW,FigH);

% Accounts for overflow
norig=n;
n = mod(n,N*M);

if all(figconfig==1)  % when only one figure, [1 1 1 1]
   
   str={'\begin{figure}[H]',...
        '  \centering',...
       ['  \includegraphics[' z ']{' epsname '}']};
   fprintf(fid,'%s\n',str{:});
   str={['  \caption{' caption '}\label{}'],...
        '\end{figure}','',''};
   fprintf(fid,'%s\n',str{:});
   return
end

if n==1
   
   str={'\begin{figure}[H]',...
        '  \centering',...
       ['  \includegraphics[' z ']{' epsname '}']};
   fprintf(fid,'%s\n',str{:});
   
elseif n==0 && norig~=NS
   
   str={['  \includegraphics[' z ']{' epsname '}\\'],...
        ['  \caption{' caption '}\label{}'],...
         '\end{figure}','',''};
   fprintf(fid,'%s\n',str{:});
   
else
   
   if mod(n,M)~=0
      str={['  \includegraphics[' z ']{' epsname '}']};
      fprintf(fid,'%s\n',str{:});
   elseif mod(n,M)==0
      str={['  \includegraphics[' z ']{' epsname '}\\']};
      fprintf(fid,'%s\n',str{:});
   end
   
end

if norig==NS
   str={['  \caption{' caption '}\label{}'],...
         '\end{figure}','',''};
   fprintf(fid,'%s\n',str{:});
end