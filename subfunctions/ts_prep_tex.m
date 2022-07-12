function [p_tex,p_fig,p_pdf,p_aux,fbase,fhead,fbody,fpdf,fpdf2] = ...
    ts_prep_tex(MEMONAME,Date,Subject)
%
%%% SETTING LATEX PRINT FILE
% Prepare DIRs, variable file names, write MemoHeader
% 2015-06-11 Tim Sonnemann: original version
% 2020-09-10 TS: use new WriteLaTeXMemoHeader, cut others,
%   remove inputs 'VERSION', 'From', 'To', simplified things.
%
% prepare folders for figures, tex, auxiliary and pdf files:
p_tex = 'tex'; % tex files path
p_fig = 'fig'; % figures file path
p_pdf = 'pdf'; % pdf files path
p_aux = 'aux'; % latex auxiliary files (compiler output)
[~,~,~] = mkdir(p_fig);
[~,~,~] = mkdir(p_tex);
[~,~,~] = mkdir(p_pdf);
[~,~,~] = mkdir(p_aux);
% file names:
FS = filesep;
fbase = MEMONAME; % base name
fhead = sprintf('%s%shead_%s.tex',p_tex,FS,fbase); % header tex file
fbody = sprintf('%s%sbody_%s.tex',p_tex,FS,fbase); % body tex file
fpdf  = sprintf('head_%s.pdf',fbase); % created pdf
fpdf2 = sprintf('%s%s%s.pdf',p_pdf,FS,fbase); % moved/renamed pdf
% write tex header:
fid=fopen(fhead,'w');
WriteLaTeXMemoHeader(fid,Subject,Date,fbody);
fclose(fid);