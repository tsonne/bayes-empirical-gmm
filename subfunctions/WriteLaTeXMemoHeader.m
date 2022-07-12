function WriteLaTeXMemoHeader(fid,Subject,Date,TexFile2)
% PURPOSE: Write the first part of a LaTeX memo.
%
% INPUT:
%   fid = file identifier of a tex file
%   Subject = titlepage title
%   Date = desired date to write under title
%   TexFile2 = name of second tex file to be included as body
%
% 2020-09-11 Tim Sonnemann
%   Based on old version from 2014,
%   rewritten to use public LaTeX documentclass instead of special one.
%
I=find(TexFile2=='\'); % to fix issues with MSWindows path backslash
if ~isempty(I), TexFile2(I)='/'; end % substitute with '/'
str={
    '\documentclass[a4paper,10pt,oneside]{book}'
    '\usepackage{graphicx}'
    '\usepackage{epstopdf}'
    '\usepackage[reqno]{amsmath}'
    '\usepackage{float}'
    '\usepackage{amsfonts}'
    '\usepackage{amssymb}'
    '\usepackage{lscape}'
    '\usepackage{natbib}'
    '\usepackage{longtable}'
    '\usepackage{hyperref}'
    '\usepackage[utf8]{inputenc}'
    '\usepackage[T1]{fontenc}'
    '\usepackage{booktabs,caption,pdflscape,afterpage}'
    '\usepackage{xcolor,array,siunitx}'
    '\setcounter{secnumdepth}{2}'
    '\usepackage[margin=1in]{geometry}'
    '\usepackage{titlesec}'
    '\titlespacing*{\chapter}{0pt}{-30pt}{20pt}'
    '\titleformat{\chapter}[block]{\normalfont\huge\bfseries}{\thechapter.}{20pt}{\huge}'
    '\pagestyle{plain}'
    '\setlength{\parindent}{0pt}'
    '\setlength\abovecaptionskip{5pt}'
    '\setlength\belowcaptionskip{10pt}'
    ''
    '\begin{document}'
    '\begin{titlepage}'
    '    \centering'
    ['    {\huge\bfseries ' Subject ' \par}']
    '    \vspace{1cm}'
    ['    {\large ' Date ' \par}']
    '%    \vspace{1cm}'
    '%    {\large to \par}'
    '%    \vspace{1cm}'
    '%    {\large from \par}'
    '\end{titlepage}'
    '\tableofcontents'
    '\newpage'
    ['\input{' TexFile2 '}']
    '\end{document}'};
 fprintf(fid,'%s\n',str{:});