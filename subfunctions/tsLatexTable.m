function tsLatexTable(fid,header,body,format,align,caption,label,varargin)
% WRITE INPUT AS TABLE TO LATEX FILE.
% INPUT: 
%  fid        = file identifier
%  **header   = first row (has \hline below, 'w1 & w2' or {'w1','w2'}) ['']
%  **body     = table values (can be matrix or cell)                   ['']
%  *format    = row format, e.g. '%d & %f' or {'%d','%f'}    [auto-num2str]
%               (auto-expands single given format for multiple columns)
%  *align     = column alignment format, e.g. ' l | cc r'         [all 'c']
%               (auto-expands single given alignment for multiple columns)
%  *caption   = table caption string                                   ['']
%  *label     = table label string                                     ['']
%  *varargin  = can be up to 3 options:
%           'literal' will detokenize text, 
%           'sci'     will use sci.notat. numbers.
%           'secure'  will only use given format with matching data types,
%                     i.e. if body has string but format is numeric, do not
%                     convert string to number and vice versa. Useful when
%                     providing mixed body cell with complex formatting.
%        These three extra settings work best when giving a format string,
%        as opposed to using the automatic formatting.
% ENTRIES WITH (*) ARE OPTIONAL, IF (**) AT LEAST ONE OF BOTH REQUIRED.
% EMPTY SETTINGS WILL BE AUTOMATICALLY SET TO DEFAULTS IN [].
% 
% LaTeX Packages required:
%   \usepackage{longtable} % if >50 lines long table
%   \usepackage{siunitx}   % if 'sci' option used
%   \usepackage{caption}   % for \captionof{table}{} command
%
%
% 2016-01-31 tsonne: created
% 2016-03-30 tsonne: header can have multiple lines (use cell, Nrows>1)
% 2016-07-05 tsonne: header can be mixed (string and double) in cell.
% 2016-10-06 tsonne: added option 'secure' for using mixed string-numeric
%                    body cells, converting only numbers by given format.
% 2017-02-27 tsonne: uses longtable if >50 lines.
% 2017-07-12 tsonne: 'secure' option now works also when table columns
%                    have input formats like '%04d%02d%02d'. (bugfix)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK ALL INPUT VARIABLES, ASSIGN DEFAULT IF EMPTY
assert(isvalidfid(fid),'Given fid is not valid');
FLAG = [1 1 1 1 1 1 0 0 0];
cINP = {'header'   , '{}' ;
        'body'     , '{}' ;
        'format'   , '{}' ;
        'align'    , '{}' ;
        'caption'  , '{}' ;
        'label'    , '{}' };
for a=1:size(cINP,1)
    if ~exist(cINP{a,1},'var') || isempty(eval(cINP{a,1}))
        eval([cINP{a,1} '=' num2str(cINP{a,2}) ';']); FLAG(a)=0;
    end
end
nv = numel(varargin); c = 1;
while c <= nv
    op = varargin{c};
    if ~isempty(op)
        switch lower(op)
            case 'literal', FLAG(7) = 1;
            case 'sci'    , FLAG(8) = 1;
            case 'secure' , FLAG(9) = 1;
            otherwise     , fprintf('Warning: Unknown option: %s',op);
        end
    end
    c = c+1;
end
% CHECK CONSISTENT INPUT, AUTOGENERATE IF NECESSARY
S4 = '\\\\'; NL = '\\\\\n';
assert(numel(body)+numel(header)>0,'Both table header and body are empty');
HD = header;
if FLAG(1) % if header exists...
    if ischar(HD) % if string: check ending
        hc = numel(strfind(HD,'&'))+1;
        HD = regexprep(HD,'\\\\\\\\\\n$',S4);
        if isempty(regexp(HD,'\\\\[ ]*$','once')), HD = [HD '\\']; end
        HD = {HD};
    else % make string
        [hr,hc] = size(HD);
        HF = repmat('%s & ',1,hc);
        HH = cell(hr,1);
        HD = cellfun(@num2str,HD,'UniformOutput',false);
        for a = 1:hr, HH{a} = sprintf([HF(1:end-2) S4],HD{a,:}); end
        HD = HH;
    end
else
    hc = 0; % no header
end
bc = size(body,2);
BF = format;
if FLAG(2) % if body exists...
    if ~iscell(body),  body = num2cell(body); end
    if FLAG(3) % if format given...
        if ischar(BF)
            nc = numel(strfind(BF,'&'))+1; % table columns
            if nc==1, BF={BF}; end % will expand if nec.
        end
        if iscell(BF)
            if numel(BF)==1, BF=repmat(BF,1,bc); end % expand
            BF = sprintf('%s & ',BF{:});
            BF = [BF(1:end-2) NL];
        end
        if ~strcmp(BF(end-5:end),NL),  BF = [BF NL]; end
        nf = numel(regexp(BF,'%[0-9*.-+# ]*[cdeEfgGisu]')); % req.inputs
        nc = numel(strfind(BF,'&'))+1; % table columns
        assert(bc==nf,'Format - body mismatch!');
        if FLAG(9) % if 'secure' option wanted...
            sf = BF(1:end-6); ns=numel(sf); % analyse format
            si = strfind(sf,'&'); % separation at each '&'...
            ii = [1 si+1; si-1 ns];
            cc = cell(1,nc);
            for a = 1:nc % get number of format inputs per table column
                ni = numel(regexp(sf(ii(1,a):ii(2,a)),...
                    '%[0-9*.-+# ]*[cdeEfgGisu]'));
                cc{a} = repmat('%s',1,ni); % new format after conversion
            end
            sf = regexprep(sf,'&',''); % remove '&'
            ns = numel(sf);
            fi = strfind(sf,'%'); % get individual format rules
            ii = [1 fi(2:end); fi(2:end)-1 ns];
            cf = cell(1,nf); % format cell: original, e.g. ' %02d km'
            cs = cell(1,nf); % format cell: to string, e.g.' %s km'
            NF = true(1,nf);
            for a = 1:nf % determine for each if numeric or string format
                cf{a} = sf(ii(1,a):ii(2,a)); % store format
                NF(a) = isempty(regexp(cf{a},'%[0-9-]*s','once'));
                cs{a} = regexprep(cf{a},'%[0-9*.-+# ]*[cdeEfgGisu]','%s');
            end
            br = size(body,1);
            for a = 1:br % for each line
                for b = 1:bc % for each column
                    % NF(b) = true if numeric format required
                    isNV = ~ischar(body{a,b}); % true if num.value
                    if (isNV && NF(b)) || (~isNV && ~NF(b))
                        body{a,b} = sprintf(cf{b},body{a,b});
                    elseif isNV && ~NF(b) % numeric but not num-format
                        body{a,b} = num2str(body{a,b});
                    else % require numeric but have string
                        body{a,b} = sprintf(cs{b},body{a,b});
                    end
                end
            end
            BF = sprintf('%s & ',cc{:}); BF = [BF(1:end-2) NL];
        end
    else % if no format given...
        body = cellfun(@num2str,body,'UniformOutput',false);
        BF = repmat('%s & ',1,bc); BF = [BF(1:end-2) NL];
        nc = bc; % as many columns as input colummns
    end
    if bc>0 && hc>0, assert(nc==hc,'Header - body mismatch!'); end
    if FLAG(7) % detokenize to print literal text
        BF = regexprep(BF,'(%[0-9-]*s)','\\\\texttt{\\\\detokenize{$1}}');
    end
    if FLAG(8) % enclose number fields with \num{} [package: siunitx]
        BF = regexprep(BF,'(%[0-9*.-+# ]*[deEfgGiu])','\\\\num{$1}');
    end
    BD = body'; % fprintf reads by columns, not rows
end
if ~FLAG(4), align = repmat('c',1,nc); end
na = numel(regexp(align,'(?<!{\\*\w*)[lcrpmbS](?!\w*})')); % alignment can
%  contain l,c,r,p{},m{},b{}; but also modifiers like >{\color{red}}
if na==1 && nc>1 % if one align given, expand to all cols:
    align=repmat(align,1,nc);
    na = nc;
end
assert(na==nc,'Alignment - columns mismatch!');
nr = 0; % total number of rows
if FLAG(1), nr = nr + size(HD,1); end
if FLAG(2), nr = nr + size(BD,2); end
%
% WRITE TABLE TO FILE:
if nr>50 % long table
    if FLAG(5), put(fid,['\captionof{table}{' caption '}']);end % caption
    if FLAG(6), put(fid,['\label{' label '}'])             ;end % label
    put(fid,['\begin{longtable}{' align '}\toprule']);          % alignment
    if FLAG(1), put(fid,HD{:},'\midrule \endhead')         ;end % header
    if FLAG(2), fprintf(fid, BF, BD{:})                    ;end % body
    put(fid,'\bottomrule');
    put(fid,'\end{longtable}');
else % fits on one page, encase as floating minipage
    put(fid,'\begin{minipage}{\linewidth}');
    put(fid,'\centering');
    if FLAG(5), put(fid,['\captionof{table}{' caption '}']);end % caption
    if FLAG(6), put(fid,['\label{' label '}'])             ;end % label
    put(fid,['\begin{tabular}{' align '}\toprule']);            % alignment
    if FLAG(1), put(fid,HD{:},'\midrule')                  ;end % header
    if FLAG(2), fprintf(fid, BF, BD{:})                    ;end % body
    put(fid,'\bottomrule');
    put(fid,'\end{tabular}');
    put(fid,'\end{minipage}','');
end
end