function WriteLaTeXmatrix(fid,parname,Matrix,format,use_siunitx)
% write an array/matrix to a tex file, such that the latex interpreter
% can transform it to a typical math matrix with brackets
% fid = file identifier to write to
% parname = string that is put before matrix, "parname = Matrix"
% Matrix = 2d-array of numbers to write
% format (optional) = string format, e.g. '%.3f'
% use_siunitx (optional) = either 1 (use it) or something else (no use),
%     when using the latex package siunitx, this allows following usage: 
%     "\num{1.234e-5}" will become the correct scientific notation in PDF,
%     so this function will wrap each input number with \num{}
% tsonne 2015-10-13,
% upd. 2015-10-26: added use_siunitx option
min_arg = 3;
assert(nargin>=min_arg,'This function requires 3 inputs');
assert(isvalidfid(fid),'Given fid is not valid');

[nr nc] = size(Matrix);
if nargin>3 ,f = format; else f = '%7.4f'; end
if nargin>4
    if use_siunitx == 1
        f = ['\\num{' f '}'];
    end
end
        

s1 = '\[';
s2 = [parname ' = \left( \begin{array}{' repmat('r',1,nc) '}'];
s3 = repmat([' ' f ' &'],1,nc); s3= ['  ' s3(1:end-2) '\\\\'];
s4 = '    \end{array} \right).';
s5 = '\]';

fprintf(fid,'%s\n',s1);
fprintf(fid,'%s\n',s2);
for a = 1:nr
    fprintf(fid,[s3 '\n'],Matrix(a,:));
end
fprintf(fid,'%s\n',s4);
fprintf(fid,'%s\n',s5);