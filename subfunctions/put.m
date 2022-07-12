function put(fid,varargin)
% write strings as they are to file, add newline on each.
% tsonne 2015-06-11
% tsonne 2016-04-01: allow without giving FID
if ~isnumeric(fid)
    fprintf('%s\n',fid);
    fprintf('%s\n',varargin{:});
else
    fprintf(fid,'%s\n',varargin{:});
end