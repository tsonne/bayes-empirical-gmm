function v = getpar(S,field,default)
% Get parameter value from struct or default.
% INPUT:
%   S       = Struct
%   field   = Field to be extracted from the S struct
%   default = Default value if field is not member of S.
%             If no default, then field must be given (non-optional).
% OUTPUT:
%   v  = field value or default
if isfield(S,field)
    v = S.(field);
elseif nargin>2
    v = default;
else
    error('Need value for struct field: %s',field);
end
end