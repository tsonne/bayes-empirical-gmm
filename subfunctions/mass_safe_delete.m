% Delete all files given in cell (each string entry)
% 2016-12-05 tsonne: created
function N = mass_safe_delete(CellOfStrings)
N = 0;
% skip if empty
if isempty(CellOfStrings), return; end
K=numel(CellOfStrings);
for a=1:K
    x = CellOfStrings{a};
    if ischar(x) % only try if string in cell element
        y = safe_delete(x);
        N = N+y;
    end
end
end