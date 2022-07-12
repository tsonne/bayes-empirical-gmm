% http://stackoverflow.com/questions/11898105/safely-delete-file-in-matlab-
% without-risking-deletion-of-other-files-or-executi
% answered Aug 10 '12 at 16:18 amro
% 2016-12-05 tsonne: copied from stackoverflow to delete exactly ONE file
%                    (or none, and no dirs), regardless of input.
%                    Added number 'n' for checks.
function n = safe_delete(filename)
    n = 0; % number of deleted files
    % listing
    d = dir(filename);
    d([d.isdir]) = [];   % only files

    % skip if more than one match or no match
    if isempty(d) || numel(d) > 1, return; end

    % delete file
    p = fileparts(filename);
    delete( fullfile(p,d(1).name) );
    n = 1;
end