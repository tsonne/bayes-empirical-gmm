function lib=get_linux_libgcc
% get most recent library directory of gcc, output as string.
% should work on different distros with different gcc versions.
% always tries to get 64bit lib, nothing else.
% will check for empty dirs (take most recent non-empty dir).
%
% tsonne 2015-12-18; created, works on openSUSE Tumbleweed 20151201
% tsonne 2017-03-09: no change, works on Ubuntu 16.06
p1 = '/usr/lib/gcc/'; % assuming 64bit default system [Debian]
p2 = '/usr/lib64/gcc/'; % other [openSUSE]
x1 = dir([p1 'x86_64*']);
x2 = dir([p2 'x86_64*']);
if ~isempty(x1)
    x = x1(1).name; p = p1;
elseif ~isempty(x2)
    x = x2(1).name; p = p2;
else
    fprintf('# ERROR: cannot find gcc lib directory on this system!\n');
    return
end
B = [p x '/'];
L = dir(B);
if isempty(L)
    fprintf('# ERROR: cannot find gcc lib directory on this system!\n');
    return
elseif numel(L)<3
    fprintf('# ERROR: gcc lib top directory %s empty!\n',B);
    return
else
    V = L(3:end); % get all versions, exclude . and ..
    % alphanumeric order means last index is most recent version
    NV = numel(V);
    c = 0; F = 0;
    while F ~= 1 && c ~= NV
        D = [B V(NV-c).name];
        ff = dir(D);
        if numel(ff) > 2 % more entries found than just . and .. (yes)
            lib = D; % this lib dir is populated with files, take it
            F = 1; % flag: got it
        else
            c = c + 1;
        end
    end
    if F ~= 1
        fprintf('# ERROR: gcc lib top dir %s : subdirs all empty!\n',B);
        return
    end
end

end