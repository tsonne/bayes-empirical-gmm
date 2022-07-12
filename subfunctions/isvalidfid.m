function valid=isvalidfid(fid)
% check whether given file identifier fid is an open file or not
% tsonne 2015-10-13
if fid>2
    if ftell(fid)~=-1
        valid=true;
    else
        valid=false;
    end
else
    valid=false;
end