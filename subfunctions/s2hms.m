function hms = s2hms(t)
% PURPOSE: Convert seconds to a vector of [hours, min, sec]
% 2019-01-15 Tim: created
hms = [fix(t/3600) fix((t-3600*fix(t/3600))/60) t-60*fix(t/60)];