function t = adapt_tuning(t,acc_rate)
%
%    Tunes the scaling parameter for the proposal distribution
%    according to the acceptance rate over the last tune_interval:
%
%    Rate    Variance adaptation
%    ----    -------------------
%    <0.001        x 0.1
%    <0.05         x 0.5
%    <0.2          x 0.9
%    >0.5          x 1.1
%    >0.75         x 2
%    >0.95         x 10
%
% 2016-03-24 tsonne: from Hannes' python routine
%
%    # Switch statement
% if acc_rate <= 0.001
%     % reduce by 90 percent
%     t = t*0.1;
% elseif acc_rate < 0.05
%     % reduce by 50 percent
%     t = t*0.5;
% elseif acc_rate < 0.2
%     % reduce by ten percent
%     t = t*0.9;
% elseif acc_rate > 0.95
%     % increase by factor of ten
%     t = t*10;
% elseif acc_rate > 0.75
%     % increase by double
%     t = t*2;
% elseif acc_rate >= 0.5
%     % increase by ten percent
%     t = t*1.1;
% end
%
% Tim's extended version to achieve acceptance rates around 0.234:
% (acr=0.234 was recommended by Gelman et al. 1996, Roberts et al. 1997)
if acc_rate <= 0.001
    % reduce by 90 percent
    t = t*0.1;
elseif acc_rate < 0.05
    % reduce by 50 percent
    t = t*0.5;
elseif acc_rate < 0.1
    % reduce by 20 percent
    t = t*0.8;
elseif acc_rate < 0.2
    % reduce by five percent
    t = t*0.95;
elseif acc_rate > 0.95
    % increase by factor of ten
    t = t*10;
elseif acc_rate > 0.75
    % increase by double
    t = t*2;
elseif acc_rate > 0.5
    % increase by 20 percent
    t = t*1.2;
elseif acc_rate > 0.4
    % increase by ten percent
    t = t*1.1;
elseif acc_rate > 0.3
    % increase by five percent
    t = t*1.05;
end