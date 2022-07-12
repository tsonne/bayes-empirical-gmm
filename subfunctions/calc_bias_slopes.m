% FUNCTION: calculate mean and stdev of bias and slopes of bias with
% distance and magnitude.
function [BiMS,RsMS,MsMS] = calc_bias_slopes(B,MRS)
% input: B = bias matrix, MRS = [Mw,dist,site] matrix
Np = size(B,2);
BiMS = [mean(B)' std(B)']; % MEAN and STDEV of bias
RsMS = nan(Np,2); % [Slope, std(Slope)] distance
MsMS = nan(Np,2); % [Slope, std(Slope)] magnitude
for a = 1:Np
    [~,slo,~,~,s_ab,~,~]= ts_lsq(log10(MRS(:,2)),B(:,a),0.05,'-q');
    RsMS(a,:) = [slo s_ab(1,2)];
    [~,slo,~,~,s_ab,~,~]= ts_lsq(    MRS(:,1)   ,B(:,a),0.05,'-q');
    MsMS(a,:) = [slo s_ab(1,2)];
end
end