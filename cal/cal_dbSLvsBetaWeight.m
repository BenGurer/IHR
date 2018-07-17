function [ ratio2Plot, level2Plot, fit , error2plot ] =  cal_dbSLvsBetaWeight(roi_av_ratio,stimulusLevels,baseLevel_dB)

ratio2Plot = [roi_av_ratio(stimulusLevels~=baseLevel_dB); mean(roi_av_ratio(stimulusLevels==baseLevel_dB))];

level2Plot = [stimulusLevels(stimulusLevels~=baseLevel_dB) baseLevel_dB];
% 
%% dB SL vs beta weight

fit = polyfit(level2Plot',ratio2Plot,1);

error2plot = [ mean(stimulusLevels(stimulusLevels==50)) , mean(roi_av_ratio(stimulusLevels==50)) , std(roi_av_ratio(stimulusLevels==50)) ];

end