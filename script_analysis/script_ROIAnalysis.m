function roiData = script_ROIAnalysis(roiData,Info,glmInfo,stimInfo,plotInfo,subjectInfo,analysisScanNum,dataType,ROInamesSelect)

%% perform ROI analysis on GLM data

%% calculate
for iAnal = 1:length(glmInfo.analysisBaseNames_Scans)/glmInfo.nScans
    %% don't want to loop over analysis names - want to just do it once for each analysis
    eval(['roiData.splitData.' glmInfo.analysisBaseNames_Scans{iAnal}  '.roiAnalysis = cal_splitAverage_roi_GLM(roiData.scanData,glmInfo.analysisBaseNames_Scans{iAnal},subjectInfo.conditionOrder,[],Info.ConATrue,dataType,stimInfo);'])
end

for iAnal = 1:length(glmInfo.analysisNames_Groups)/length(glmInfo.groupNames)
%     data = cal_average_roi_GLM(roiData, roiName, groupNames, analName,restrictIndex,conATrue,dataType,stimInfo)
    eval(['roiData.concatData.' glmInfo.analysisNames_Groups{iAnal} '.roiAnalysis = cal_average_roi_GLM(roiData,glmInfo.groupNames,glmInfo.analysisNames_Groups{iAnal},[],Info.ConATrue,dataType,stimInfo);'])
end

%% plot
for iAnal = 1:length(glmInfo.analysisBaseNames_Scans)/glmInfo.nScans
    %     plot_average_roi_GLM(data,analName,stimInfo,plotInfo,figureName)
    figureName = ['ROI: ' roiData.roi.name, ' Anal: ', glmInfo.analysisBaseNames_Scans{iAnal}];
    plot_splitAverage_roi_GLM(eval(['roiData.splitData.'  glmInfo.analysisBaseNames_Scans{iAnal} '.roiAnalysis']),  glmInfo.analysisBaseNames_Scans{iAnal} ,stimInfo,plotInfo,figureName);
end
for iAnal = 1:length(glmInfo.analysisNames_Groups)/length(glmInfo.groupNames)
%     plot_average_roi_GLM(data,analName,stimInfo,plotInfo,figureName)
    figureName = ['ROI: ' roiData.roi.name, ' Anal: ', glmInfo.analysisBaseNames_Scans{iAnal}];
    plot_average_roi_GLM(eval(['roiData.concatData.' glmInfo.analysisNames_Groups{iAnal} '.roiAnalysis']), glmInfo.analysisNames_Groups{iAnal},stimInfo,plotInfo,figureName);
end



%% need to add concat level ROI analysis


% %% Plot ROI analysis
%     for iPlot = 1:length(plotInfo.ROIplotList)
%         plot_splitAverage_roi_GLM(eval(['roiData.scanData.' ROInames{iROI}]),plotInfo.ROIplotList{iPlot},stimInfo,plotInfo,ROInames{iROI});
%     end

% case 'Groups'
    




% 
% %% Split GLM ROI analysis
% % Perform ROI analysis
% if isempty(ROInamesSelect)
%     ROInames = Info.ROInames;
% else
%     if iscell(ROInamesSelect)
%         ROInames = ROInamesSelect;
%     else
%         ROInames = {ROInamesSelect};
%     end
% end
% 
% % switch GroupType
% %     case 'Scans'
% %% perform ROI analysis on GLM data
% for iROI = 1:length(ROInames)
%     for iAnal = 1:length(glmInfo.analysisBaseNames_Scans)
%         %% don't want to loop over analysis names - want to just do it once for each analysis
%         eval(['roiData.groupData.' ROInames{iROI} '.roiAnalysis_' glmInfo.analysisBaseNames_Scans{iAnal}  ' = cal_splitAverage_roi_GLM(roiData.scanData.' ROInames{iROI} ',glmInfo.analysisBaseNames_Scans{iAnal},subjectInfo.conditionOrder,[],Info.ConATrue,dataType,stimInfo);'])
%     end
% end
% %% calculate ratio between dB senstation level and BOLD fMRI activity
% 
% %% need to add concat level ROI analysis
% 
% 
% %% Plot ROI analysis
% for iROI = 1:length(ROInames)
%     for iPlot = 1:length(plotInfo.ROIplotList)
%         plot_splitAverage_roi_GLM(eval(['roiData.scanData.' ROInames{iROI}]),plotInfo.ROIplotList{iPlot},stimInfo,plotInfo,ROInames{iROI});
%     end
% end
% % case 'Groups'
%     
% for iROI = 1:length(ROInames)
%     for iAnal = 1:length(glmInfo.analysisNames_Groups)
%         eval(['roiData.groupData.' ROInames{iROI} '.roiAnalysis_' glmInfo.analysisNames_Groups{iAnal} '= cal_average_roi_GLM(roiData,ROInames{iROI},glmInfo.groupNames,glmInfo.analysisNames_Groups{iAnal},[],Info.ConATrue,dataType,stimInfo);'])
%     end
% end
% 
%  for iROI = 1:length(ROInames)
%     for iAnal = 1:length(glmInfo.analysisNames_Groups)/length(glmInfo.groupNames)
%         plot_average_roi_GLM(eval(['roiData.groupData.' ROInames{iROI}]),['roiAnalysis_' glmInfo.analysisNames_Groups{iAnal+1}],stimInfo,plotInfo,ROInames{iROI});
%     end
% end   
end

