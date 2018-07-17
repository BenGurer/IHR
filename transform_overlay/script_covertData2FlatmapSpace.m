function [thisView, analysisData] = script_covertData2FlatmapSpace(thisView,groupName,analysisName,iScan,overlays,flatmapName)

% compare thisView current group to intended group


if viewGet(thisView,'curgroup') ~= viewGet(thisView,'groupNum',groupName)
    thisView = viewSet(thisView,'curgroup',groupName);
end
if ~isempty(iScan)
    %     thisView = viewSet(thisView,'curgroup',groupName,['curScan=' mat2str(iScan)]);
    %     thisView = viewSet(thisView,'curgroup',groupName,['scannum=' mat2str(iScan)]);
    thisView = viewSet(thisView,'curScan', iScan);
end

thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',analysisName));
analysisData = viewGet(thisView,'analysis',viewGet(thisView,'analysisNum',analysisName));
if isempty(overlays)
    overlays = 1:length(analysisData.overlays);
end

%% check if glm still works without below code
% data.conditionNames = analysisData.params.EVnames;
% for i = 1:length(analysisData.params.EVnames)
%     if ~isempty(iScan)
%         data.overlayConditionNames{i} = ['Scan ', mat2str(iScan), ' -  (' analysisData.params.EVnames{i} ',0)'];
%     else
%         data.overlayConditionNames{i} = [groupName ' (' analysisData.params.EVnames{i} ',0)'];
%     end
% end

% a = viewGet(thisView,'Overlay',data.overlayConditionNames{1})

% for iSide=1:2
% copy overlays to flat volume
thisView = viewSet(thisView,'curbase',viewGet(thisView,'basenum',flatmapName));
[thisView,params] = combineTransformOverlays(thisView,[],'justGetParams=1','defaultParams=1',['overlayList=' mat2str(overlays)]);
params.combineFunction='User Defined';
params.customCombineFunction = 'plus'; %add 0
params.combinationMode = 'Apply function to each overlay';
params.additionalArgs = '0';
if ~isempty(iScan)
    params.outputName = [analysisName ' '];
else
    params.outputName = [groupName '_' analysisName ' '];
end
params.baseSpace = 1; %export result to base space (flat map)
params.exportToNewGroup=1; %export to volume in a new group
[thisView,params] = combineTransformOverlays(thisView,params);
end