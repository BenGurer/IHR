function conNames = get_analysisConditionNames(thisView,analysisName,groupName,iScan)

thisView = viewSet(thisView,'curgroup',groupName);
if ~isempty(iScan)
    thisView = viewSet(thisView,'curScan', iScan);
end

thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',analysisName));
analysisData = viewGet(thisView,'analysis',viewGet(thisView,'analysisNum',analysisName));

conNames = analysisData.params.EVnames;

end