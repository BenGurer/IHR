function data = get_analysisData(thisView,analysisName)
%
%   usage: get_analysisData(thisView)
%      by: Ben Gurer
%    date: 19/03/2018
% purpose: get analysis data from mrTools view
%   input: mrTools view, analysis name, ROI, scan number (if defined)
%

% get analysis data for specified group
thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',analysisName));
% get analysis data for specified analysis
data = viewGet(thisView,'analysis',viewGet(thisView,'analysisNum',analysisName));
end