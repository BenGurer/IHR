function overlayData = get_overlayData(thisView,overlayNames)
    %
    %   usage: get_overlayData(thisView,overlayNames)
    %      by: Ben Gurer
    %    date: 26/10/2017
    % purpose: Get overlay data from group
    %   input: mrTools view, names of overlays, 
    %          assumes thisView is the correct group and analysis
    %
if iscell(overlayNames)
    overlayStruct = cell(1,length(overlayNames));
    overlayData = struct;
    for i = 1:length(overlayNames)
        overlayStruct{i} = viewGet(thisView,'Overlay',overlayNames{i});
        overlayData.data{i} = overlayStruct{i}.data{1};
        overlayData.name{i} = overlayStruct{i}.name;
    end
else
        overlayStruct = viewGet(thisView,'Overlay',overlayNames);
        overlayData.data = overlayStruct.data;
        overlayData.name = overlayStruct.name;
end