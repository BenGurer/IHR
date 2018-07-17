
function [relativeDistances, overlayRoiData, pathDistances, verticesScanCoords] = cal_tonotopicMagnification(thisView,roiList,flatmapname,saveName)

%   usage: corticalMagnificationAuditory(thisView,saveName)
%      by: Ben Gurer based on Julien Besle "corticalMagnificationAuditory" based on "calcDist" by eli merriam, denis schluppeck, etc...
%    date: 18/06/2018
% purpose: estimates cortical magnification from best freqeuncy estiamtes and overlayData maps
%             by plotting the best freqeuncy value (corresponding to frequency that elicites the greatest response) 
%             as a function of distance from the low frequency gradient
%             reversal to the high frequency gradient reversal

%               - load a flat map 
%               -define one single-point rois representing the likely location of the foveal singularity
%               -define an ROI for the target cortical area (e.g. V1 or V2 dorsal)
%               -load in corAnal to display ph (average) / co (average) data
%               -then run corticalMagnificationAuditory().

%   input: mrView, pRF analysis information
%  output: updated mrView, pRF data
%

% corticalMagnificationAuditory.m
%
%        $Id$
%      usage: corticalMagnificationAuditory([view])
%         by: Julien Besle
%             based on calcDist by eli merriam, denis schluppeck, etc...
%
%       date: 2012-08-22
%    purpose: estimates cortical magnification from corAnal rings phase and overlayData maps
%             by plotting the phase value (corresponding to location in retinal/visual space) 
%             as a function of distance from the foveal singularity within a given cortical area
%             
%               - load a flat map 
%               -define one single-point rois representing the likely location of the foveal singularity
%               -define an ROI for the target cortical area (e.g. V1 or V2 dorsal)
%               -load in corAnal to display ph (average) / co (average) data
%               -then run corticalMagnificationAuditory().


% check arguments
% if ~any(nargin == [0 1 2])
%     help cal_tonotopicMagnification
%     return
% end

invFnirt=1;


% allow this to be run from command line
%  pick the first mrLoadRet view that exists
if ieNotDefined('thisView');
    vnums = viewGet([],'viewnums');
    if isempty(vnums)
        disp('uhoh - no views. just returning')
        return
    elseif numel(vnums) > 1
        disp('uhoh - more than one view open. check that everything is ok')
    end
    thisView = viewGet([],'view',vnums(1));
end

nRois = 3;
if isempty(roiList)
    roiList = viewGet(thisView,'curroi');
    %choose ROIs    
    while length(roiList)~=nRois
        roiList = selectInList(thisView,'roi','Select 2 reversal and 1 gradient ROIs');
        if isempty(roiList)
            disp('User pressed cancel');
            return;
        elseif  length(roiList)~=nRois
            mrWarnDlg('Please only select three ROIs')
        end
    end
    rois = viewGet(thisView,'roi',roiList);
else
    
    rois = viewGet(thisView,'roi',roiList);
end

%----------------------------------------------------------------------------------------------- Load surface coordinates
% baseCoords contains the mapping from pixels in the displayed slice to
%  voxels in the current base volume.
baseCoords = viewGet(thisView,'cursliceBaseCoords');
if isempty(baseCoords)
    mrWarnDlg('Load base anatomy before drawing an ROI');
    return;
end


% get the base CoordMap for the current flat patch
corticalDepth = mean(viewGet(thisView, 'corticalDepth'));
baseCoordMap = viewGet(thisView,'baseCoordMap');
if isempty(baseCoordMap)
    mrWarnDlg('(corticalMagnificationAuditory) You cannot use this function unless you are viewing a flat patch with a baseCoordMap');
    return;
end
baseHdr = viewGet(thisView, 'basehdr');
%load the flat patch 
flat = loadSurfOFF(fullfile(baseCoordMap.path, baseCoordMap.flatFileName));


% load the appropriate surface files
disp(sprintf('Loading %s', baseCoordMap.innerCoordsFileName));
surface.inner = loadSurfOFF(fullfile(baseCoordMap.path, baseCoordMap.innerCoordsFileName));
if isempty(surface.inner)
    mrWarnDlg(sprintf('(calcDist) Could not find surface file %s',fullfile(baseCoordMap.path, baseCoordMap.innerCoordsFileName)));
    return
end

surface.inner = xformSurfaceWorld2Array(surface.inner, baseHdr);

disp(sprintf('Loading %s', baseCoordMap.outerCoordsFileName));
surface.outer = loadSurfOFF(fullfile(baseCoordMap.path, baseCoordMap.outerCoordsFileName));
if isempty(surface.outer)
    mrWarnDlg(sprintf('(calcDist) Could not find surface file %s',fullfile(baseCoordMap.path, baseCoordMap.outerCoordsFileName)));
    return
end
surface.outer = xformSurfaceWorld2Array(surface.outer, baseHdr);
surface.filename = baseCoordMap.outerCoordsFileName;

% %restrict surface to flat patch
% flatInnerCoords = reshape(baseCoordMap.innerCoords,size(baseCoordMap.innerCoords,1)*size(baseCoordMap.innerCoords,2),3); %take flat patch coordinates (indices into the canonical base volume?)
% flatInnerCoords(any(flatInnerCoords==0,2),:)=[]; %remove any zeroes
% 
% %find vertices of the surface corresponding to the flat map
% nnz(ismember(round(flatInnerCoords),round(surface.inner.vtcs),'rows'))

% build up a mrMesh-style structure, taking account of the current corticalDepth
restrictToPatch=1;
if ~restrictToPatch
  m.vertices = surface.inner.vtcs+corticalDepth*(surface.outer.vtcs-surface.inner.vtcs);
  m.faceIndexList  = surface.inner.tris;
else
% %(restricting to flat patch)
  %re-order face vertices in patch to match that of the original surface
  %(this is only necessary for displaying the patch, as apparently, the
  %orientation of the faces matters and the order of vertices (for each face) 
  %differs between the OFF flat patch and the surface)
  patch2parent=flat.patch2parent(:,2);
  flatTrisSurf = patch2parent(flat.tris); %patch faces with surface vertices
  permutations = perms(1:3);
  flatTris = flat.tris;
  for i=1:size(permutations,1)
    whichTris = ismember(flatTrisSurf(:,permutations(i,:)),surface.inner.tris,'rows');
    flatTris(whichTris,:) = flatTris(whichTris,permutations(i,:));
  end
  m.vertices = surface.inner.vtcs(flat.patch2parent(:,2),:)+corticalDepth*(surface.outer.vtcs(flat.patch2parent(:,2),:)-surface.inner.vtcs(flat.patch2parent(:,2),:));
  m.faceIndexList  = flatTris;
end
[m.uniqueVertices,m.vertsToUnique,m.UniqueToVerts] = unique(m.vertices,'rows');


%find vertices at the current cortical depth that are in the ROIs
%We assume that both the vertices and the ROI coordinates are in the canonical base space
%at least make sure that this is the case for the ROIs
if any(any((viewGet(thisView,'base2roi',roiList(1))-eye(4))>eps)) ...
    || any(any((viewGet(thisView,'base2roi',roiList(2))-eye(4))>eps)) ...
    || any(any((viewGet(thisView,'base2roi',roiList(3))-eye(4))>eps))
  keyboard %need to convert roi coords to base space
end

%find all vertices in areal ROI
for i=1:nRois
  [~,m.roiVertices{i}] = intersect(round(m.uniqueVertices),round(rois(i).coords(1:3,:))','rows');
  m.nRoiVertices(i) = size(m.roiVertices{i},1);
end
%put gradient ROI last (assuming it is the largest)
[~,gradientRoi] = max(m.nRoiVertices);
m.roiVertices = m.roiVertices([setdiff(1:3,gradientRoi) gradientRoi]);
m.nRoiVertices = m.nRoiVertices([setdiff(1:3,gradientRoi) gradientRoi]);
%remove gradient ROI voxles that in the other 2 ROIs
m.roiVertices{3} = setdiff(m.roiVertices{3},union(m.roiVertices{1},m.roiVertices{2}));
m.nRoiVertices(3) = size(m.roiVertices{3},1);
rois = rois([setdiff(1:3,gradientRoi) gradientRoi]);


%----------------------------------------------------------------------------------------------- Find overlay values for areal ROI vertices
scanList = viewGet(thisView,'curscan');
if isempty(scanList)
  %select scans
  scanList = selectInList(thisView,'scans','Select scans to process');
end

%convert  vertices coords to scan coords (we assume that the current base is a surface of flat map)
surf2scan = viewGet(thisView,'base2scan');
for iRoi = 1:nRois
  verticesScanCoords{iRoi} = surf2scan*[m.uniqueVertices(m.roiVertices{iRoi},:) ones(m.nRoiVertices(iRoi),1)]';
  %interpolate values using chosen method
  cScan=0;
  interpMethod = 'nearest'; %TRY LINEAR
  interpMethod = 'linear';
  for iScan=scanList
    cScan=cScan+1;
    scanNames{cScan} = viewGet(thisView,'description',iScan);
    overlayData = viewGet(thisView,'overlaydata',iScan,viewGet(thisView,'curOverlay'));

    %note: we use interpn here instead of interp3 because that avoids having to swap X and Y coordinates
    overlayRoiData{cScan,iRoi} = interpn(overlayData,verticesScanCoords{iRoi}(1,:),verticesScanCoords{iRoi}(2,:),verticesScanCoords{iRoi}(3,:),interpMethod);
  end
end

%----------------------------------------------------------------------------------------------- Compute Dijkstra distance from fovea vertex to roi vertices
if invFnirt%load the non-fnirted surface
  %get the undistorted surface files
  if isempty(flatmapname)
      [filename,pathname] = uigetfile([baseCoordMap.path '/*.off'],'Select undistorted inner surface file');
      % load the appropriate surface files
      disp(sprintf('Loading %s', filename));
      surface.inner = loadSurfOFF(fullfile(pathname, filename));
      surface.inner = xformSurfaceWorld2Array(surface.inner, baseHdr); %we assume they're in the same whole-head anatomy space
      
      [filename,pathname] = uigetfile([pathname '/*.off'],'Select undistorted outer surface file');
      disp(sprintf('Loading %s', filename));
      surface.outer = loadSurfOFF(fullfile(pathname, filename));
      surface.outer = xformSurfaceWorld2Array(surface.outer, baseHdr);
      
      
      surface.filename = filename;
  else
      
      % load the appropriate surface files
      surface.inner = loadSurfOFF(fullfile(baseCoordMap.path, flatmapname.wm));
      surface.inner = xformSurfaceWorld2Array(surface.inner, baseHdr); %we assume they're in the same whole-head anatomy space
      
      surface.outer = loadSurfOFF(fullfile(baseCoordMap.path, flatmapname.gm));
      surface.outer = xformSurfaceWorld2Array(surface.outer, baseHdr);
      
      surface.filename = flatmapname.gm;
      
  end
  
  % replace the old vertices by the undistorted ones
  m.distortedVertices = m.vertices;
  m.vertices = surface.inner.vtcs(flat.patch2parent(:,2),:)+corticalDepth*(surface.outer.vtcs(flat.patch2parent(:,2),:)-surface.inner.vtcs(flat.patch2parent(:,2),:));
%   [uniqueVertices,vertsToUnique,UniqueToVerts] = unique(m.vertices,'rows'); %This gives different results, so use unique vector from distorted mesh
  m.uniqueVertices = m.vertices(m.vertsToUnique,:);
end

% calculate the connection matrix
m.uniqueFaceIndexList = findUniqueFaceIndexList(m);
m.connectionMatrix = findConnectionMatrix(m);

% find the distance of all vertices from their neighbours
%  D matrix contains all the distances from vertex i to j in one big sparse
%  matrix.
disp('(corticalMagnificationAuditory) calculating distances between neighbour nodes...')
scaleFactor = baseHdr.pixdim(2:4)'; % use the actual voxel dimensions of the base instead of assuming 1 mm isotropic
D = find3DNeighbourDists(m,scaleFactor);

% using the D matrix for the graph calculated above, find path between each
% reversal vertex and _all_ other vertices in the mesh
% use dijkstrap to also spit out the predecessor matrix that will allow use
% to define the actual paths.
disp('(corticalMagnificationAuditory) Getting shortest paths between ROIs w/ dijkstra ...')
for i=1:2
  [m.dist{i}, m.pred{i}] = dijkstrap( D, m.roiVertices{i}  );
  m.pred{i}(~isinf(m.dist{i})&m.pred{i} == -1) = 0; % set self-distances to 0 (JB: not sure if m.dist shouldn't be used as an input to pred2path instead)
  pathDistances{i} = m.dist{i}(:,m.roiVertices{3});   % distance to all vertices of gradient ROI

  %get shortest distance
  [pathDistances{i}, whichEnd{i}] = min(pathDistances{i});
  
  %this is needed only to get the details of the paths and can be sipped if
  %only the distances are needed
  % rte will be a cell array numel(src) by numel(trg); each entry in the cell
  % array contain a vector of vertices that define the path from the
  % corresponding vertex in src to trg
  m.rte{i} = pred2path(m.pred{i}, m.roiVertices{i}, m.roiVertices{3}); %JB: note that this requires a modified version of pre2path that DOES NOT reorder the sources
end

pathDistances = [pathDistances{1};pathDistances{2}];
whichEnd = [whichEnd{1};whichEnd{2}];

%calculate relative distance from each vertex in the gradient ROI to
%closest points in the 2 reversal ROIs
relativeDistances = pathDistances./repmat(sum(pathDistances),2,1);

[~,lowFrequencyBorder] = min([nanmean(overlayRoiData{1}) nanmean(overlayRoiData{2})]);

%(average distance estimate for unique overlay voxels?)

p=polyfit(relativeDistances(lowFrequencyBorder,:),overlayRoiData{3},2);
pf = @(x)p(1)*x.^2 + p(2)*x + p(3);

if ~ieNotDefined('saveName')
  save(saveName,'relativeDistances','overlayRoiData');
end

figure;
if ~ieNotDefined('saveName')
  set(gcf,'name',saveName);
end
plot(relativeDistances(lowFrequencyBorder,:),overlayRoiData{3},'.');
hold on
plot(0:0.01:1,pf(0:0.01:1),'r');
plot(zeros(size(overlayRoiData{lowFrequencyBorder})),overlayRoiData{lowFrequencyBorder},'.k');
plot(ones(size(overlayRoiData{3-lowFrequencyBorder})),overlayRoiData{3-lowFrequencyBorder},'.k');

% figure;
% if ~ieNotDefined('saveName')
%   set(gcf,'name',saveName);
% end
% semilogy(relativeDistances(lowFrequencyBorder,:),nErb2f(overlayRoiData{3},0.25,6,7),'ok');
% hold on
% semilogy(0:0.01:1,nErb2f(pf(0:0.01:1),0.25,6,7),'m','linewidth',2);
% semilogy(zeros(size(overlayRoiData{lowFrequencyBorder})),nErb2f(overlayRoiData{lowFrequencyBorder},0.25,6,7),'ob');
% semilogy(ones(size(overlayRoiData{3-lowFrequencyBorder})),nErb2f(overlayRoiData{3-lowFrequencyBorder},0.25,6,7),'or');
% xlim([-0.05 1.05]);


% 
% keyboard
% return

% %---------------------------------------------define shortest paths between border ROIs
% 
% %find the shortest paths that involve one of the two border ROIs (this is
% %equivalent to an 'OR' version of unique applied to the distance-sorted
% %end-start pairs)
% 
% %first create end-start pairs and sort them by distance
% pathStarts = repmat(m.roiVertices{1},1,m.nRoiVertices(2));
% pathEnds = repmat(m.roiVertices{2}',m.nRoiVertices(1),1);
% 
% pathDistances = pathDistances(:);
% pathStarts = pathStarts(:);
% pathEnds = pathEnds(:);
% 
% [sortedDistances,pathSortingIndices] = sort(pathDistances);
% sortedStarts=pathStarts(pathSortingIndices);
% sortedEnds = pathEnds(pathSortingIndices);
% sortedEndpoints = [sortedStarts sortedEnds];
% 
% %then go through all the distance-sorted start-end pairs and select the
% %first instance of each start and end point until all points have been found
% foundStarts = [];
% foundEnds = [];
% toKeep = [];
% c=0;
% while ~isempty(setdiff(m.roiVertices{1},foundStarts)) || ~isempty(setdiff(m.roiVertices{2},foundEnds))
%   c=c+1;
%   if ~ismember(sortedStarts(c), foundStarts)
%     toKeep = [toKeep;c];
%     foundStarts = [foundStarts; sortedStarts(c)];
%   end
%   if ~ismember(sortedEnds(c), foundEnds)
%     toKeep = [toKeep;c];
%     foundEnds = [foundEnds; sortedEnds(c)];
%   end
% end
% toKeep = unique(toKeep);
% 
% pathDistances = sortedDistances(toKeep);
% pathStarts = sortedStarts(toKeep);
% pathEnds = sortedEnds(toKeep);
% paths = m.rte(pathSortingIndices(toKeep));
%   
% 
% 
% m.vertices(:,2) = -1*m.vertices(:,2); %change orientation for display
% 
% % ROI patch
% % roiVertices = m.vertsToUnique(m.roiVertices{3});
% % roiFaces = m.faceIndexList((all(ismember(m.faceIndexList,roiVertices),2)),:);
% 
% roiVertices1 = m.vertices(m.vertsToUnique(m.roiVertices{1}),:);
% roiVertices2 = m.vertices(m.vertsToUnique(m.roiVertices{2}),:);
% roiVertices3 = m.vertices(m.vertsToUnique(m.roiVertices{3}),:);
% myMrPrintSurf(m.vertices,m.faceIndexList);
% % patch('vertices', m.vertices, 'faces',roiFaces,'FaceVertexCData', .6*ones(size(m.vertices,1),3),'facecolor','interp','edgecolor',[0 0 0]);
% hold on
% % plot3(roiVertices1(:,1),roiVertices1(:,2),roiVertices1(:,3),'.','color',color2RGB(rois(1).color));
% % plot3(roiVertices2(:,1),roiVertices2(:,2),roiVertices2(:,3),'.','color',color2RGB(rois(2).color));
% plot3(roiVertices1(:,1),roiVertices1(:,2),roiVertices1(:,3),'.r');
% plot3(roiVertices2(:,1),roiVertices2(:,2),roiVertices2(:,3),'.b');
% % plot3(roiVertices3(:,1),roiVertices3(:,2),roiVertices3(:,3),'.r');
% 
% point = find(abs(relativeDistances(1,:)-1/3)<0.0005);
% thisPoint = m.vertices(m.vertsToUnique(m.roiVertices{3}(point)),:);
% thisPath = m.vertices(m.vertsToUnique(m.rte{1}{whichEnd(1,point),point}),:);
% plot3(thisPath(:,1),thisPath(:,2),thisPath(:,3),'m','linewidth',2);
% plot3(thisPoint(:,1),thisPoint(:,2),thisPoint(:,3),'.g');
% thisPath = m.vertices(m.vertsToUnique(m.rte{2}{whichEnd(2,point),point}),:);
% plot3(thisPath(:,1),thisPath(:,2),thisPath(:,3),'c','linewidth',2);
% 
% 
% % nColors=length(paths);
% % cmap = hsv(nColors);
% % for i=1:length(paths)
% %   thisPath = m.vertices(m.vertsToUnique(paths{i}),:);
% %   plot3(thisPath(:,1),thisPath(:,2),thisPath(:,3),'color',cmap(rem(i,nColors)+1,:),'linewidth',2);
% % end
% % 
% keyboard,
% return

%----------------------------------------------------------------------------------------------- plot data and fit cortical magnification parameters

% 
% coherenceThreshold = 0.5;
% maxEccentricity = 5.5;
% minEccentricity = 0;
% twoThirdDistance = (max(m.pathDistances)+min(m.pathDistances))*.6;
% minEccDistantVertices = 2;
% cScan=0;
% figure('name',surface.filename);
% numPlotsPerScan=3;
% for iScan=scanList
%   cScan=cScan+1;
%   
% 
%   fprintf(1,'converting phase to eccentricity with eccentricityMode=''%s'', minRadius = %.2f, maxRadius=%.2f, dutyCycle=%.2f\n', eccentricityMode, minRadius, maxRadius, dutyCycle);
%   eccentricity = phase2eccentricity(m.phase{cScan}, eccentricityMode, minRadius, maxRadius, dutyCycle);
%   %threshold the data
%   aboveThreshold = overlayRoiData{cScan}>=coherenceThreshold;
%   eccThres = double(eccentricity(aboveThreshold));
%   distThres = double(m.pathDistances(aboveThreshold));
%   % stop at max eccentricity for fit
%   distThres = distThres(eccThres<maxEccentricity & eccThres>minEccentricity);
%   eccThres = eccThres(eccThres<maxEccentricity & eccThres>minEccentricity);
%   %also exclude voxels with low phase values that are unreasonably far from the fovea
%   %for example, for now, all eccentricites less than 1 in the furthest third of remaining vertices
%   toKeep = distThres<twoThirdDistance | eccThres >minEccDistantVertices;
%   distData{cScan} = distThres(toKeep);
%   eccData{cScan} = eccThres(toKeep);
% 
%   %fit logarithmic function with a parameters for x and y offsets
%   h=subplot(length(scanList),numPlotsPerScan,(cScan-1)*numPlotsPerScan+1);
%   scatter(eccentricity,m.pathDistances,30,[.7 .7 .7]);
%   hold on
%   scatter(eccentricity(aboveThreshold),m.pathDistances(aboveThreshold),30,'g');
%   scatter(eccThres,distThres,30)
%   scatter(eccData{cScan},distData{cScan},30,'b')
%   ylabel('Cortical distance from ''fovea''');
%   xlabel('Eccentricity')
%   [sortedEccData,sortingVector] = sort(eccData{cScan});
%   sortedDistData = distData{cScan}(sortingVector);
%   
%   func1 = @(params,sortedEccData)params(1)+params(2).*log(sortedEccData+params(3));
% %   lowerbound = [-5 5 0];
% %   upperbound = [50 30 10];
% %   params1 = lsqcurvefit(func1,[0 1 0],sortedEccData,sortedDistData,lowerbound,upperbound);
%    params1 = lsqcurvefit(func1,[0 1 0],sortedEccData,sortedDistData);
%   %compute derivative values at different eccentricities
%   eccentricities=[1 3 5]';
%   fprintf('\n\tEcc\t Dist (With b parameter)\n');
%   distances = params1(2)./(eccentricities+params1(3));
%   [eccentricities distances]
%   
%   plot(sortedEccData,func1(params1,sortedEccData),'k','linewidth',2);
% %   legend({sprintf('Data points (max eccentricty = %.1f deg)',maxEccentricity),sprintf('Fit: D = %.2f + %.2f * log(E + %.2f)',params1(1),params1(2),params1(3) )},'location','SouthEast');
%   title('D = b + k * log(E + a)');
% 
%   %without b parameter
%   func2 = @(params,sortedEccData)params(1).*log(sortedEccData+params(2));
%   params2 = lsqcurvefit(func2,[1 0],sortedEccData,sortedDistData);
%   %compute derivative values at different eccentricities
%   eccentricities=[1 3 5]';
%   fprintf('\n\tEcc\t Dist (Without b parameter)\n');
%   distances = params2(1)./(eccentricities+params2(2));
%   [eccentricities distances]
%   
% 
%   plot(sortedEccData,func2(params2,sortedEccData),'r','linewidth',2);
%   set(h,'Ylim',[0 50]);
%   set(h,'Xlim',[0 12]);
%   legend({sprintf('Coherence > %.1f',coherenceThreshold),...
%           sprintf('Eccentricty > %.1f or < %.1f deg',maxEccentricity,minEccentricity),...
%           sprintf('Distance > %.2f & Eccentricity < %.1f',twoThirdDistance,minEccDistantVertices),...
%           'Data points used for fit',...
%           sprintf('Fit: D = %.2f + %.2f * log(E + %.2f)',params1(1),params1(2),params1(3) ),...
%           sprintf('Fit: D = %.2f * log(E + %.2f)',params2(1),params2(2) )},...
%           'location','SouthEast');
%  
%   p=get(h,'position');
%   uicontrol('style','text','unit','normalized','String',scanNames{cScan},'position',[0.02 p(2)+p(4)/3 p(1)-0.05 p(4)/3]);
%   
%   %%%%%%%%%%%%%%% Fit exponential like in Larsson & Heeger 2006 J. Neurosci.
% 
%   %set cortical distance to 0 at 3 deg eccentricity
%   %first estimate cortical distance at 3 deg eccentricity using previous fit
%   dist3deg = func1(params1,3);
% 
%   [sortedDistData,sortingVector] = sort(distData{cScan}-dist3deg);
%   sortedEccData = eccData{cScan}(sortingVector);
% 
%   func = @(params,sortedDistData) exp(params(1).*(sortedDistData+params(2)));
%   params = lsqcurvefit(func,[0 0],sortedDistData,sortedEccData);
% 
%   % %this is the same using lsqnonlin
%   % func = @(params) sortedEccData - exp(params(1).*(sortedDistData+params(2)));
%   % params = lsqnonlin(func,[0 0]);
% 
%   subplot(length(scanList),numPlotsPerScan,(cScan-1)*numPlotsPerScan+2);
%   scatter(m.pathDistances-dist3deg,eccentricity,30,[.7 .7 .7]);
%   hold on
%   scatter(m.pathDistances(aboveThreshold)-dist3deg,eccentricity(aboveThreshold),30,'g');
%   scatter(distThres-dist3deg,eccThres,30)
%   scatter(distData{cScan}-dist3deg,eccData{cScan},30,'b')
%   xlabel(sprintf('Cortical distance from 3 deg (%.0f mm)',dist3deg));
%   ylabel('Eccentricity')
%   plot(sortedDistData,func(params,sortedDistData),'k','linewidth',2);
%   legend({sprintf('Coherence > %.1f',coherenceThreshold),...
%           sprintf('Eccentricty > %.1f or < %.1f deg',maxEccentricity,minEccentricity),...
%           sprintf('Distance > %.2f & Eccentricity < %.1f',twoThirdDistance,minEccDistantVertices),...
%           'Data points used for fit',...
%           sprintf('Fit: E = exp[%.2f * (D + %.2f)]',params(1),params(2))},...
%           'location','NorthWest');
%   title('E = exp[c*(D + d)] (Larsson & Heeger 2006)');
% 
% 
% 
%   %%%%%%%%%%%%%%% Derive and fit inverse function
%   [sortedEccData,sortingVector] = sort(eccData{cScan});
%   sortedDistData = distData{cScan}(sortingVector);
% 
%   %bin values by 10 
%   numPerBins=20;
%   numBins = floor(length(sortedEccData)/numPerBins);
%   sortedEccData = reshape(sortedEccData(1:numBins*numPerBins),numPerBins,numBins);
%   sortedDistData = reshape(sortedDistData(1:numBins*numPerBins),numPerBins,numBins);
%   meanSortedEccData = mean(sortedEccData);
%   meanSortedDistData = mean(sortedDistData);
%   %numerically derive
%   derivDist = diff(meanSortedDistData)./diff(meanSortedEccData);
%   centreBinEcc = (meanSortedEccData(1:end-1)+meanSortedEccData(2:end))/2;
%   %fit inverse function
%   func = @(params,eccentricity) params(1)./(eccentricity+params(2));
%   params = lsqcurvefit(func,[1 0],centreBinEcc,derivDist);
% 
%   subplot(length(scanList),numPlotsPerScan,(cScan-1)*numPlotsPerScan+3);
%   plot(centreBinEcc,derivDist,'o');
%   ylabel('Cortical distance');
%   xlabel('Eccentricity')
%   hold on
%   plot(centreBinEcc,func(params,centreBinEcc),'k','linewidth',2);
%   legend({sprintf('Numerical derivative (bins of %d points)',numPerBins),sprintf('Fit: D = %.2f/(E + %.2f)',params(1),params(2) )},'location','NorthEast');
%   title('D'' = k / (E + a)');
% 
% end
% 
% 
% % fit exp function to eccentricity = f(cortical distance)
% 
% % fit 1/(x+a) to deriv(cortical distance) = f(eccentricity)
% % for this need to bin the data somehow
% 

%%%%%%%%%%%%%%%% Modified pred2path function %%%%%%%%%%%%%%%%%%%%
% does not reorder the source vertices in the output


function rte = pred2path(P,s,t)
%PRED2PATH Convert predecessor indices to shortest paths from 's' to 't'.
%   rte = pred2path(P,s,t)
%     P = |s| x n matrix of predecessor indices (from DIJK)
%     s = FROM node indices
%       = [] (default), paths from all nodes
%     t = TO node indices
%       = [] (default), paths to all nodes
%   rte = |s| x |t| cell array of paths (or routes) from 's' to 't', where
%         rte{i,j} = path from s(i) to t(j)
%                  = [], if no path exists from s(i) to t(j)
%
% (Used with output of DIJK)

% Copyright (c) 1994-2006 by Michael G. Kay
% Matlog Version 9 13-Jan-2006 (http://www.ie.ncsu.edu/kay/matlog)

% Input Error Checking ****************************************************
error(nargchk(1,3,nargin));

[rP,n] = size(P);

if nargin < 2 || isempty(s), s = (1:n)'; else s = s(:); end
if nargin < 3 || isempty(t), t = (1:n)'; else t = t(:); end

if any(P < 0 | P > n)
   error(['Elements of P must be integers between 1 and ',num2str(n)]);
elseif any(s < 1 | s > n)
   error(['"s" must be an integer between 1 and ',num2str(n)]);
elseif any(t < 1 | t > n)
   error(['"t" must be an integer between 1 and ',num2str(n)]);
end
% End (Input Error Checking) **********************************************

rte = cell(length(s),length(t));

[idx_i,idxs] = find(P==0);
idxs(idx_i) = idxs;   %JB: Matlab returns i and j indices sorted along the second dimension (j)
                      %    in order to keep the ordering of the input, need to re-order j indices by i indices

for i = 1:length(s)
%    if rP == 1
%       si = 1;
%    else
%       si = s(i);
%       if si < 1 | si > rP
%          error('Invalid P matrix.')
%       end
%    end
   si = find(idxs == s(i));
   for j = 1:length(t)
      tj = t(j);
      if tj == s(i)
         r = tj;
      elseif P(si,tj) == 0
         r = [];
      else
         r = tj;
         while tj ~= 0
            if tj < 1 || tj > n
               error('Invalid element of P matrix found.')
            end
            r = [P(si,tj) r];
            tj = P(si,tj);
         end
         r(1) = [];
      end
      rte{i,j} = r;
   end
end

if length(s) == 1 && length(t) == 1
   rte = rte{:};
end

%rte = t;
while 0%t ~= s
   if t < 1 || t > n || round(t) ~= t
      error('Invalid "pred" element found prior to reaching "s"');
   end
   rte = [P(t) rte];
   t = P(t);
end

% ********** lcfNErb **********
function nerb = lcfNErb(f)
  nerb = 21.4*log10(4.37*f+1);

% ***** lcfInvNErb *****
function f = lcfInvNErb(nerb)
  f = 1/4.37*(10.^(nerb/21.4)-1);
  

function f = nErb2f(nerb,lowFreq,highFreq,nFreqs) 

f = lcfInvNErb( lcfNErb(lowFreq) + (nerb-1) / (nFreqs-1) * (lcfNErb(highFreq)-lcfNErb(lowFreq)));




