function thisView = script_importPSIRasOverlay(thisView,Info,subjectInfo)
    %
    %   usage: importPSIRasOverlay(thisView,Info,subjectInfo)
    %      by: Ben Gurer
    %    date: 05/22/2018
    % purpose: cimport PSIR as an overlay and base anatomy
    %   input: mrView, Study information and subject information
    %  output: mrView with PSIR overlay and anotmay
    %
    
    % move to anatomy folder
    cd(fullfile(Info.dataDir,'Anatomy/originals/',subjectInfo.freeSurferName))
    % rename PSIR file so it works with matlab
    eval (['!cp ', [subjectInfo.freeSurferName, '_PSIR_pos_-.7_thr.nii'], ' ' [subjectInfo.freeSurferName, '_PSIR_mrTools.nii']])
    
    cd(fullfile(Info.dataDir,Info.studyDir,subjectInfo.subjectID));
    
    % import PSIR and PD
    thisView = viewSet(thisView,'newGroup','PSIR');
    thisView = viewSet(thisView,'curGroup','PSIR');
    thisView = importTSeries(thisView,[],'defaultParams=1',['pathname=' fullfile(Info.dataDir,'Anatomy/originals/',subjectInfo.freeSurferName,[subjectInfo.freeSurferName, '_PSIR_mrTools.nii'])]);
    thisView = newAnalysis(thisView,'dummy');
    thisView = importOverlay(thisView,[],'defaultParams=1',['pathname=' fullfile(Info.dataDir,'Anatomy/originals/',subjectInfo.freeSurferName,[subjectInfo.freeSurferName, '_PD_smooth7.nii'])]);
    thisView = importOverlay(thisView,[],'defaultParams=1',['pathname=' fullfile(Info.dataDir,'Anatomy/originals/',subjectInfo.freeSurferName,[subjectInfo.freeSurferName, '_PSIR_mrTools.nii'])]);
    thisView = viewSet(thisView,'overlaycolorrange',[.5 1]);
    %load PSIR as anatomy
    thisView = loadAnat(thisView,[subjectInfo.freeSurferName, '_PSIR_mrTools.nii'],fullfile(Info.dataDir,'Anatomy/originals/',subjectInfo.freeSurferName));
    
    return
    
    %%%%%%%%%%%%%%%%%%%%!!!!!!!!!!!!!!!!!!!!!!!!%%%%%%%%%%%%%%%%%%%
    
      %curvature and thickness flat maps
  thisView = viewSet(thisView,'curbase',viewGet(thisView,'basenum',[freeSurferName{iSubj} '_' sides{iSide} '_Flat_curvature']));
  saveAnat(thisView,'temp','OverWrite',0,fullfile(dataDir,studyDir,subjects{iSubj}))
  [thisView,params] = importOverlay(thisView,[],'defaultParams=1','justGetParams=1',['pathname=' fullfile(dataDir,studyDir,subjects{iSubj},'temp.nii')]);
  params.nameFrame1 = 'Curvature';
  [thisView,params] = importOverlay(thisView,params);
  thisView = viewSet(thisView,'curbase',viewGet(thisView,'basenum',[freeSurferName{iSubj} '_' sides{iSide} '_Flat_thickness']));
  saveAnat(thisView,'temp','OverWrite',0,fullfile(dataDir,studyDir,subjects{iSubj}))
  [thisView,params] = importOverlay(thisView,[],'defaultParams=1','justGetParams=1',['pathname=' fullfile(dataDir,studyDir,subjects{iSubj},'temp.nii')]);
  params.nameFrame1 = 'Thickness';
  [thisView,params] = importOverlay(thisView,params);
  thisView = viewSet(thisView,'curbase',viewGet(thisView,'basenum',[freeSurferName{iSubj} '_' sides{iSide} '_Flat_invFNIRT_' subjects{iSubj} 'Volume']));

 %regress curvature and thickness out of PSIR
  psirOverlay = viewGet(thisView,'overlayNum','PSIR (temp.nii)');
  fprintf('CORRELATION PSIR CURVATURE THICKNESS %s side\n',sides{iSide});
  thisView = regressOverlayOut(thisView,psirOverlay+[0 [4 3]+logical(B1field(iSubj))],viewGet(thisView,'curScan'));
  %average across depths
  curOverlay = viewGet(thisView,'curOverlay');
  [thisView,params] = combineTransformOverlays(thisView,[],'justGetParams=1','defaultParams=1',['overlayList=' num2str(curOverlay-1)]);
  params.combineFunction='averageDepthVol';
  [thisView,params] = combineTransformOverlays(thisView,params);
  thisView = viewSet(thisView,'overlayColorRange',[.5 1]);
  %smooth
  [thisView,params] = combineTransformOverlays(thisView,[],'justGetParams=1','defaultParams=1');
  params.combineFunction='spatialSmooth';
  params.additionalArgs = '[6 6 1]';
  [thisView,params] = combineTransformOverlays(thisView,params);
  thisView = viewSet(thisView,'overlayColorRange',[.5 1]);
  %compute gradient magnitude
  [thisView,params] = combineTransformOverlays(thisView,[],'justGetParams=1','defaultParams=1');
  params.combineFunction='imGradient2D';
  [thisView,params] = combineTransformOverlays(thisView,params);

end

%%% Flat volume
% for iSide=1:2
%   %first copy overlays to flat volume
%   thisView = viewSet(thisView,'curgroup',concatenationGroup);
%   thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',functionalAnalysis));
%   thisView = viewSet(thisView,'curbase',viewGet(thisView,'basenum',[freeSurferName{iSubj} '_' sides{iSide} '_Flat_invFNIRT_' subjects{iSubj}]));
%   [thisView,params] = combineTransformOverlays(thisView,[],'justGetParams=1','defaultParams=1',['overlayList=' mat2str(mainOverlays(iSubj,[8:10 10+iSide]))]);
%   params.combineFunction='User Defined';
%   params.customCombineFunction = 'plus'; %add 0
%   params.combinationMode = 'Apply function to each overlay';
%   params.additionalArgs = '0';
%   params.outputName=' ';
%   params.baseSpace = 1; %export result to base space (flat map)
%   params.exportToNewGroup=1; %export to volume in a new group
%   [thisView,params] = combineTransformOverlays(thisView,params);
%   curOverlay=viewGet(thisView,'curOverlay');
%   thisView = viewSet(thisView,'overlayColorRange',[0 5],curOverlay-[0 1 3]);
%   thisView = viewSet(thisView,'overlayColorRange',[0 8],curOverlay-2);
%   thisView = viewSet(thisView,'curSlice',6);
%   curOverlay=viewGet(thisView,'curOverlay');
%   thisView = viewSet(thisView,'overlaycolorrange',[0 8],curOverlay-2);
%   thisView = viewSet(thisView,'overlaycolorrange',[0 5],curOverlay-1);
%   thisView = viewSet(thisView,'overlaycolorrange',[0 5],curOverlay);
%   thisView = viewSet(thisView,'overlayMin',1,curOverlay-3);
%   thisView = viewSet(thisView,'alphaOverlay',curOverlay-3,curOverlay-(0:2));
%   thisView = viewSet(thisView,'alphaOverlayExponent',0,curOverlay-(0:2));
%  
%   %copy PSIR overlay
%   thisView = viewSet(thisView,'curgroup',psirGroup);
%   thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',psirAnalysis));
%   thisView = viewSet(thisView,'curOverlay',[1 2]);
%   thisView = viewSet(thisView,'curbase',viewGet(thisView,'basenum',[freeSurferName{iSubj} '_' sides{iSide} '_Flat']));
%   mrExport2SR(thisView.viewNum,fullfile(dataDir,studyDir,subjects{iSubj},'temp.nii'));
%   thisView = viewSet(thisView,'curgroup',[freeSurferName{iSubj} '_' sides{iSide} '_Flat_invFNIRT_' subjects{iSubj} 'Volume']);
%   thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum','combineTransformOverlays'));
%   [thisView,params] = importOverlay(thisView,[],'defaultParams=1','justGetParams=1',['pathname=' fullfile(dataDir,studyDir,subjects{iSubj},'temp.nii')]);
%   params.nameFrame1 = 'PSIR';
%   params.nameFrame2 = 'TFE PD';
%   [thisView,params] = importOverlay(thisView,params);
%   curOverlay=viewGet(thisView,'curOverlay');
%   thisView = viewSet(thisView,'overlayColorRange',[0.5 1],curOverlay-1);
%   
%   %copy MTR overlay
%   thisView = viewSet(thisView,'curgroup',mtrGroup);
%   thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',mtrAnalysis));
%   thisView = viewSet(thisView,'curOverlay',1);
%   thisView = viewSet(thisView,'curbase',viewGet(thisView,'basenum',[freeSurferName{iSubj} '_' sides{iSide} '_Flat']));
%   mrExport2SR(thisView.viewNum,fullfile(dataDir,studyDir,subjects{iSubj},'temp.nii'));
%   thisView = viewSet(thisView,'curgroup',[freeSurferName{iSubj} '_' sides{iSide} '_Flat_invFNIRT_' subjects{iSubj} 'Volume']);
%   thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum','combineTransformOverlays'));
%   [thisView,params] = importOverlay(thisView,[],'defaultParams=1','justGetParams=1',['pathname=' fullfile(dataDir,studyDir,subjects{iSubj},'temp.nii')]);
%   params.nameFrame1 = 'MTR';
%   [thisView,params] = importOverlay(thisView,params);
%   thisView = viewSet(thisView,'overlayColorRange',[0.2 .5]);
%   
%   %copy B1field overlay
%   if B1field(iSubj)
%     thisView = viewSet(thisView,'curgroup',b1Group);
%     thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',b1Analysis));
%     thisView = viewSet(thisView,'curOverlay',2);
%     thisView = viewSet(thisView,'curbase',viewGet(thisView,'basenum',[freeSurferName{iSubj} '_' sides{iSide} '_Flat']));
%     mrExport2SR(thisView.viewNum,fullfile(dataDir,studyDir,subjects{iSubj},'temp.nii'));
%     thisView = viewSet(thisView,'curgroup',[freeSurferName{iSubj} '_' sides{iSide} '_Flat_invFNIRT_' subjects{iSubj} 'Volume']);
%     thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum','combineTransformOverlays'));
%     [thisView,params] = importOverlay(thisView,[],'defaultParams=1','justGetParams=1',['pathname=' fullfile(dataDir,studyDir,subjects{iSubj},'temp.nii')]);
%     params.nameFrame1 = 'B1 error';
%     [thisView,params] = importOverlay(thisView,params);
%   end
%   
%   %curvature and thickness flat maps
%   thisView = viewSet(thisView,'curbase',viewGet(thisView,'basenum',[freeSurferName{iSubj} '_' sides{iSide} '_Flat_curvature']));
%   saveAnat(thisView,'temp','OverWrite',0,fullfile(dataDir,studyDir,subjects{iSubj}))
%   [thisView,params] = importOverlay(thisView,[],'defaultParams=1','justGetParams=1',['pathname=' fullfile(dataDir,studyDir,subjects{iSubj},'temp.nii')]);
%   params.nameFrame1 = 'Curvature';
%   [thisView,params] = importOverlay(thisView,params);
%   thisView = viewSet(thisView,'curbase',viewGet(thisView,'basenum',[freeSurferName{iSubj} '_' sides{iSide} '_Flat_thickness']));
%   saveAnat(thisView,'temp','OverWrite',0,fullfile(dataDir,studyDir,subjects{iSubj}))
%   [thisView,params] = importOverlay(thisView,[],'defaultParams=1','justGetParams=1',['pathname=' fullfile(dataDir,studyDir,subjects{iSubj},'temp.nii')]);
%   params.nameFrame1 = 'Thickness';
%   [thisView,params] = importOverlay(thisView,params);
%   thisView = viewSet(thisView,'curbase',viewGet(thisView,'basenum',[freeSurferName{iSubj} '_' sides{iSide} '_Flat_invFNIRT_' subjects{iSubj} 'Volume']));
%   
% end

% %DRAW INSULA EXCLUSION ROIS
% insulaROI{1} = 'leftInsulaFlatVol';
% insulaROI{2} = 'rightInsulaFlatVol';
% thisView = getMLRView;
% thisView = viewSet(thisView,'showRois','selected perimeter');
% 
% for iSide=1:2
%   thisView = viewSet(thisView,'curgroup',[freeSurferName{iSubj} '_' sides{iSide} '_Flat_invFNIRT_' subjects{iSubj} 'Volume']);
%   thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum','combineTransformOverlays'));
%   thisView = viewSet(thisView,'curbase',viewGet(thisView,'basenum',[freeSurferName{iSubj} '_' sides{iSide} '_Flat_invFNIRT_' subjects{iSubj} 'Volume']));
%   thisView = viewSet(thisView,'curROI',viewGet(thisView,'roinum',insulaROI{iSide}));
%   
%   %regress B1 error (and PD) out of MTR
%   mtrOverlay = viewGet(thisView,'overlayNum','MTR (temp.nii)');
%   thisView = viewSet(thisView,'overlayMax',.95,mtrOverlay);
%   if B1field(iSubj)
%     fprintf('CORRELATION MTR - B1error %s side\n',sides{iSide});
%     thisView = regressOverlayOut(thisView,mtrOverlay+[0 1],viewGet(thisView,'curScan'),[],[],[],[],0,0);
%     thisView = viewSet(thisView,'overlayColorRange',[.2 .6]);
%   end
%   fprintf('CORRELATION MTR - TFE-PD %s side\n',sides{iSide});
%   thisView = regressOverlayOut(thisView,mtrOverlay-[0 1],viewGet(thisView,'curScan'),[],[],[],[],0,2);
%   thisView = viewSet(thisView,'overlayColorRange',[.2 .6]);
% 
%   %regress curvature and thickness out of PSIR
%   psirOverlay = viewGet(thisView,'overlayNum','PSIR (temp.nii)');
%   fprintf('CORRELATION PSIR CURVATURE THICKNESS %s side\n',sides{iSide});
%   thisView = regressOverlayOut(thisView,psirOverlay+[0 [4 3]+logical(B1field(iSubj))],viewGet(thisView,'curScan'));
%   %average across depths
%   curOverlay = viewGet(thisView,'curOverlay');
%   [thisView,params] = combineTransformOverlays(thisView,[],'justGetParams=1','defaultParams=1',['overlayList=' num2str(curOverlay-1)]);
%   params.combineFunction='averageDepthVol';
%   [thisView,params] = combineTransformOverlays(thisView,params);
%   thisView = viewSet(thisView,'overlayColorRange',[.5 1]);
%   %smooth
%   [thisView,params] = combineTransformOverlays(thisView,[],'justGetParams=1','defaultParams=1');
%   params.combineFunction='spatialSmooth';
%   params.additionalArgs = '[6 6 1]';
%   [thisView,params] = combineTransformOverlays(thisView,params);
%   thisView = viewSet(thisView,'overlayColorRange',[.5 1]);
%   %compute gradient magnitude
%   [thisView,params] = combineTransformOverlays(thisView,[],'justGetParams=1','defaultParams=1');
%   params.combineFunction='imGradient2D';
%   [thisView,params] = combineTransformOverlays(thisView,params);
%   
%   if B1field(iSubj)
%     %regress curvature and thickness out of MTR
%     mtrOverlay = viewGet(thisView,'overlayNum','MTR (temp.nii) regress out B1 error (temp.nii), ');
%     fprintf('CORRELATION MTR(-B1error) CURVATURE THICKNESS %s side\n',sides{iSide});
%     thisView = regressOverlayOut(thisView,mtrOverlay-[0 1 2],viewGet(thisView,'curScan'));
%     %average across depths
%     curOverlay = viewGet(thisView,'curOverlay');
%     [thisView,params] = combineTransformOverlays(thisView,[],'justGetParams=1','defaultParams=1',['overlayList=' num2str(curOverlay-1)]);
%     params.combineFunction='averageDepthVol';
%     [thisView,params] = combineTransformOverlays(thisView,params);
%     thisView = viewSet(thisView,'overlayColorRange',[.2 .6]);
%     %smooth
%     [thisView,params] = combineTransformOverlays(thisView,[],'justGetParams=1','defaultParams=1');
%     params.combineFunction='spatialSmooth';
%     params.additionalArgs = '[6 6 1]';
%     [thisView,params] = combineTransformOverlays(thisView,params);
%     thisView = viewSet(thisView,'overlayColorRange',[.2 .6]);
%     %compute gradient magnitude
%     [thisView,params] = combineTransformOverlays(thisView,[],'justGetParams=1','defaultParams=1');
%     params.combineFunction='imGradient2D';
%     [thisView,params] = combineTransformOverlays(thisView,params);
%   end
% 
%   %regress curvature and thickness out of MTR
%   mtrOverlay = viewGet(thisView,'overlayNum','MTR (temp.nii) regress out TFE PD (temp.nii), ');
%   fprintf('CORRELATION MTR(-TFE-PD) CURVATURE THICKNESS %s side\n',sides{iSide});
%   thisView = regressOverlayOut(thisView,mtrOverlay-[0 1 2],viewGet(thisView,'curScan'));
%   %average across depths
%   curOverlay = viewGet(thisView,'curOverlay');
%   [thisView,params] = combineTransformOverlays(thisView,[],'justGetParams=1','defaultParams=1',['overlayList=' num2str(curOverlay-1)]);
%   params.combineFunction='averageDepthVol';
%   [thisView,params] = combineTransformOverlays(thisView,params);
%   thisView = viewSet(thisView,'overlayColorRange',[.2 .6]);
%   %smooth
%   [thisView,params] = combineTransformOverlays(thisView,[],'justGetParams=1','defaultParams=1');
%   params.combineFunction='spatialSmooth';
%   params.additionalArgs = '[6 6 1]';
%   [thisView,params] = combineTransformOverlays(thisView,params);
%   thisView = viewSet(thisView,'overlayColorRange',[.2 .6]);
%   %compute gradient magnitude
%   [thisView,params] = combineTransformOverlays(thisView,[],'justGetParams=1','defaultParams=1');
%   params.combineFunction='imGradient2D';
%   [thisView,params] = combineTransformOverlays(thisView,params);
% end
%   
% for iSide=1:2
%   % gradient reversals
%   thisView = viewSet(thisView,'curgroup',concatenationGroup);
%   thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',functionalAnalysis));
%   thisView = viewSet(thisView,'curbase',viewGet(thisView,'basenum',[freeSurferName{iSubj} '_' sides{iSide} '_Flat_invFNIRT_' subjects{iSubj}]));
%   [thisView,params] = combineTransformOverlays(thisView,[],'justGetParams=1','defaultParams=1',['overlayList=' mat2str(mainOverlays(iSubj,9))]);
%   params.combineFunction='gradientReversal';
%   params.additionalArgs = '[18 18 21]';
%   params.baseSpaceInterp = 'linear';
%   params.nOutputOverlays=7;
%   params.baseSpace = 1;
%   params.exportToNewGroup=1;
%   [thisView,params] = combineTransformOverlays(thisView,params);
%   curOverlay=viewGet(thisView,'curOverlay');
%   thisView = viewSet(thisView,'overlayMin',15,curOverlay-1);
%   thisView = viewSet(thisView,'overlayMax',180,curOverlay-1);
%   thisView = viewSet(thisView,'overlaycolorRange',[45 180],curOverlay-1);
%   thisView = viewSet(thisView,'overlayMax',75);
%   thisView = viewSet(thisView,'overlaycolorRange',[0 90]);
% end
