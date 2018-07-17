function thisView = script_importFlatmaps(thisView,Info,subjectInfo)
    %
    %   usage: script_importFlatmaps(thisView,Info,subjectInfo)
    %      by: Ben Gurer
    %    date: 05/22/2018
    % purpose: import flatmaps
    %   input: mrView, Study information and subject information
    %  output: mrView with PSIR overlay and anotmay
    %
    
%import flat maps
params.path = fullfile(Info.dataDir,'Anatomy/freesurfer/subjects/',subjectInfo.freeSurferName,'/surfRelax');
for iSide=1:2
  params.anatFileName = [subjectInfo.freeSurferName '_mprage_pp.nii'];
  params.flatRes=3;
  params.threshold = 1;
  params.baseValues = 'curvature';
  params.flatFileName = [subjectInfo.flatmapNames{iSide} '.off' ];
  params.outerCoordsFileName = [subjectInfo.freeSurferName '_' Info.sides{iSide} '_GM.off'];
  params.innerCoordsFileName = [subjectInfo.freeSurferName '_' Info.sides{iSide} '_WM.off'];
  params.curvFileName = [subjectInfo.freeSurferName '_' Info.sides{iSide} '_Curv.vff'];
  base = importFlatOFF(params);
  thisView = viewSet(thisView, 'newbase', base);
  thisView = viewSet(thisView,'rotate',subjectInfo.flatmapRotation{iSide});
  thisView = viewSet(thisView,'corticalDepth',[0.3 0.7]);
end
refreshMLRDisplay(thisView);

end