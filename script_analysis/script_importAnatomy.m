function thisView = script_importAnatomy(thisView,Info,subjectInfo)


%% load reference EPI as anatomy
thisView = loadAnat(thisView,'lastFrameEPI.nii',fullfile(Info.dataDir,Info.studyDir,subjectInfo.subjectID,'FNIRT'));

%% Load HighRes inplane T2*
thisView = loadAnat(thisView,[subjectInfo.niftiBaseName, subjectInfo.T2star, '_1_modulus_stripped.nii'],fullfile(Info.dataDir,Info.studyDir,subjectInfo.subjectID,'Anatomy'));

%load skull-stripped  EPI as overlay
% [thisView,params] = importOverlay(thisView,[],'defaultParams=1',['pathname=' fullfile(dataDir,studyDir,subjects{iSubj},'/FNIRT/lastFrameEPI_stripped.nii')]);

% script running glm on each scan individually and performing split
% analysis

%IMPORT  FREESURFER SURFACES
if exist(fullfile(Info.dataDir,'Anatomy/freesurfer/subjects/',subjectInfo.freeSurferName,'/surfRelax')) ~= 7
cd(fullfile(Info.dataDir,'Anatomy/freesurfer/subjects/',subjectInfo.freeSurferName));
% Below needs scripting to only be performed if files not made
% % If 7T
% mlrImportFreeSurfer('defaultParams=1','volumeCropSize=[336 336 336]'); % this will depend on the resolution of the PSIR
% elseif 3T
% % mlrImportFreeSurfer('defaultParams=1','volumeCropSize=[240 240 175]');
% end
% 

%apply FNIRT warping coefficient to surfaces
cd(fullfile(Info.dataDir,Info.studyDir,subjectInfo.subjectID));
fslApplyWarpSurfOFF(fullfile(Info.dataDir,Info.studyDir,subjectInfo.subjectID,'FNIRT/',[subjectInfo.niftiBaseName subjectInfo.T2star '_1_modulus_crop_resampled_2_lastFrameEPI_warpcoef.nii']),...
    fullfile(Info.dataDir,Info.studyDir,subjectInfo.subjectID,'FNIRT/','lastFrameEPI.nii'),...
    fullfile(Info.dataDir,'Anatomy/freesurfer/subjects/',subjectInfo.freeSurferName, 'surfRelax',[subjectInfo.freeSurferName '_mprage_pp.nii']),...
    subjectInfo.subjectID);

end
%import surfaces (this step has not been made scriptable yet) 

for iSide=1:2
  base = importSurfaceOFF(fullfile(Info.dataDir,'Anatomy/freesurfer/subjects/',subjectInfo.freeSurferName,'surfRelax',...
   [subjectInfo.freeSurferName '_' Info.sides{iSide} '_GM.off']));
  thisView = viewSet(thisView, 'newbase', base);
  thisView = viewSet(thisView,'corticalDepth',[0.3 0.7]);
end

refreshMLRDisplay(thisView.viewNum);


% save view and quit
mrSaveView(thisView);
end