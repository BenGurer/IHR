function sHL_preprocess(Info, subjectInfo,distortionCorrection)
%
%   usage: sHL_preprocess(Info,subjectInfo)
%      by: Ben Gurer
%    date: 07/12/2017
% purpose: pre-process data for hearing loss simulation study
%   input: Info: information about study; subjectInfo: information
%   about subjects
%

%% Pre-process
% move to subject folder
cd(fullfile(Info.dataDir,Info.studyDir,subjectInfo.subjectID));

%% copy freesurfer code to here

%% fsl topup Distortion Correction

% move to Distorted folder
cd Distorted

if distortionCorrection
    
    % select fsltopup distortion information file
    if subjectNumber >= 5
        distInfo = fullfile(Info.dataDir,Info.studyDir,'base_fat_echo_shifted_params.txt');
        % distInfo = 'base_fat_echo_shifted_params.txt'
    else
        distInfo = fullfile(Info.dataDir,Info.studyDir,'base_fat_shifted_AP_params.txt');
    end
    
    % stick distortion correction scans together and use as input to fsltop
    system(['fslmerge -t topup_input.nii ', ' ', [subjectInfo.niftiBaseName 'WIP_5DYN_AP_*_1_modulus.nii']])
    system(['topup --imain=topup_input.nii --datain=' distInfo ' --subsamp=1 --fwhm=0  --out=topup_usable_output --iout=topup_check_output']);
    
    % apply fsltop up to fMRI EPI scans
    for iScan = 1:subjectInfo.nScans
        %     fileName = [subjectInfo.niftiBaseName 'WIP_73DYN_fMRI_0' num2str(iScan) '_' num2str(subjectInfo.fMRIScans{iScan}) '_1_modulus'];
        fileName = [subjectInfo.niftiBaseName 'WIP_73DYN_fMRI_0' num2str(iScan)];
        system(['applytopup --imain=' fileName '_' num2str(subjectInfo.fMRIScans{iScan}) '_1_modulus.nii --datain=' distInfo ' --inindex=1 --topup=topup_usable_output --out=' fileName '_unwarped.nii --method=jac']);
    end
    
    % copy to Raw group time series folder for mrLoadRet
    system('cp *unwarped.nii ../Raw/TSeries/');
    
else
    system('cp *73DYN*modulus.nii ../Raw/TSeries/');
end


%% Align scans
cd ../Raw/TSeries/

% crop last frame of reference EPI
system(['fslroi *' subjectInfo.refScan '*.nii lastFrameEPI ' num2str(Info.epiDims(4)-1) ' 1']);


% move last frame of ref EPI to FNIRT
!mv lastFrameEPI.nii ../../FNIRT

% move to anatomy folder
cd ../../Anatomy/

% align T2* in-plane to T1 wholehead using mrAlign
mrAlign
keyboard
% if subject has been scanned before - first align MPRAGE to PSIR
% align hi res t2* to MPRAGE
% Load MPRAGE as destination and ' Set base cooridinate frame'
% Load hi res t2* as source
% manually align then go to
% Compute Alignment->Advanced Alignment Menu
% then select Reverse Contrast (T2*)
% Click 'Compute Coarse Aligment' then
% Click 'Compute fine Aligment'
% save alignment
% doing it here means cropped versions with have the correct alignement -
% skull stripping removes the s-form matrix so using set alignment to
% identity and save (t2*_stripped to t2*)

% align EPI to T2* in-plane
% crop T2* in-plane to use with FNIRT and FLIRT
T2starFile = [subjectInfo.niftiBaseName subjectInfo.T2star '_1_modulus'];
T2StarDims_Index = [400 400 29];
EPIDims_Index = [128 128 24];

T2StarDims_mm = [0.55 0.55 1.5];
EPIDims_mm = [1.359375  1.359375 1.5];

T2StarDims = T2StarDims_Index .* T2StarDims_mm;
EPIDims = EPIDims_Index .* EPIDims_mm;

stackDiff = T2StarDims - EPIDims;

T2StartStart = round((stackDiff ./ T2StarDims_mm)/2);
T2StartEnd = round(EPIDims ./ T2StarDims_mm);

system(sprintf('fslroi %s %s_crop %d %d %d %d %d %d 0 1',T2starFile,T2starFile,T2StartStart(1), T2StartEnd(1),T2StartStart(2), T2StartEnd(2),T2StartStart(3), T2StartEnd(3))); %% cropping t2 inplane structural image to match size of functional data
system(sprintf('cp %s_crop.nii ../FNIRT',T2starFile));

system(['bet ' T2starFile ' ' T2starFile '_stripped -f .2 -Z'])

mrAlign
keyboard
% skull stripping removes the s-form matrix so using set alignment to
% identity and save (t2*_stripped to t2*)

% Use FLIRT AND FNIRT to register functional to structural
cd ../FNIRT
%run FNIRT
system(sprintf('fnirtEpi2T2star lastFrameEPI %s_crop -separateFLIRT',T2starFile));

lastFrameNames = struct2cell(dir('lastFrameEPI*.nii'));
lastFrameNames = lastFrameNames(1,:);
alignFunctional2HighResT2Star(sprintf('lastFrameEPI2%s_crop_resampled.affmat',T2starFile),[T2starFile '_crop.nii'],lastFrameNames);

functionalNames = struct2cell(dir('../Raw/TSeries/*.nii'));
functionalNames = functionalNames(1,:);
for i = 1:length(functionalNames)
    functionalNames{i} = ['../Raw/TSeries/' functionalNames{i}];
end
alignFunctional2HighResT2Star(sprintf('lastFrameEPI2%s_crop_resampled.affmat',T2starFile),[T2starFile '_crop.nii'],functionalNames);

%% Pre-processing DONE!

cd ../

end