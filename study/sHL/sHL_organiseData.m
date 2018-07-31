function sHL_organiseData(Info, subjectInfo)
    %
    %   usage: sHL_organiseData(Info,subjectInfo)
    %      by: Ben Gurer
    %    date: 07/12/2017
    % purpose: organise and move raw data for hearing loss simulation study
    %   input: Info: information about study; subjectInfo: information
    %   about subjects
    %


%% Move data from scanner to study/subject folders
cd([Info.dataDir '/scanner/' subjectInfo.subjectID]) 

% convert from Par/Rec to niffti format
!ptoa -f -q -nii *.PAR

if ~isempty(subjectInfo.wholeheadMPRAGE)
    % Move whole head PSIR
    mkdir(fullfile(Info.dataDir, 'Anatomy','originals', subjectInfo.freeSurferName));
    movefile(fullfile(Info.dataDir, 'scanner', subjectInfo.subjectID, [subjectInfo.niftiBaseName subjectInfo.wholeheadMPRAGE '*.nii']), fullfile(Info.dataDir,'Anatomy','originals',subjectInfo.freeSurferName));
    cd(fullfile(Info.dataDir, 'Anatomy/originals/', subjectInfo.freeSurferName));        
    %RUN recon-all DIRECTLY IN TERMINAL 
    % if not in folder
    fprintf(['recon-all -subjid ' subjectInfo.freeSurferName ' -i ' fullfile(Info.dataDir, 'Anatomy','originals', subjectInfo.freeSurferName) '/' [subjectInfo.niftiBaseName subjectInfo.wholeheadMPRAGE '_1.nii'] ' -all\n']);
    % if in folder
    fprintf(['recon-all -subjid ' subjectInfo.freeSurferName ' -i ' [subjectInfo.niftiBaseName subjectInfo.wholeheadMPRAGE '_1.nii'] ' -all']);
end

% make subject directory
mkdir(fullfile(Info.dataDir,Info.studyDir,subjectInfo.subjectID));
cd(fullfile(Info.dataDir,Info.studyDir,subjectInfo.subjectID));

mkdir('Etc')
mkdir('Anatomy')
mkdir('Raw')
mkdir('Raw/TSeries')
!mkdir FNIRT
!mkdir Distorted

% Move in-plane T2 star
movefile(fullfile(Info.dataDir,'scanner',subjectInfo.subjectID,[subjectInfo.niftiBaseName subjectInfo.T2star '*.nii']),fullfile(Info.dataDir,Info.studyDir,subjectInfo.subjectID,'Anatomy'))

% check for filenames below 10 and add a zero before the number - make sure linking of stim files works later
cd Raw/TSeries/
% scanFiles = dir;
% for id = 1:length(scanFiles)
%     str = scanFiles(id).name;
%     strParts = strsplit(str,'_');
%     if length(strParts) > 1
%     checkScanNum = char(strParts(4));
% %         numStrParts = length(strParts);
%         if 10>str2num(checkScanNum) && numel(checkScanNum)==1
%             newName = [subjects{iSubj} '_0' char(strParts(4)) '.nii']
%             movefile(scanFiles(id).name,newName);
%         else
%                         newName = [subjects{iSubj} '_' char(strParts(4)) '.nii']
%             movefile(scanFiles(id).name,newName);
%         end
%    
%     end
% end

%% move fMRI EPIs
movefile(fullfile(Info.dataDir,'scanner',subjectInfo.subjectID,[subjectInfo.niftiBaseName 'WIP_73DYN_fMRI_*.nii']), fullfile(Info.dataDir,Info.studyDir,subjectInfo.subjectID,'Distorted'))

%% move Distortion Correction Scans
movefile(fullfile(Info.dataDir,'scanner',subjectInfo.subjectID,[subjectInfo.niftiBaseName 'WIP_5DYN_AP_*.nii']), fullfile(Info.dataDir,Info.studyDir,subjectInfo.subjectID,'Distorted'))

cd ../..