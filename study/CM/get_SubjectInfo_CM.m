function subjectInfo = get_SubjectInfo_CM(iSub)
    %
    %   usage: getSubjectInfo_CM(iSub)
    %      by: Ben Gurer
    %    date: 23/11/2017
    % purpose: get subject information for cortical magnificaiton
    %   input: subject number (iSub)
    %
    

% define subject info
% Subject info
subject{1} = '03644_012';
% niftiBaseName{1} = 'pRFpilot2_';
niftiBaseName{1} = '03644_012_';
T2star{1} = '16';
refScan{1} = '15'; % scan before t2 structural
wholeheadPSIR{1} = [];
distCorrectionRefSparse{1} = {'17','18'};
distCorrectionRefCont{2} = {'20','21'};
% scanList = [13,15,19,22];
% freeSurferName{1} = 'kkPSIR_reorient_p7';
freeSurferName{1} = '03644_010';
% distCorrectionRef = {'20','21'};
flatmapName{1} = {'x03644_010_left_WM_invFNIRT_03644_012_Flat_37_131_130_Rad55','x03644_010_right_WM_invFNIRT_03644_012_Flat_140_132_128_Rad55'};
flatmapRotation{1} = {120,340};
% kkPSIR_reorient_p7_right_WM_invFNIRT_03644_012_Flat_249_155_173.off
conditionOrder{1} = {[1,3],[2,4]};

subject{2} = '12013_001';
niftiBaseName{2} = 'cm_12013_001_';
psirNiftiBaseName{2} = 'cm_12013_001';
T2star{2} = '13';
refScan{2} = '8'; % scan before t2 structural
wholeheadPSIR{2} = '16';
distCorrectionRefSparse{2} = {'9','10'};
distCorrectionRefCont{2} = {'11','12'};
freeSurferName{2} = '12013_001';
sparseScans{2} =  {'7','14'};
contScans{2} =  {'8','15'};
flatmapName{2} = {'x12013_001_left_WM_invFNIRT_12013_001_Flat_56_109_74_Rad55','x12013_001_right_WM_invFNIRT_12013_001_Flat_187_124_69_Rad55'};
flatmapRotation{2} = {250,300};
conditionOrder{2} = {[1,3],[2,4]};

subject{3} = '12022_001';
psirSubject{3} = '12022_001';
niftiBaseName{3} = 'cm_12022_001_';
psirNiftiBaseName{3} = 'cm_12022_001';
T2star{3} = '13';
refScan{3} = '08'; % scan before t2 structural
wholeheadPSIR{3} = '17';
distCorrectionRefSparse{3} = {'09','10'};
distCorrectionRefCont{3} = {'11','12'};
freeSurferName{3} = '12022_001';
sparseScans{3} =  {'7','15'};
contScans{3} =  {'8','16'};
flatmapName{3} = {'x12022_001_left_WM_invFNIRT_12022_001_Flat_64_114_84_Rad55','x12022_001_right_WM_invFNIRT_12022_001_Flat_191_114_85_Rad55'};
flatmapRotation{3} = {50,340};
conditionOrder{3} = {[1,3],[2,4]};

subject{4} = '12023_001';
psirSubject{4} = '12023_001';
niftiBaseName{4} = 'cm_12023_001_';
psirNiftiBaseName{4} = 'cm_12023_001';
T2star{4} = '14';
refScan{4} = '9'; % scan before t2 structural
wholeheadPSIR{4} = '17';
distCorrectionRefSparse{4} = {'10','11'};
distCorrectionRefCont{4} = {'12','13'};
freeSurferName{4} = '12023_001';
sparseScans{4} =  {'8','15'};
contScans{4} =  {'9','16'};
flatmapName{4} = {'x12023_001_left_WM_invFNIRT_12023_001_Flat_61_126_66_Rad55','x12023_001_right_WM_invFNIRT_12023_001_Flat_196_129_70_Rad55'};
flatmapRotation{4} = {[],[]};
conditionOrder{4} = {[1,3],[2,4]};

subject{5} = '11108_006';
niftiBaseName{5} = 'cm_11108_006_';
psirNiftiBaseName{5} = 'cm_11108_006';
T2star{5} = '17';
refScan{5} = '12'; % scan before t2 structural
wholeheadPSIR{5} = '20';
distCorrectionRefSparse{5} = {'13','14'};
distCorrectionRefCont{5} = {'15','16'};
freeSurferName{5} = '11108_006';
sparseScans{5} =  {'08','18'};
contScans{5} =  {'12','19'};
% flatmapName{5} = {'x11108_006_left_WM_Flat_56_111_77_Rad55','x11108_006_right_WM_Flat_179_127_72_Rad55'};
flatmapName{5} = {'x11108_006_left_WM_invFNIRT_Flat_56_111_77_Rad55','x11108_006_right_WM_invFNIRT_Flat_179_127_72_Rad55'};
flatmapRotation{5} = {310,60};
conditionOrder{5} = {[1,3],[2,4]};

subject{6} = '11020_002';
niftiBaseName{6} = 'cm_11020_002_';
psirNiftiBaseName{6} = 'cm_11020_002';
T2star{6} = '14';
refScan{6} = '09'; % scan before t2 structural
wholeheadPSIR{6} = '19';
distCorrectionRefSparse{6} = {'10','11'};
distCorrectionRefCont{6} = {'12','13'};
freeSurferName{6} = '11020_002';
sparseScans{6} =  {'08','15'};
contScans{6} =  {'09','18'};
flatmapName{6} = {'x11020_002_left_WM_invFNIRT_11020_002_Flat_74_112_61_Rad55','x11020_002_right_WM_invFNIRT_11020_002_Flat_185_125_75_Rad55'};
flatmapRotation{6} = {0,330};
conditionOrder{6} = {[1,3],[2,4]};

% Import raw scans from subject 9 post motion correction
% Import to one group
subject{7} = '08773';
niftiBaseName{7} = 'cm_08773_007_';
psirNiftiBaseName{7} = 'cm_08773_007';
T2star{7} = '16';
refScan{7} = '10'; % scan before t2 structural
wholeheadPSIR{7} = '17';
distCorrectionRefSparse{7} = {'14','15'};
distCorrectionRefCont{7} = {'21','22'};
freeSurferName{7} = '08773_007';
sparseScans{7} =  {'10','23'};
contScans{7} =  {'20',[]};
flatmapName{7} = {'x08773_007_left_WM_invFNIRT_08773_007_Flat_71_105_84_Rad55','x08773_007_right_WM_invFNIRT_08773_007_Flat_186_111_86_Rad55',};
flatmapRotation{7} = {200,340};
conditionOrder{7} = {[1,3],[4,6]};

% subject 9 is subject 7's repeat session of 3 functional scans
subject{9} = '08773_008';
niftiBaseName{9} = 'cm_08773_008_';
T2star{9} = '11';
refScan{9} = '10'; % scan before t2 structural
distCorrectionRefSparse{9} = {'12','13'};
distCorrectionRefCont{9} = {'14','15'};
sparseScans{9} =  {'10',[]};
contScans{9} =  {'09','18'};
flatmapName{9} = {};

subject{8} = '09933_005';
niftiBaseName{8} = 'cm_09933_005_';
psirNiftiBaseName{8} = 'cm_09933_005';
T2star{8} = '15';
refScan{8} = '14'; % scan before t2 structural
wholeheadPSIR{8} = '14_scan1';
distCorrectionRefSparse{8} = {'07','08'};
distCorrectionRefCont{8} = {'09','10'};
freeSurferName{8} = '09933_005';
sparseScans{8} =  {'07_scan1','13'};
contScans{8} =  {'06','14'};
flatmapName{8} = {'x09933_005_left_WM_invFNIRT_09933_005_Flat_71_105_79_Rad55','x09933_005_right_WM_invFNIRT_09933_005_Flat_190_115_90_Rad55',};
flatmapRotation{8} = {[],[]};
conditionOrder{8} = {[1,3],[2,4]};

%% output selected subjects daya
subjectInfo = struct();
subjectInfo.subjectID = subject{iSub};
subjectInfo.freeSurferName = freeSurferName{iSub};
subjectInfo.flatmapNames = flatmapName{iSub};
subjectInfo.flatmapRotation = flatmapRotation{iSub};
subjectInfo.niftiBaseName = niftiBaseName{iSub};
subjectInfo.T2star = T2star{iSub};

% subjectInfo = struct();
% subjectInfo.subjectID = subject{iSub};
% subjectInfo.niftiBaseName = niftiBaseName{iSub};
% subjectInfo.wholeheadMPRAGE = wholeheadMPRAGE{iSub};
% subjectInfo.fMRIScans = fMRIScans{iSub};
% subjectInfo.nScans = nScans(iSub);
% subjectInfo.refScan = refScan{iSub};
% subjectInfo.T2star = T2star{iSub};
subjectInfo.conditionOrder = conditionOrder{iSub};
% subjectInfo.freeSurferName = freeSurferName{iSub};
% subjectInfo.flatmapNames = flatmapName{iSub};
