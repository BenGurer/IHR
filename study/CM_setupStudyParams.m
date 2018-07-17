function [stimInfo, glmInfo, pRFInfo, Info, plotInfo] = CM_setupStudyParams
    %
    %   usage: sHL_setupStudyParams
    %      by: Ben Gurer
    %    date: 11/10/2017
    % purpose: setup study info and params for hearing loss simulation study
    %   input: n/a
    %  output: stimInfo, stimulus information; glmInfo, glm analysis information
    %
%% Define output structures
stimInfo = struct();
glmInfo = struct();
pRFInfo = struct();
Info = struct();
plotInfo = struct();
    
%% Get stimulus properties
% get stimulus frequencies
stimInfo.lowFreqkHz = 0.1;
stimInfo.highFreqkHz = 8;
stimInfo.nStim = 32;

[stimInfo.stimFreqs, stimInfo.stimFreqs_bin, stimInfo.stimFreqs_mv, stimInfo.stimNERBs, stimInfo.stimNERBs_bin, stimInfo.stimNERBs_mv] = convertStimIDtoFrequency(stimInfo.lowFreqkHz,stimInfo.highFreqkHz,stimInfo.nStim);

stimInfo.sizes = [8 29 32];

% below is too allow older code to still work
stimInfo.stimNames.all = stimInfo.stimFreqs;
stimInfo.stimNames.bin = stimInfo.stimFreqs_bin;
stimInfo.stimNames.mvstim = stimInfo.stimFreqs_mv;

%% get stimulus senssation level
% [stimLevel_SL, maskingLevel] = calStimulusSensationLevel(stimInfo.stimNames.all);
% % bin sensation level
% binSize = 4;
% stimLevel_SL_bin = zeros(1,length(stimLevel_SL)/binSize);
% c = 1;
% for i = 1:length(stimLevel_SL)/binSize
%     stimLevel_SL_bin(i) = mean(stimLevel_SL(c:c + (binSize-1)));
%     c = c + binSize;
% end
% 
% % moving average sensation level
% nBins = 8;
% windowAvSize = length(stimLevel_SL)/nBins;
% 
% if isreal(windowAvSize) && rem(windowAvSize,1)==0
%     loopLength = (length(stimLevel_SL) - windowAvSize) +1;
%     stimLevel_SL_mv = zeros(1,loopLength);
%     for i = 1:loopLength
%         stimLevel_SL_mv(i) = mean(stimLevel_SL(i:i+windowAvSize-1));
%     end
% else
%     error('Moving average window not an integer')
% end
% 
% % save in structure
% stimInfo.stimLevel_SL = stimLevel_SL;
% stimInfo.stimLevel_SL_bin = stimLevel_SL_bin;
% stimInfo.stimLevel_SL_mv = stimLevel_SL_mv;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Setup GLM analysis
% save in glmInfo structure
glmInfo.hrfModel = {'hrfBoxcar', 'hrfDoubleGamma'};
glmInfo.groupNames = {'ConcatenationSparse', 'ConcatenationCont'};
% define vocel property overlay names - index max, centriod, spread, julienCentriod, julienTuningWidth
glmInfo.voxelPropertyNames = {'Centriod','Spread','julien_pCF','julien_pTW','indexMax'};
glmInfo.nScans = 4;
glmInfo.nStim = [32, 8];
% glmInfo.analysisNames_Scans = cell(1,(glmInfo.nScans.*length(glmInfo.nStim)).*length(glmInfo.hrfModel));
c = 0;
d = 0;
f = 0;
for iScan = 1:glmInfo.nScans
    for iStim = 1:length(glmInfo.nStim)
        for iHRF = 1:length(glmInfo.hrfModel)
            %         iHRF = 1; % only perform boxcar analysis on individual scans
           c = c + 1;
%             glmInfo.analysisNames_Scans{c+d+f} = ['glm_' glmInfo.hrfModel{iHRF} '_nCons_' mat2str(glmInfo.nStim(iStim)) '_Scan_' mat2str(iScan)];
%             glmInfo.analysisBaseNames_Scans{c+d+f} = ['glm_' glmInfo.hrfModel{iHRF} '_nCons_' mat2str(glmInfo.nStim(iStim))];
%             glmInfo.analysisScanNum {c+d+f} = iScan;
            glmInfo.analysisNames_Scans{c + d} = ['glm_' glmInfo.hrfModel{iHRF} '_nCons_' mat2str(glmInfo.nStim(iStim)) '_Scan_' mat2str(iScan)];
            glmInfo.analysisBaseNames_Scans{c + d}= ['glm_' glmInfo.hrfModel{iHRF} '_nCons_' mat2str(glmInfo.nStim(iStim))];
            glmInfo.analysisScanNum{c + d} = iScan;
            glmInfo.analysisNStim{c + d} = glmInfo.nStim(iStim);
        end
%         c = 0;
%         d = d + length(glmInfo.hrfModel); % remove when only using boxcar
%                 d = d + 1;
    end
    c = 0;
%     f = f + length(glmInfo.nStim);
d = d + length(glmInfo.hrfModel) + length(glmInfo.nStim);
end

glmInfo.analysisNames_Groups = cell(1,length(glmInfo.groupNames)*length(glmInfo.groupNames));
c = 0;
d = 0;
for iGroup = 1:length(glmInfo.groupNames)
    for iHRF = 1:length(glmInfo.hrfModel)
        c = c + 1;
        glmInfo.analysisNames_Groups{c+d} = ['glm_' glmInfo.hrfModel{iHRF}];
    end
    c = 0;
    d = d + length(glmInfo.hrfModel);
end

glmInfo.scanGroupName = 'MotionComp';

c = 0;
for iStim = 1:length(glmInfo.nStim)
    for iHRF = 1:length(glmInfo.hrfModel)
        c = c + 1;
        glmInfo.analysisNames_nCons{c}= ['glm_' glmInfo.hrfModel{iHRF} '_nCons_' mat2str(glmInfo.nStim(iStim))];
    end
end

glmInfo.analysisNames = cell(1,length(glmInfo.hrfModel));
for iHRF = 1:length(glmInfo.hrfModel)
    glmInfo.analysisNames{iHRF} = ['glm_' glmInfo.hrfModel{iHRF}];
end

%% pRF information
% what info do I need?

% pRFInfo.name = pRF;

% pRFinfo.analysisNames_Groups{iGroup}{iAnal}
% stimulusWeighting = {'None','SL_level','BOLD','fit'};
% pRFinfo.analysisNames_Groups{1}{1} = {'pRF'};
% pRFInfo.stimulusWeighting = {'None','SL_level','BOLD','fit'};
% for iWeight = 1:length(stimulusWeighting)
%     pRFinfo.analysisNames_Groups{2}{iWeight} = {['pRF_' stimulusWeighting{iWeight}]};
% end
    
pRFInfo.analysisNames_Groups{1}{1} = 'pRF';
pRFInfo.analysisNames_Groups{2}{1} = 'pRF';
% pRFInfo.stimulusWeighting{1} = {'None'};
% % pRFInfo.stimulusWeighting{2} = {'None','SL_level','BOLD','fit'};
% 
% pRFInfo.stimulusWeighting{2} = {'BOLD'};
% % pRFInfo.stimulusWeighting{2} = {'None','SL_level'};
% for iWeight = 1:length(pRFInfo.stimulusWeighting{2})
%     pRFInfo.analysisNames_Groups{2}{iWeight} = ['pRF_' pRFInfo.stimulusWeighting{2}{iWeight}];
% end
pRFInfo.pRFrestrictROI = 'ARexp';
pRFInfo.pRFrois = {'Left_AR_exp','Right_AR_exp'}; % ROIs to restrict pRF analysis - expanded around AC
pRFInfo.pRFOverlayNames = {'r2','PrefCentreFreq','rfHalfWidth'};
pRFInfo.pRFgradientReversalOverlay = 2;



%% Info - save general info needed to struct

if ispc
    Info.dataDir = 'N:/data';
elseif isunix
    Info.dataDir = '/home/beng/data';
end
Info.studyDir = 'CorticalMagnification';
%% Order of condition runs {[ConA Run1,ConA Run2],[ConB Run1,ConB Run2]}

Info.conditionRunIndex = {[1,3],[2,4]};
% Info.conditionRunIndex = {[2,4],[6,8]};

Info.ConATrue = 1;

%% Define ROI names to create
% Info.ROInames = {'RightAC','RightPosAC','RightAntAC','LeftAC','LeftPosAC','LeftAntAC','AC'};

% Must match order of sides
Info.ROInames = {'LeftGR','RightGR'};
Info.ROIbasenames = {'GR'};
Info.LeftROInames = {'LeftGR_GLM','LeftGRa_GLM','LeftGRp_GLM'};
Info.RightROInames = {'RightGR_GLM','RightGRa_GLM','RightGRp_GLM'};
% Info.ROInames_flat = {'RightAC_flat','LeftAC_flat'};

Info.epiDims = [128 128 24 73]; % dims of functional scans

Info.sides = {'left','right'};
Info.Sides = {'Left','Right'};

Info.gradReversalInfo.analysisBase = glmInfo.analysisNames_Groups{1};
Info.gradReversalInfo.groupBase = glmInfo.groupNames{1};
Info.gradReversalInfo.overlayBase = 41; % 39

%% Define what to plot
% cell array of analysis to plot
% what plots to present
% use structure to group correct stim names with it
% plotInfo = struct();
% list of analysis to plot from
plotInfo.ROIplotList = {['roiAnalysis_' glmInfo.analysisBaseNames_Scans{1}], ['roiAnalysis_' glmInfo.analysisBaseNames_Scans{2}]};

% plotLogic = [roiAv, roiTWav, ratio]
% use logicals to tell function what to plot
plotInfo.plotLOGIC.ROI_bin = [0, 1, 0];
plotInfo.plotLOGIC.ROI_mv = [1, 0, 0];
plotInfo.plotLOGIC.ROI_all = [0, 0, 0];