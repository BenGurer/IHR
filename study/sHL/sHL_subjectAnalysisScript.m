%% sHL_subjectAnalysisScript
% Scripted subject analysis for simulated hearing loss study
% Study dates: 2017 - 2018

% Info:
% Uses mrTools functions for data analysis, processing and visulisation.
% Calls subfunctions from IHR repository for data analysis, processing and visulisation.

% Aims:
% Measure cortical tonotopic properties in the presense, and absense, of hearing loss
% Compare normal hearing and hearing loss estimates
% Compare analysis methods

%% Tasks:
% load subject data
% perform GLM analysis - General Linear Modelling
% perform pRF analysis - Population Receptive Field modelling 
% calculate gradient reversals
% define ROIs
% export data from volume to flatmap space
% average data over cortical depths
% get cortical distances for tonotopic estimates between gradient reversals
% save data to structure for group analysis

%% save structure
% data.roi.dataset.analysis.voxelParameter
% ie. data.Left.scan1.glm_boxcar_nCons8.R2

%% close and clear everything
clear all; close all; clc

%% Get study parameters
[stimInfo, glmInfo, pRFInfo, Info, plotInfo] = sHL_setupStudyParams;
% stimulus info, condition names, nummber of subjects etc

% use ispc to set data directory
if ispc
    Info.dataDir = 'E:\OneDrive - The University of Nottingham\data';
else
    Info.dataDir = '/Volumes/data_PSY/OneDrive - The University of Nottingham/data';
end

% set variable to be ' to use with eval
q = char(39);

%% set logicals to what we want to do
doLoadView = 1;
doAntomy = 0;
doPSIR = 0;
doGLMbc = 0;
doMakeFlatmaps = 0;
doMakeARrois = 0;
doDeconv = 0;
doHRF = 0;
doGLMdg = 0;
doConvertOverlays = 0;
doDeleteOverlays = 0;
doConvertvol2FlatAvDepth_GLM = 0;
doGradientReversal_GLM = 0;
doROIsGLM = 0;
doGetDATA_GLM = 0;
doROIRestrict_GLM = 0;
dopRF = 0;
doConvertvol2FlatAvDepth_pRF = 1;
doGradientReversal_pRF = 1;
doROIspRF = 0;
doGetAndRestrictDATA_pRF = 1;
doCorticalMagnification = 0;
doGRrois = 0;

%% define subjects
iSubs2Run = [1,2,3,4,5,6,7,8];
iSubs2Run = 3

%% loop over subjects
% use for loop to repeat analysis for each subject
for iSub = 1:length(iSubs2Run)
    
    %% Setting up subject
    % Get subject info
    disp('Getting subject info...')
    subjectInfo = get_SubjectInfo_CM(iSubs2Run(iSub));
    % Subject ID, flatmap names etc
    
    % define file name for data
    saveName = [subjectInfo.subjectID '_data.mat'];
    
    % move to subject folder
    cd(fullfile(Info.dataDir,Info.studyDir,subjectInfo.subjectID));
    
    % load data if it exists
    if exist(saveName, 'file')
        disp('Loading data...')
        load(saveName)
        if isfield(data, 'hrf')
            glmInfo.hrfParamsDoubleGamma = data.hrf.x_doubleGamma;
            pRFInfo.hrfParamsGamma = data.hrf.x_Gamma;
            pRFInfo.hrfParamsDiffofGamma = data.hrf.x_dGamma;
        end
    else
        data = struct;
    end
    
    %% open mrLoadRet and get view
    if doLoadView
        disp('Loading mrView...')
        
        % clear any other views
        clear thisView
        
        % open view in mrLoadRet
        mrLoadRet
        
        % get view
        thisView = getMLRView;
        
        % refresh mrLoadRet view
        refreshMLRDisplay(thisView.viewNum);
        
    end
    
    %% Import anatomy
    if doAntomy
        thisView = script_importAnatomy(thisView,Info,subjectInfo);
        % load in:
        % reference EPI
        % High resolution in-plane T2*
        % surfaces
        
    end
    
       %% Make flatmaps
    if doMakeFlatmaps
        disp('Make flatmaps for each hemisphere...')
        disp('centre on HG (use R2 and f-test to guide)')
        disp('use default name')
        disp('Params: Radius=55; Resolution=3; Method=mrFlatMesh')
        disp('rotate flatmaps for easy viewing')
        disp('Press F5 when done...')
        keyboard
        % radius = 55
        % centre on HG (use R2 and f-test to guide)
        % use default names
        % resolution = 3; method = mrFlatMesh
        % rotate flatmaps for easy viewing (do before exporting to flatmap space)
    end
    
    %% GLM Analysis with Double Gamma HRF
    % Estiamting tonotopic properties
    if doGLMdg
        disp('Performing GLM analysis (HRF=double gamma)')
        % Perform GLM analysis with double gamma model of HRF
        % hrf = Box Car
        % All runs
        % [thisView] = script_glmAnalysis(thisView,glmInfo,groupNames,hrfModel,runSplitHalf,runTonotopic)
        glmInfo.hrfParamsDoubleGamma = data.hrf.x_doubleGamma;
        thisView = script_glmAnalysis(thisView,glmInfo,glmInfo.groupNames,{'hrfDoubleGamma'},1,1);
    end
    
    %% Auditory Responsive ROI CREATION
    if doMakeARrois
        % create ROIs with the names:
        % LeftAR and RightAR: Group=Sparse; Analysis=glm_hrfdoublegamma; Overlay=f-test (FDR adjusted)- set min to 0.05
        % Create ROI - continuous voxels;
        % ROIs>transform>expandROI([3 3 3])(convolves ROI with a sphere with a diameter of 3^3 voxels)
        % name=(Left or Right)ARexp
        % Project through depths 0.3 to 0.7 to remove voxels outside of grey matter for ALL ROIs
        % combine LeftAR and RightAR = AR & combine LeftARexp and RightARexp = ARexp
        
        % set group
        groupName = glmInfo.groupNames{1};
        if viewGet(thisView,'curgroup') ~= viewGet(thisView,'groupNum',groupName)
            thisView = viewSet(thisView,'curgroup',groupName);
        end
        % set analysis
        analysisName = glmInfo.analysisNames{1};
        if viewGet(thisView,'curAnalysis') ~= viewGet(thisView,'analysisNum',analysisName)
            thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',analysisName));
        end
        % set base
        flatmapName = subjectInfo.flatmapNames{1};
        if viewGet(thisView,'curbase') ~= viewGet(thisView,'basenum',flatmapName)
            thisView = viewSet(thisView,'curbase',viewGet(thisView,'basenum',flatmapName));
        end
        
        % set overlay
        overlayNum = viewGet(thisView,'overlayNum','FDR-adjusted P [F (fTest - all conditions)]');
        thisView = viewSet(thisView,'curOverlay',overlayNum);
        
        refreshMLRDisplay(thisView.viewNum);
        
        disp('create ROIs with the names:')
        disp('LeftAR and RightAR: Group=Sparse; Analysis=glm_hrfboxcar; Overlay=f-test (FDR adjusted)- set min to 0.05')
        disp('Create ROI - continuous voxels;')
        disp('hit F5 when done')
        
        keyboard
        
        thisView = getMLRView;
        
        % get ROI numbers
        roiList(1) = viewGet(thisView,'roinum','LeftAR');
        roiList(2) = viewGet(thisView,'roinum','RightAR');
        
        % get rois
        rois = viewGet(thisView,'roi',roiList);
        
        
        % combine ROIs
        newName = 'AR';
        roi1 = 'LeftAR';
        roi2 = 'RightAR';
        thisView = combineROIs(thisView,roi1,roi2,'union',newName);
        
        for iSide = 1:2
            clear outputRoi
            outputRoi = expandROI(rois(iSide),[3 3 3],'sphere');
            outputRoi.name = [outputRoi.name 'exp'];
            thisView = viewSet(thisView,'newROI',outputRoi);
        end
        
        refreshMLRDisplay(thisView.viewNum);
        
        disp('check ROIs include HG and then project between 0 - 1 (ALL) cortical depths.')
        disp('Restrict to any tonotopic overlay (scan data)(cmd+x)')
        disp('type return in command line when done')
        
        keyboard
        
        % get view with latest rois
        thisView = getMLRView;
        
        % combine ROIs
        newName = 'ARexp';
        roi1 = 'LeftARexp';
        roi2 = 'RightARexp';
        thisView = combineROIs(thisView,roi1,roi2,'union',newName);
        
        % save view
        mrSaveView(thisView);
        
    end