%% CM_subjectAnalysisScript
% Scripted subject analysis for Comparisions at 7T and Cortical Magnification (CM) studies
% Study dates: 2015 - 2018

% Aims:
% Measure cortical magnification
% Compare acquisition protocols
% Compare analysis methods

%% Tasks:
% load subject data
% perform GLM analysis
% perform pRF analysis
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
[stimInfo, glmInfo, pRFInfo, Info, plotInfo] = CM_setupStudyParams;
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
doCorticalMagnification = 1;
doGRrois = 1;

%% define subjects
% iSubs2Run = [1,2,3,4,5,6,7,8];
iSubs2Run = [1,2,3,6,7,8];

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
    
    %% Import PSIR
    if doPSIR
        disp('Importing PSIR information')
        % rotate surfaces and then get view
        % save rotations to subjectInfo function
        thisView = getMLRView;
        
        % Import PSIR
        thisView = importPSIRasOverlay(thisView,Info,subjectInfo);
        
        % import flatmaps without FNIRT coords (after flatmaps have been created)
        thisView = script_importFlatmaps(thisView,Info,subjectInfo);
        
    end
    
    %% GLM analysis %%
    % voxelwise analysis using general linear modelling
    
    % HRF = Boxcar
    % Determine auditory responsive voxels
    
    % HRF = Deconvolution
    % Determine hemodynamic response function
    
    % HRF = Double gamma
    % Determine tonotopic properites
    
    %% GLM Analysis with Boxcar HRF
    % Estiamting auditory responsive voxels
    if doGLMbc
        disp('Performing GLM analysis (HRF=boxcar)')
        % Perform GLM analysis with boxcar model of HRF
        % hrf = Box Car
        % concatenated runs only
        % [thisView] = script_glmAnalysis(thisView,glmInfo,groupNames,hrfModel,runSplitHalf,runTonotopic)
        thisView = script_glmAnalysis(thisView,glmInfo,glmInfo.groupNames,{'hrfBoxcar'},1,1);
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
    
    %% Auditory Responsive ROI CREATION
    if doMakeARrois
        % create ROIs with the names:
        % LeftAR and RightAR: Group=Sparse; Analysis=glm_hrfboxcar; Overlay=f-test (FDR adjusted)- set min to 0.05
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
    
    %% GLM Analysis with deconvolved HRF
    % Estiamting hemodynamic response function
    if doDeconv
        % estimate hrf using deconvolution
        disp('Performing GLM analysis (HRF=deconvolution)')
        % get view so we have the ROIs
        thisView = getMLRView;
        thisView = script_hrfAnalysis(thisView,glmInfo.groupNames{2});
    end
    
    %% get average HRF estimate for Auditory Responsive(AR) ROI
    if doHRF
        % use results for GLM
        % save result to data
        % get view so we have the ROIs and data
        thisView = getMLRView;
        [ thisView, data.hrf.x_doubleGamma, data.hrf.x_Gamma, data.hrf.x_dGamma, data.hrf.estimate, data.hrf.deconv, data.hrf.deconvTW] = script_hrfROIAnalysis(thisView,'AR',glmInfo);
        % save to analysis structures to pass to functions
        glmInfo.hrfParamsDoubleGamma = data.hrf.x_doubleGamma;
        pRFInfo.hrfParamsGamma = data.hrf.x_Gamma;
        pRFInfo.hrfParamsDiffofGamma = data.hrf.x_dGamma;
        
        % save data
        % save so don't need to load again
        disp('saving HRF estimate data...')
        save(saveName,'data','-v7.3');
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
    
    %% get condition names
    % Get condition names to allow us to access and save data
    getConditionNames = cell(1,length(glmInfo.nStim));
    % save condition names
    if exist('data') && isfield(data, 'conditions')
        conditionNames = data.conditions;
    else
        analysisName = glmInfo.analysisNames_Scans{1};
        conditionNames{1} = get_analysisConditionNames(thisView,analysisName,glmInfo.scanGroupName,1);
        analysisName = glmInfo.analysisNames_Scans{3};
        conditionNames{2} = get_analysisConditionNames(thisView,analysisName,glmInfo.scanGroupName,1);
        data.conditions = conditionNames;
    end
    
    %% Convert overlay values from stimulus scale to nERB
    % convert pCF overlays to nERB
    % all tonotopic estimates - groups and scans
    if doConvertOverlays
        thisView = getMLRView;
        % Groups
        for iGroup = 1:length(glmInfo.groupNames)
            % select group
            groupName = glmInfo.groupNames{iGroup};
            if viewGet(thisView,'curgroup') ~= viewGet(thisView,'groupNum',groupName)
                thisView = viewSet(thisView,'curgroup',groupName);
            end
            % Select analysis
            for iAnal = 1:length(glmInfo.analysisNames)
                analysisName = glmInfo.analysisNames{iAnal};
                thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',analysisName));
                % get condition names and create concatenate string
                conNamesString = [];
                for iCon = 1:length(conditionNames{1})
                    stimOne = stimInfo.stimNERBs(1);
                    stimTwo = stimInfo.stimNERBs(2);
                    if iCon == 1
                        conNamesString  = [conNamesString conditionNames{1}{iCon}];
                    else
                        conNamesString  = [conNamesString ',' conditionNames{1}{iCon}];
                    end
                end
                
                % get the overlay numbers of the data we want
                % glmInfo.voxelPropertyNames = {'Centriod','Spread','julien_pCF','julien_pTW','indexMax'};
                overlayNum = zeros(1,length(glmInfo.voxelPropertyNames));
                for i = 1:length(glmInfo.voxelPropertyNames)-1
                    overlayNum(i) = viewGet(thisView,'overlayNum',['Ouput ' num2str(i) ' - weightedMeanStd_CM(' conNamesString ')']);
                end
                overlayNum(end) = viewGet(thisView,'overlayNum',['Ouput 1 - indexMax(' conNamesString ')']);
                
                % convert overlays to nERB
                for i = 1:length(overlayNum)
                    overlayIN = viewGet(thisView,'overlay',overlayNum(i));
                    [ thisView , ~ ] = convertOverlay_GLMCF2NERB(thisView,overlayIN,stimOne,stimTwo,[glmInfo.voxelPropertyNames{i} '_nERB']);
                end
                
            end
            
        end
        
        % Scans
        % Select scan group
        groupName = glmInfo.scanGroupName;
        if viewGet(thisView,'curgroup') ~= viewGet(thisView,'groupNum',groupName)
            thisView = viewSet(thisView,'curgroup',groupName);
        end
        for iScan = 1:glmInfo.nScans
            thisView = viewSet(thisView,'curScan', iScan);
            for iAnal = 1:length(glmInfo.analysisNames_nCons)
                % select analysis
                analysisName = [glmInfo.analysisNames_nCons{iAnal} '_Scan_' mat2str(iScan)];
                thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',analysisName));
                % get condition names and create concatenate string
                if glmInfo.analysisNStim{iAnal} == length(conditionNames{1})
                    overlayFlatNames = cell(1,length(conditionNames{1}));
                    conNamesString = [];
                    for iCon = 1:length(conditionNames{1})
                        stimOne = stimInfo.stimNERBs(1);
                        stimTwo = stimInfo.stimNERBs(2);
                        if iCon == 1
                            conNamesString  = [conNamesString conditionNames{1}{iCon}];
                        else
                            conNamesString  = [conNamesString ',' conditionNames{1}{iCon}];
                        end
                    end
                else
                    stimOne = stimInfo.stimNERBs_bin(1);
                    stimTwo = stimInfo.stimNERBs_bin(2);
                    conNamesString = [];
                    for iCon = 1:length(conditionNames{2})
                        if iCon == 1
                            conNamesString  = [conNamesString conditionNames{2}{iCon}];
                        else
                            conNamesString  = [conNamesString ',' conditionNames{2}{iCon}];
                        end
                    end
                end
                
                % get the overlay numbers of the data we want
                overlayNum = zeros(1,length(glmInfo.voxelPropertyNames));
                for i = 1:length(glmInfo.voxelPropertyNames)-1
                    overlayNum(i) = viewGet(thisView,'overlayNum',['Ouput ' num2str(i) ' - weightedMeanStd_CM(' conNamesString ')']);
                end
                overlayNum(end) = viewGet(thisView,'overlayNum',['Ouput 1 - indexMax(' conNamesString ')']);
                
                % convert overlays to nERB
                for i = 1:length(overlayNum)
                    overlayIN = viewGet(thisView,'overlay',overlayNum(i));
                    [ thisView , ~ ] = convertOverlay_GLMCF2NERB(thisView,overlayIN,stimOne,stimTwo,[glmInfo.voxelPropertyNames{i} '_nERB']);
                end
            end
        end
        mrSaveView(thisView);
    end
    
    %% Convert GLM data to flatmap space and average over cortical depth
    % [thisView, analysisData] = script_covertData2FlatmapSpace(thisView,groupName,analysisName,iScan,overlays,flatmapName)
    
    % delete overlays if no longer needed
    if doConvertvol2FlatAvDepth_GLM
        
        % auto delete all overlays because its a pain to do manually
        if doDeleteOverlays
            thisView = getMLRView;
            for iSide = 1:length(subjectInfo.flatmapNames)
                thisView = viewSet(thisView,'curgroup',[subjectInfo.flatmapNames{iSide}, 'Volume']);
                thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum','combineTransformOverlays'));
                analysisData = viewGet(thisView,'analysis',viewGet(thisView,'analysisNum','combineTransformOverlays'));
                overlayNumbers = 1:length(analysisData.overlays);
                thisView = viewSet(thisView,'deleteoverlay',overlayNumbers);
            end
        end
        
        % export scan data
        for iScan = 1:glmInfo.nScans
            for iAnal = 1:length(glmInfo.nStim)*length(glmInfo.hrfModel)
                analysisName = [glmInfo.analysisBaseNames_Scans{iAnal}, '_Scan_' mat2str(iScan)];
                for iSide = 1:length(subjectInfo.flatmapNames)
                    thisView = script_covertData2FlatmapSpace(thisView,glmInfo.scanGroupName,analysisName,iScan,[],subjectInfo.flatmapNames{iSide});
                end
            end
        end
        
        % export group data
        for iGroup = 1:length(glmInfo.groupNames)
            for iSide = 1:length(subjectInfo.flatmapNames)
                for iAnal = 1:length(glmInfo.analysisNames_Groups)
                    thisView = script_covertData2FlatmapSpace(thisView,glmInfo.groupNames{iGroup},glmInfo.analysisNames_Groups{iAnal},[],[],subjectInfo.flatmapNames{iSide});
                end
            end
        end
        
        % average over cortical depth
        for iSide = 1:length(subjectInfo.flatmapNames)
            thisView = script_averageAcrossDepths(thisView,[],[subjectInfo.flatmapNames{iSide}, 'Volume'],1);
        end
        
        % save view
        mrSaveView(thisView);
        
    end
    
    %% GLM grandient reversals
    % calculate gradient reversals using GLM (double gamma) analysis.
    % pCF estimation = Juliens debiased method
    if doGradientReversal_GLM
        groupName = glmInfo.groupNames{1}; % 'ConcatenationSparse';
        if viewGet(thisView,'curgroup') ~= viewGet(thisView,'groupNum',groupName)
            thisView = viewSet(thisView,'curgroup',groupName);
        end
        analysisName = glmInfo.analysisNames{2}; % 'glm_hrfDoubleGamma';
        thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',analysisName));
        conNamesString = [];
        
        % get overlay number for debiased prefered frequency in number of ERBS overlay
        overlayNum = viewGet(thisView,'overlayNum','julien_pCF_nERB');
        
        % calcuate gradient reversals
        thisView = script_flatMapAnalysis(thisView,Info,subjectInfo,Info.gradReversalInfo.groupBase, glmInfo.analysisNames_Groups{2},overlayNum,'[18 18 21]');
        
        % save view
        mrSaveView(thisView)
    end
    
    %% Define GLM Gradient Reversal ROIs
    if doROIsGLM
        
        % create ROIs with the names:
        % LeftGR_GLM, LeftGRa_GLM, LeftGRp_GLM, RightGR_GLM, RightGRa_GLM, RightGRp_GLM based on gradient reversals, unsmoothed tonotopic maps and f-test maps
        
        % Also, line ROIs for each  reversal with the names:
        % LeftHighRevA_GLM, LeftLowRev_GLM, LeftHighRevP_GLM, RightHighRevA_GLM, RightLowRev_GLM, RightHighRevP_GLM
        
        % set view to what we need for ROI creation
        % set group
        groupName = [subjectInfo.flatmapNames{1} 'Volume'];
        if viewGet(thisView,'curgroup') ~= viewGet(thisView,'groupNum',groupName)
            thisView = viewSet(thisView,'curgroup',groupName);
        end
        % set analysis
        analysisName = 'combineTransformOverlay';
        if viewGet(thisView,'curAnalysis') ~= viewGet(thisView,'analysisNum',analysisName)
            thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',analysisName));
        end
        % set base
        flatmapName = [subjectInfo.flatmapNames{1} 'Volume'];
        if viewGet(thisView,'curbase') ~= viewGet(thisView,'basenum',flatmapName)
            thisView = viewSet(thisView,'curbase',viewGet(thisView,'basenum',flatmapName));
        end
        
        % set overlay
        overlayNum = viewGet(thisView,'overlayNum','Ouput 6 - gradientReversal_left_glm_hrfDoubleGamma(julien_pCF_nERB,[18 18 21])');
        thisView = viewSet(thisView,'curOverlay',overlayNum);
        
        refreshMLRDisplay(thisView.viewNum);
        
        disp('create ROIs with the names:')
        disp('LeftGRa_GLM, LeftGRp_GLM, , RightGRa_GLM, RightGRp_GLM based on gradient reversals, anatomy, unsmoothed tonotopic maps and f-test maps')
        disp('Also, line ROIs for each  reversal with the names: LeftHighRevA_GLM, LeftLowRev_GLM, LeftHighRevP_GLM, RightHighRevA_GLM, RightLowRev_GLM, RightHighRevP_GLM')
        disp('hit F5 when done')
        
        keyboard
        
        thisView = getMLRView;
        
        for iSide = 1:length(Info.Sides)
            newName = [Info.Sides{iSide} 'GR' '_GLM'];
            roi1 = [Info.Sides{iSide} 'GRa' '_GLM'];
            roi2 = [Info.Sides{iSide} 'GRp' '_GLM'];
            thisView = combineROIs(thisView,roi1,roi2,'union',newName);
        end
        
        disp('Check there are no holes in GR rois')
        disp('hit F5 when done')
        
        keyboard
        
        % get view in case we filled any holes
        thisView = getMLRView;
        
        % save view
        mrSaveView(thisView);
        
    end
    
    %% get GLM data
    % data has now be converted to flatmap space and averaged across cortical depth,
    % the following gets data from flatmap group, in flatmap space, from one depth (middle) that is the average across cortical depth.
    % create names to get data from overlays and save using structure side.Group.anal.data{iScan}
    if doGetDATA_GLM
        %% get data from SCANS
        analysisName = 'combineTransformOverlays';
        % voxelPropertyNames = {'Centriod','Spread','julien_pCF','julien_pTW','indexMax'};
        for iSide = 1:length(subjectInfo.flatmapNames)
            if viewGet(thisView,'curgroup') ~= viewGet(thisView,'groupNum',[subjectInfo.flatmapNames{iSide}, 'Volume'])
                thisView = viewSet(thisView,'curgroup',[subjectInfo.flatmapNames{iSide}, 'Volume']);
            end
            if viewGet(thisView,'curAnalysis') ~= viewGet(thisView,'analysisNum','combineTransformOverlays')
                thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum','combineTransformOverlays'));
            end
            groupName = [subjectInfo.flatmapNames{iSide} 'Volume'];
            thisView = viewSet(thisView,'curgroup',groupName);
            thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',analysisName));
            for iScan = 1:glmInfo.nScans
                for iAnal = 1:length(glmInfo.analysisNames_nCons)
                    % define overlay names
                    % First, deteremine how many stimuli
                    if glmInfo.analysisNStim{iAnal} == length(conditionNames{1})
                        overlayFlatNames = cell(1,length(conditionNames{1}));
                        conNamesString = [];
                        for iCon = 1:length(conditionNames{1})
                            overlayFlatNames{iCon} = ['averageDepthVol(Scan ' mat2str(iScan) ' - ' glmInfo.analysisNames_nCons{iAnal} '_Scan_' mat2str(iScan) ' (' conditionNames{1}{iCon} ',0))'];
                            if iCon == 1
                                conNamesString  = [conNamesString conditionNames{1}{iCon}];
                            else
                                conNamesString  = [conNamesString ',' conditionNames{1}{iCon}];
                            end
                        end
                    else
                        overlayFlatNames = cell(1,length(conditionNames{2}));
                        conNamesString = [];
                        for iCon = 1:length(conditionNames{2})
                            overlayFlatNames{iCon} = ['averageDepthVol(Scan ' mat2str(iScan) ' - ' glmInfo.analysisNames_nCons{iAnal} '_Scan_' mat2str(iScan)  ' (' conditionNames{2}{iCon} ',0))'];
                            if iCon == 1
                                conNamesString  = [conNamesString conditionNames{2}{iCon}];
                            else
                                conNamesString  = [conNamesString ',' conditionNames{2}{iCon}];
                            end
                        end
                    end
                    
                    % beta weights
                    overlayData = get_overlayData(thisView,overlayFlatNames);
                    eval(['data.' Info.Sides{iSide}, '.scan_',  num2str(iScan) '.', glmInfo.analysisNames_nCons{iAnal}, '.betas =  overlayData;']);
                    
                    % r2 overlay name
                    r2OverlayName = ['averageDepthVol(Scan ' mat2str(iScan) ' - ' glmInfo.analysisNames_nCons{iAnal} '_Scan_' mat2str(iScan) ' (r2,0))'];
                    
                    % voxel property estmate overlay names (change here if overlay names change)
                    voxelPropertyOverlayName = cell(1,length(glmInfo.voxelPropertyNames));
                    for iName = 1:length(glmInfo.voxelPropertyNames) - 1
                        voxelPropertyOverlayName{iName} = ['averageDepthVol(Scan ' mat2str(iScan) ' - ' glmInfo.analysisNames_nCons{iAnal} '_Scan_' mat2str(iScan) ' (' glmInfo.voxelPropertyNames{iName} '_nERB,0))'];
                    end
                    voxelPropertyOverlayName{end} = ['averageDepthVol(Scan ' mat2str(iScan) ' - ' glmInfo.analysisNames_nCons{iAnal} '_Scan_' mat2str(iScan) ' (' glmInfo.voxelPropertyNames{end} '_nERB,0))'];
                    
                    
                    % get overlay data using get_overlayData - outputs a structure with fields: .data & .name
                    % R2
                    clear tempData
                    tempData = get_overlayData(thisView,r2OverlayName);
                    eval(['data.' Info.Sides{iSide}, '.scan_', num2str(iScan), '.' , glmInfo.analysisNames_nCons{iAnal}, '.r2  = tempData;']);
                    
                    % voxel property estiamtes
                    for iName = 1:length(glmInfo.voxelPropertyNames)
                        clear tempData
                        tempData = get_overlayData(thisView,voxelPropertyOverlayName{iName});
                        eval(['data.' Info.Sides{iSide}, '.scan_', num2str(iScan), '.' , glmInfo.analysisNames_nCons{iAnal}, '.', glmInfo.voxelPropertyNames{iName},' = tempData;']);
                    end
                end
            end
        end
        
        %% get data from GROUPs
        % glmInfo.voxelPropertyNames = {'Centriod','Spread','julien_pCF','julien_pTW','indexMax'};
        for iSide = 1:length(subjectInfo.flatmapNames)
            if viewGet(thisView,'curgroup') ~= viewGet(thisView,'groupNum',[subjectInfo.flatmapNames{iSide}, 'Volume'])
                thisView = viewSet(thisView,'curgroup',[subjectInfo.flatmapNames{iSide}, 'Volume']);
            end
            if viewGet(thisView,'curAnalysis') ~= viewGet(thisView,'analysisNum','combineTransformOverlays')
                thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum','combineTransformOverlays'));
            end
            baseNum = viewGet(thisView,'baseNum',[subjectInfo.flatmapNames{iSide} 'Volume']);
            thisView = viewSet(thisView,'currentbase',baseNum);
            for iGroup = 1:length(glmInfo.groupNames)
                for iAnal = 1:length(glmInfo.analysisNames)
                    % define overlay names
                    % betas overlay names
                    overlayFlatNames = cell(1,length(conditionNames{1}));
                    conNamesString = [];
                    for iCon =1:length(conditionNames{1})
                        overlayFlatNames{iCon} = ['averageDepthVol(' glmInfo.groupNames{iGroup} '_' glmInfo.analysisNames{iAnal} ' (' conditionNames{1}{iCon} ',0))'];
                        if iCon == 1
                            conNamesString  = [conNamesString conditionNames{1}{iCon}];
                        else
                            conNamesString  = [conNamesString ',' conditionNames{1}{iCon}];
                        end
                    end
                    
                    % r2 overlay name
                    r2OverlayName = ['averageDepthVol(' glmInfo.groupNames{iGroup} '_' glmInfo.analysisNames{iAnal} ' (r2,0))'];
                    
                    % voxel property estmate overlay names (change here if overlay names change)
                    voxelPropertyOverlayName = cell(1,length(glmInfo.voxelPropertyNames));
                    for iName = 1:length(glmInfo.voxelPropertyNames) - 1
                        voxelPropertyOverlayName{iName} = ['averageDepthVol(' glmInfo.groupNames{iGroup} '_' glmInfo.analysisNames{iAnal} ' (' glmInfo.voxelPropertyNames{iName} '_nERB,0))'];
                    end
                    voxelPropertyOverlayName{end} = ['averageDepthVol(' glmInfo.groupNames{iGroup} '_' glmInfo.analysisNames{iAnal} ' (' glmInfo.voxelPropertyNames{end} '_nERB,0))'];
                    
                    
                    % get overlay data using get_overlayData - outputs a structure with fields: .data & .name
                    % R2
                    clear tempData
                    tempData = get_overlayData(thisView,r2OverlayName);
                    eval(['data.' Info.Sides{iSide}, '.' glmInfo.groupNames{iGroup} '.', glmInfo.analysisNames{iAnal}, '.r2 = tempData;']);
                    
                    % beta weights
                    clear tempData
                    tempData = get_overlayData(thisView,overlayFlatNames);
                    eval(['data.' Info.Sides{iSide}, '.' glmInfo.groupNames{iGroup} '.', glmInfo.analysisNames{iAnal}, '.betas = tempData;']);
                    
                    % voxel property estiamtes
                    for iName = 1:length(glmInfo.voxelPropertyNames)
                        clear tempData
                        tempData = get_overlayData(thisView,voxelPropertyOverlayName{iName});
                        eval(['data.' Info.Sides{iSide}, '.' glmInfo.groupNames{iGroup} '.', glmInfo.analysisNames{iAnal}, '.', glmInfo.voxelPropertyNames{iName},' = tempData;']);
                    end
                    
                end
            end
        end
    end
    
    %% Restrict data by ROIs
    % Now we need to restrict the data by the ROIs
    if doROIRestrict_GLM
        
        % first get view so we have the ROIS
        thisView = getMLRView;
        
        q = char(39);
        for iSide = 1:length(Info.Sides)
            if viewGet(thisView,'curgroup') ~= viewGet(thisView,'groupNum',[subjectInfo.flatmapNames{iSide}, 'Volume'])
                thisView = viewSet(thisView,'curgroup',[subjectInfo.flatmapNames{iSide}, 'Volume']);
            end
            if viewGet(thisView,'curAnalysis') ~= viewGet(thisView,'analysisNum','combineTransformOverlays')
                thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum','combineTransformOverlays'));
            end
            eval(['ROInames = Info.' Info.Sides{iSide} 'ROInames;']);
            baseNum = viewGet(thisView,'baseNum',[subjectInfo.flatmapNames{iSide} 'Volume']);
            thisView = viewSet(thisView,'currentbase',baseNum);
            
            for iROI = 1:length(ROInames)
                
                % get ROI struct (overwrite old roi definitions)
                eval(['data.' ROInames{iROI} ' = struct;']);
                eval(['data.' ROInames{iROI} '.roi = viewGet(thisView,' q 'roi' q ',ROInames{iROI});']);
                
                % Restrict Group data
                for iGroup = 1:length(glmInfo.groupNames)
                    for iAnal = 1:length(glmInfo.analysisNames_Groups)/length(glmInfo.groupNames)
                        analysisName = glmInfo.analysisNames_Groups{iAnal};
                        
                        % restrict r2
                        eval(['groupDataVar = data.' Info.Sides{iSide} '.' glmInfo.groupNames{iGroup} '.' glmInfo.analysisNames_Groups{iAnal} '.r2.data;']);
                        eval(['data.'  ROInames{iROI} '.' glmInfo.groupNames{iGroup} '.' glmInfo.analysisNames_Groups{iAnal} '.r2  = get_ROIdata(groupDataVar,data.' ROInames{iROI} '.roi);']);
                        clear groupDataVar
                        
                        % restrict betas
                        eval(['groupDataVar = data.' Info.Sides{iSide} '.' glmInfo.groupNames{iGroup} '.' glmInfo.analysisNames_Groups{iAnal} '.betas.data;']);
                        eval(['data.' ROInames{iROI} '.' glmInfo.groupNames{iGroup} '.' glmInfo.analysisNames_Groups{iAnal} '.betas  = get_ROIdata(groupDataVar,data.' ROInames{iROI} '.roi);']);
                        clear groupDataVar
                        
                        % restrict tonotopic estimates
                        for iName = 1:length(glmInfo.voxelPropertyNames)
                            eval(['groupDataVar = data.' Info.Sides{iSide} '.' glmInfo.groupNames{iGroup} '.' glmInfo.analysisNames_Groups{iAnal} '.' glmInfo.voxelPropertyNames{iName} '.data;']);
                            eval(['data.' ROInames{iROI} '.' glmInfo.groupNames{iGroup} '.' glmInfo.analysisNames_Groups{iAnal} '.' glmInfo.voxelPropertyNames{iName} '  = get_ROIdata(groupDataVar,data.' ROInames{iROI} '.roi);']);
                            clear groupDataVar
                        end
                        
                    end
                end
                
                % Restrict Scan data
                for iAnal = 1:length(glmInfo.analysisBaseNames_Scans)/glmInfo.nScans
                    for iScan = 1:glmInfo.nScans
                        
                        % restrict r2
                        tempData = [];
                        eval(['scanDataVar = data.' Info.Sides{iSide} '.scan_' num2str(iScan) '.' glmInfo.analysisBaseNames_Scans{iAnal} '.r2.data;']);
                        eval(['tempData = get_ROIdata(scanDataVar,data.' ROInames{iROI} '.roi);']);
                        eval(['data.' ROInames{iROI} '.scan_'  num2str(iScan) '.' glmInfo.analysisBaseNames_Scans{iAnal} '.r2 = tempData{:};']);
                        clear scanDataVar
                        
                        eval(['scanDataVar = data.' Info.Sides{iSide} '.scan_' num2str(iScan) '.' glmInfo.analysisBaseNames_Scans{iAnal} '.betas.data;']);
                        
                        eval(['data.' ROInames{iROI} '.scan_' num2str(iScan) '.' glmInfo.analysisBaseNames_Scans{iAnal} '.betas = get_ROIdata(scanDataVar,data.' ROInames{iROI} '.roi);']);
                        clear scanDataVar
                        
                        % restrict tonotopic estimates
                        for iName = 1:length(glmInfo.voxelPropertyNames)
                            
                            tempData = [];
                            eval(['scanDataVar = data.' Info.Sides{iSide} '.scan_' num2str(iScan) '.' glmInfo.analysisBaseNames_Scans{iAnal} '.' glmInfo.voxelPropertyNames{iName} '.data;']);
                            eval(['tempData = get_ROIdata(scanDataVar,data.' ROInames{iROI} '.roi);']);
                            eval(['data.' ROInames{iROI} '.scan_' num2str(iScan) '.' glmInfo.analysisBaseNames_Scans{iAnal} '.' glmInfo.voxelPropertyNames{iName} ' = tempData{:};']);
                            clear scanDataVar
                        end
                        
                    end
                end
            end
        end
        
        
        %% Check data / GLM ROI analysis
        % check data using GLM ROI analysis
        % quick and dirty plots to check the data
        for iSide = 1:length(Info.Sides)
            eval(['ROInames = Info.' Info.Sides{iSide} 'ROInames;']);
            %% %%%%%%%%%% comment what these are %%%%%%%%%%%%%
            iROI = 1;
            iGroup = 1;
            iAnal = 2;
            iName = 3;
            
            eval(['roipCFdataA = data.' ROInames{iROI} '.' glmInfo.groupNames{1} '.' glmInfo.analysisNames_Groups{iAnal} '.' glmInfo.voxelPropertyNames{iName} ';']);
            eval(['roipCFdataB = data.' ROInames{iROI} '.' glmInfo.groupNames{2} '.' glmInfo.analysisNames_Groups{iAnal} '.' glmInfo.voxelPropertyNames{iName} ';']);
            
            figure
            subplot(2,1,1)
            histogram(cell2mat(roipCFdataA))
            hold on
            histogram(cell2mat(roipCFdataB))
            
            eval(['roiBetadataA = data.' ROInames{iROI} '.' glmInfo.groupNames{1} '.' glmInfo.analysisNames_Groups{iAnal} '.betas;']);
            eval(['roiBetadataB = data.' ROInames{iROI} '.' glmInfo.groupNames{2} '.' glmInfo.analysisNames_Groups{iAnal} '.betas;']);
            
            betas_mv_A = cal_movingAverage(cell2mat(roiBetadataA'));
            betas_mv_B = cal_movingAverage(cell2mat(roiBetadataB'));
            
            subplot(2,1,2)
            plot(mean(betas_mv_A,2))
            hold on
            plot(mean(betas_mv_B,2))
            
            for iScan = 1:4
                eval(['roiSplitBetas{iScan} = data.' ROInames{iROI} '.scan_' num2str(iScan) '.' glmInfo.analysisBaseNames_Scans{4} '.betas;']);
                splitBetas{iScan} = cell2mat(roiSplitBetas{iScan}');
            end
            [splitMeanA, ROI_data, Voxel_data, totalROIpCF] = cal_splitMean(splitBetas{1},splitBetas{3});
            [splitMeanB, ROI_data, Voxel_data, totalROIpCF] = cal_splitMean(splitBetas{2},splitBetas{4});
            
            figure
            for i =1:8
                subplot(2,4,i)
                plot(splitMeanA(i,:))
                hold on
                plot(splitMeanB(i,:))
            end
            
        end
        
        %% save data
        % save so don't need to load mrView again
        % this data structure is loaded by CM_groupAnalysisScript for group analysis
        save(saveName,'data','-v7.3');
        
        %% save view
        mrSaveView(thisView);
        
    end
    
    %% pRF Analysis %%
    % voxelwise analysis using popualtion receptive field modelling
    
    if dopRF
        % first get view so we have the ROIs
        thisView = getMLRView;
        
        %% pRF analysis
        % perform pRF analysis (restricted to auditory responsive voxels * [3 3 3] sphere (ARexp ROI))
        % make sure analysis roi is resctricted (cmd + x) to glm overlay data (so we only perform analysis for voxels with functional data)
        pRFrestrictROI = 'ARexp';
        pRFanalysisName = ['pRF_', pRFrestrictROI];
        
        % use HRF params from GLM deconvolution estimate
        pRFInfo.hrfParamsGamma = data.hrf.x_Gamma;
        pRFInfo.hrfParamsDiffofGamma = data.hrf.x_dGamma;
        
        [thisView, pRFParams] = script_pRFAnalysis(thisView,pRFInfo,glmInfo,pRFrestrictROI,1,0);
    end
    
    %% pRF grandient reversals
    if doGradientReversal_pRF
        % thisView = script_flatMapAnalysis(thisView,Info,subjectInfo,groupBase,analysisBase,overlayNumber,smoothingParams)
        thisView = script_flatMapAnalysis(thisView,Info,subjectInfo,Info.gradReversalInfo.groupBase,pRFanalysisName,pRFInfo.pRFgradientReversalOverlay,'[18 18 21]');
    end
    
    %% Define pRF Gradient Reversal ROIs
    if doROIspRF
        
        % create ROIs with the names:
        % LeftGR_pRF, LeftGRa_pRF, LeftGRp_pRF, RightGR_pRF,
        % RightGRa_pRF, RightGRp_pRF based on gradient reversals, unsmoothed tonotopic maps and f-test maps.
        
        % Also, line ROIs for each  reversal with the names:
        % LeftHighRevA_pRF, LeftLowRev_pRF, LeftHighRevP_pRF, RightHighRevA_pRF, RightLowRev_pRF, RightHighRevP_pRF
        
        % set group
        groupName = [subjectInfo.flatmapNames{1} 'Volume']
        if viewGet(thisView,'curgroup') ~= viewGet(thisView,'groupNum',groupName)
            thisView = viewSet(thisView,'curgroup',groupName);
        end
        % set analysis
        analysisName = 'combineTransformOverlay';
        if viewGet(thisView,'curAnalysis') ~= viewGet(thisView,'analysisNum',analysisName)
            thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',analysisName));
        end
        % set base
        flatmapName = [subjectInfo.flatmapNames{1} 'Volume'];
        if viewGet(thisView,'curbase') ~= viewGet(thisView,'basenum',flatmapName)
            thisView = viewSet(thisView,'curbase',viewGet(thisView,'basenum',flatmapName));
        end
        
        % set overlay
        overlayNum = viewGet(thisView,'overlayNum','Ouput 6 - gradientReversal_left_pRF(pCF,[18 18 21])');
        thisView = viewSet(thisView,'curOverlay',overlayNum);
        
        disp('create ROIs with the names:')
        disp('LeftGRa_pRF, LeftGRp_pRF, , RightGRa_pRF, RightGRp_pRF based on gradient reversals, unsmoothed tonotopic maps and f-test maps')
        disp('Also, line ROIs for each reversal with the names:')
        disp('LeftHighRevA_pRF, LeftLowRev_pRF, LeftHighRevP_pRF, RightHighRevA_pRF, RightLowRev_pRF, RightHighRevP_pRF')
        disp('hit F5 when done')
        
        keyboard
        
        thisView = getMLRView;
        
        for iSide = 1:length(Info.Sides)
            newName = [Info.Sides{iSide} 'GR' '_pRF'];
            roi1 = [Info.Sides{iSide} 'GRa' '_pRF'];
            roi2 = [Info.Sides{iSide} 'GRp' '_pRF'];
            thisView = combineROIs(thisView,roi1,roi2,'union',newName);
        end
        
        disp('Check there are no holes in GR rois')
        disp('hit F5 when done')
        
        keyboard
        
        % get view in case we filled any holes
        thisView = getMLRView;
        
        % save view
        mrSaveView(thisView);
        
    end
    
    %% Convert pRF data to flatmap space then average over depth cortical depth
    % pRFOverlayNames = {'r2','PrefCentreFreq','rfHalfWidth'};
    if doConvertvol2FlatAvDepth_pRF
        for iSide = 1:length(subjectInfo.flatmapNames)
            % Groups
            if viewGet(thisView,'curgroup') ~= viewGet(thisView,'groupNum',[subjectInfo.flatmapNames{iSide}, 'Volume'])
                thisView = viewSet(thisView,'curgroup',[subjectInfo.flatmapNames{iSide}, 'Volume']);
            end
            if viewGet(thisView,'curAnalysis') ~= viewGet(thisView,'analysisNum','combineTransformOverlays')
                thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum','combineTransformOverlays'));
            end
            
            for iGroup = 1:length(glmInfo.groupNames)
                groupName = glmInfo.groupNames{iGroup};
                
                for iAnal = 1:length(pRFInfo.analysisNames_Groups{iGroup})
                    pRFanalysisName = [pRFInfo.analysisNames_Groups{iGroup}{iAnal}, '_',  pRFrestrictROI];
                    % get overlay names and numbers
                    overlayFlatNames = cell(1,length(pRFInfo.pRFOverlayNames));
                    overlayNum = zeros(1,length(pRFInfo.pRFOverlayNames));
                    overlay2Get = cell(1,length(pRFInfo.pRFOverlayNames));
                    for iOverlay = 1:length(pRFInfo.pRFOverlayNames)
                        overlayFlatNames{iOverlay} = [groupName '_' pRFanalysisName ' (' pRFInfo.pRFOverlayNames{iOverlay} ',0)'];
                        overlayNum(iOverlay) = iOverlay;
                        overlay2Get{iOverlay} = ['averageDepthVol(' overlayFlatNames{iOverlay} ')'];
                    end
                    
                    % export group data from volumetric to flatmap space
                    % [thisView, analysisData] = script_covertData2FlatmapSpace(thisView,groupName,analysisName,iScan,overlays,flatmapName)
                    thisView = script_covertData2FlatmapSpace(thisView,glmInfo.groupNames{iGroup},pRFanalysisName,[],overlayNum,subjectInfo.flatmapNames{iSide});
                    
                    % average over cortical depth
                    thisView = script_averageAcrossDepths(thisView,overlayFlatNames,[subjectInfo.flatmapNames{iSide}, 'Volume'],1);
                end
            end
        end
        for iSide = 1:length(subjectInfo.flatmapNames)
            % Scans
            for iScan = 1:glmInfo.nScans
                if viewGet(thisView,'curgroup') ~= viewGet(thisView,'groupNum',[subjectInfo.flatmapNames{iSide}, 'Volume'])
                    thisView = viewSet(thisView,'curgroup',[subjectInfo.flatmapNames{iSide}, 'Volume']);
                end
                if viewGet(thisView,'curAnalysis') ~= viewGet(thisView,'analysisNum','combineTransformOverlays')
                    thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum','combineTransformOverlays'));
                end
                pRFanalysisName = ['pRF_',  pRFrestrictROI, '_Scan_' num2str(iScan)];
                overlayNames = cell(1,length(pRFInfo.pRFOverlayNames));
                overlayFlatNames = cell(1,length(pRFInfo.pRFOverlayNames));
                overlayNum = zeros(1,length(pRFInfo.pRFOverlayNames));
                overlay2Get = cell(1,length(pRFInfo.pRFOverlayNames));
                for iOverlay = 1:length(pRFInfo.pRFOverlayNames)
                    
                    overlayNames{iOverlay} = [pRFanalysisName ' (' pRFInfo.pRFOverlayNames{iOverlay} ',0)'];
                    overlayFlatNames{iOverlay} = ['Scan ' num2str(iScan) ' - '  pRFanalysisName ' (' pRFInfo.pRFOverlayNames{iOverlay} ',0)'];
                    overlayNum(iOverlay) = iOverlay;
                    overlay2Get{iOverlay} = ['averageDepthVol(' overlayFlatNames{iOverlay} ')'];
                end
                
                % export scan data from volumetric to flatmap space
                % [thisView, analysisData] = script_covertData2FlatmapSpace(thisView,groupName,analysisName,iScan,overlays,flatmapName)
                thisView = script_covertData2FlatmapSpace(thisView,glmInfo.scanGroupName,pRFanalysisName,iScan,overlayNum,subjectInfo.flatmapNames{iSide});
                
                % average over cortical depth
                thisView = script_averageAcrossDepths(thisView,overlayFlatNames,[subjectInfo.flatmapNames{iSide}, 'Volume'],1);
                
            end
        end
        
    end
    
    %% Get pRF data and restrict by ROIs
    if doGetAndRestrictDATA_pRF
        
        % get group data and restrict by ROIs
        for iSide = 1:length(Info.Sides)
            if viewGet(thisView,'curgroup') ~= viewGet(thisView,'groupNum',[subjectInfo.flatmapNames{iSide}, 'Volume'])
                thisView = viewSet(thisView,'curgroup',[subjectInfo.flatmapNames{iSide}, 'Volume']);
            end
            if viewGet(thisView,'curAnalysis') ~= viewGet(thisView,'analysisNum','combineTransformOverlays')
                thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum','combineTransformOverlays'));
            end
            eval(['roiNames = Info.' Info.Sides{iSide} 'ROInames;']);
            baseNum = viewGet(thisView,'baseNum',[subjectInfo.flatmapNames{iSide} 'Volume']);
            thisView = viewSet(thisView,'currentbase',baseNum);
            
            for iGroup = 1:length(glmInfo.groupNames)
                groupName = glmInfo.groupNames{iGroup};
                for iAnal = 1:length(pRFInfo.analysisNames_Groups{iGroup})
                    pRFanalysisName = [pRFInfo.analysisNames_Groups{iGroup}{iAnal}, '_',  pRFrestrictROI];
                    overlayFlatNames = cell(1,length(pRFInfo.pRFOverlayNames));
                    overlay2Get = cell(1,length(pRFInfo.pRFOverlayNames));
                    for iOverlay = 1:length(pRFInfo.pRFOverlayNames)
                        % get overlay names
                        overlayFlatNames{iOverlay} = [groupName '_' pRFanalysisName ' (' pRFInfo.pRFOverlayNames{iOverlay} ',0)'];
                        overlay2Get{iOverlay} = ['averageDepthVol(' overlayFlatNames{iOverlay} ')'];
                        clear tempData
                        tempData = get_overlayData(thisView,overlay2Get{iOverlay});
                        % get data
                        eval(['data.' Info.Sides{iSide}, '.', glmInfo.groupNames{iGroup}, '.', pRFInfo.analysisNames_Groups{iGroup}{iAnal}, '.', pRFInfo.pRFOverlayNames{iOverlay}, ' = tempData;']);
                        
                        % restrict by ROI
                        for iROI = 1:length(roiNames)
                            eval(['roi = data.' Info.Sides{iSide} '.' roiNames{iROI} '.roi;']);
                            eval(['data.', roiNames{iROI}, '.', glmInfo.groupNames{iGroup}, '.', pRFInfo.analysisNames_Groups{iGroup}{iAnal}, '.', pRFInfo.pRFOverlayNames{iOverlay}, '  = get_ROIdata(tempData.data,roi);']);
                            
                        end
                        
                    end
                end
            end
        end
        
        % get scan data and restrict by ROIs
        for iSide = 1:length(Info.Sides)
            if viewGet(thisView,'curgroup') ~= viewGet(thisView,'groupNum',[subjectInfo.flatmapNames{iSide}, 'Volume'])
                thisView = viewSet(thisView,'curgroup',[subjectInfo.flatmapNames{iSide}, 'Volume']);
            end
            if viewGet(thisView,'curAnalysis') ~= viewGet(thisView,'analysisNum','combineTransformOverlays')
                thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum','combineTransformOverlays'));
            end
            eval(['roiNames = Info.' Info.Sides{iSide} 'ROInames;']);
            baseNum = viewGet(thisView,'baseNum',[subjectInfo.flatmapNames{iSide} 'Volume']);
            thisView = viewSet(thisView,'currentbase',baseNum);
            
            for iScan = 1:glmInfo.nScans
                for iAnal = 1:length(pRFInfo.analysisNames_Groups{iGroup})
                    pRFanalysisName = ['pRF_',  pRFrestrictROI, '_Scan_' num2str(iScan)];
                    overlayFlatNames = cell(1,length(pRFInfo.pRFOverlayNames));
                    overlay2Get = cell(1,length(pRFInfo.pRFOverlayNames));
                    for iOverlay = 1:length(pRFInfo.pRFOverlayNames)
                        % get overlay names
                        overlayFlatNames{iOverlay} = ['Scan ' num2str(iScan) ' - '  pRFanalysisName ' (' pRFInfo.pRFOverlayNames{iOverlay} ',0)'];
                        overlay2Get{iOverlay} = ['averageDepthVol(' overlayFlatNames{iOverlay} ')'];
                        clear tempData
                        tempData = get_overlayData(thisView,overlay2Get{iOverlay});
                        % get data
                        eval(['data.', Info.Sides{iSide}, '.scan_', num2str(iScan), '.', pRFInfo.analysisNames_Groups{iGroup}{iAnal}, '.', pRFInfo.pRFOverlayNames{iOverlay}, ' = tempData;']);
                        
                        % restrict by ROI
                        for iROI = 1:length(roiNames)
                            clear tempROIdata
                            eval(['roi = data.' Info.Sides{iSide} '.' roiNames{iROI} '.roi;']);
                            %  eval(['data.', Info.Sides{iSide}, '.', roiNames{iROI}, '.scanData.', pRFInfo.analysisNames_Groups{iGroup}{iAnal}, '.', pRFInfo.pRFOverlayNames{iOverlay}, '{iScan}  = get_ROIdata(tempData.data,roi);']);
                            eval(['tempROIdata = get_ROIdata(tempData.data,roi);']);
                            
                            eval(['data.', roiNames{iROI}, '.scan_', num2str(iScan), '.', pRFInfo.analysisNames_Groups{iGroup}{iAnal}, '.', pRFInfo.pRFOverlayNames{iOverlay}, '  = tempROIdata{:};'])
                        end
                        
                    end
                end
            end
        end
        
        
        %% pRF analysis
        % create function or add to: voxel comparisions;
        % pCF scatter, pCF distribution, pCF correlation
        % r = corr2(A,B)
        
        %% check data
        % quick plots to make sure its worked
        for iSide = 1:length(Info.Sides)
            eval(['ROInames = Info.' Info.Sides{iSide} 'ROInames;']);
            iROI = 1;
            iGroup = 1;
            iAnal = 1;
            iOverlay = 2;
            
            pRFanalysisName = [pRFInfo.analysisNames_Groups{iGroup}{iAnal}, '_',  pRFrestrictROI];
            
            eval(['roipCFdataA = data.' ROInames{iROI} '.' glmInfo.groupNames{1} '.' pRFInfo.analysisNames_Groups{iGroup}{iAnal} '.' pRFInfo.pRFOverlayNames{iOverlay} ';']);
            eval(['roipCFdataB = data.' ROInames{iROI} '.' glmInfo.groupNames{2} '.' pRFInfo.analysisNames_Groups{iGroup}{iAnal} '.' pRFInfo.pRFOverlayNames{iOverlay} ';']);
            
            figure
            histogram(cell2mat(roipCFdataA))
            hold on
            histogram(cell2mat(roipCFdataB))
            
            for iScan = 1:4
                eval(['roiSplitpCF{iScan} = data.' ROInames{iROI} '.scan_' num2str(iScan) '.' pRFInfo.analysisNames_Groups{iGroup}{iAnal} '.' pRFInfo.pRFOverlayNames{iOverlay} ';']);
            end
            figure
            subplot(2,1,1)
            scatter(roiSplitpCF{1},roiSplitpCF{3})
            subplot(2,1,2)
            scatter(roiSplitpCF{2},roiSplitpCF{4})
            
            % need to remove nans before - also findout how many
            corrcoef(double(roiSplitpCF{1}'),double(roiSplitpCF{3}'))
            corr(roiSplitpCF{1}',roiSplitpCF{3}')
        end
        
        
        % save data
        save(saveName,'data','-v7.3');
        
        % save view
        mrSaveView(thisView);
        
    end
    
    %% Cortical Magnification %%
    if doCorticalMagnification
        % get cortical distances
        % get curvature?
        % get pCF estiamte
        % plot to check
        
        % Define gradient reversal ROIs
        % 	High frequency - line
        % 	Low frequency - line
        % 	Gradient - polygon
        
        % loop over:
        %   sides
        %   rois / anterior/ posterior
        %   group
        %   anlysis
        
        % ROIs could have already been create, if not:
        % Create ROIs for each gradient reversal and name:
        %   'LeftHa','LeftHp','LeftLow,'LeftGRa''LeftGRp'
        %   'RightHa','RightHp','RightLow','RightGRa''RightGRp'
        %   Post fix with _GLM or _pRF
        
        thisView = getMLRView;
        
        % convert rois from flatmap to volume space
        pRFrestrictROI = 'ARexp';
        pRFanalysisName = ['pRF_', pRFrestrictROI];
        analysisNames = {'glm_hrfDoubleGamma',pRFanalysisName};
        AP = {'a','p'};
        %         analName = {'GLM', 'pRF'};
        analName = {'GLM', 'GLM'};
        
        if doGRrois
            for iSide = 1:length(Info.Sides)
                
                % define gradient reversals ROIs
                
                % set group
                groupName = [subjectInfo.flatmapNames{iSide} 'Volume'];
                if viewGet(thisView,'curgroup') ~= viewGet(thisView,'groupNum',groupName)
                    thisView = viewSet(thisView,'curgroup',groupName);
                end
                % set analysis
                analysisName = 'combineTransformOverlay';
                if viewGet(thisView,'curAnalysis') ~= viewGet(thisView,'analysisNum',analysisName)
                    thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',analysisName));
                end
                % set base
                flatmapName = [subjectInfo.flatmapNames{iSide} 'Volume'];
                if viewGet(thisView,'curbase') ~= viewGet(thisView,'basenum',flatmapName)
                    thisView = viewSet(thisView,'curbase',viewGet(thisView,'basenum',flatmapName));
                end
                
                % set overlay
                overlayNum = viewGet(thisView,'overlayNum',['Ouput 6 - gradientReversal_' Info.sides{iSide} '_glm_hrfDoubleGamma(julien_pCF_nERB,[18 18 21])']);
                thisView = viewSet(thisView,'curOverlay',overlayNum);
                
                disp('create line ROIs for each reversal with the names:')
                disp('LeftHighRevA_pRF, LeftLowRev_pRF, LeftHighRevP_pRF, RightHighRevA_pRF, RightLowRev_pRF, RightHighRevP_pRF')
                disp('hit F5 when done')
                
                keyboard
                
                % get view
                thisView = getMLRView;
                
                % make sure we have the correct group and analysis
                if viewGet(thisView,'curgroup') ~= viewGet(thisView,'groupNum',[subjectInfo.flatmapNames{iSide}, 'Volume'])
                    thisView = viewSet(thisView,'curgroup',[subjectInfo.flatmapNames{iSide}, 'Volume']);
                end
                if viewGet(thisView,'curAnalysis') ~= viewGet(thisView,'analysisNum','combineTransformOverlays')
                    thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum','combineTransformOverlays'));
                end
                
                for iAnal = 1:length(analysisNames)
                    for iAP = 1:length(AP)
                        clear rois
                        
                        % define ROI names then get ROI numbers then put in a list
                        HighROIname = [Info.Sides{iSide}, 'High' AP{iAP} '_' analName{iAnal}];
                        roiList(1) = viewGet(thisView,'roinum',HighROIname);
                        
                        LowROIname = [Info.Sides{iSide}, 'Low', '_' analName{iAnal}];
                        roiList(2) = viewGet(thisView,'roinum',LowROIname);
                        
                        GRROIname = [Info.Sides{iSide}, 'GR' AP{iAP} '_' analName{iAnal}];
                        roiList(3) = viewGet(thisView,'roinum',GRROIname);
                        
                        % get rois
                        rois = viewGet(thisView,'roi',roiList);
                        
                        % Project ROIs across all cortical depths
                        for iROI = 1:3
                            clear outputRoi
                            outputRoi = expandROI(rois(iROI),[1 1 6],[]);
                            outputRoi.name = [outputRoi.name '_ex'];
                            thisView = viewSet(thisView,'newROI',outputRoi);
                        end
                        
                        clear roiList
                        
                        % Define expanded ROI names and get numbers then put in a list
                        HighROIname = [Info.Sides{iSide}, 'High' AP{iAP} '_' analName{iAnal} '_ex'];
                        roiList(1) = viewGet(thisView,'roinum',HighROIname);
                        
                        LowROIname = [Info.Sides{iSide}, 'Low', '_' analName{iAnal} '_ex'];
                        roiList(2) = viewGet(thisView,'roinum',LowROIname);
                        
                        GRROIname = [Info.Sides{iSide}, 'GR' AP{iAP} '_' analName{iAnal} '_ex'];
                        roiList(3) = viewGet(thisView,'roinum',GRROIname);
                        
                        clear rois_ex
                        
                        % get expanded rois
                        rois_ex = viewGet(thisView,'roi',roiList);
                        
                        % convert from flatmap space to volumetric
                        for iROI = 1:3
                            clear outputRoi
                            outputRoi = convertFromFlatVolumeToBase(rois_ex(iROI));
                            outputRoi.name = [outputRoi.name '_Vol'];
                            thisView = viewSet(thisView,'newROI',outputRoi);
                        end
                        
                        % delete expanded rois we don't need
                        thisView = viewSet(thisView,'deleteroi',roiList);
                    end
                end
            end
        end
        
        % set view to FNIRTed flatmaps
        % change view to the group / analysis / overlay of interest **** loop over this ****
        % get undistorted flatmap names and send to analysis
        % Perform AverageDepthVol on overlay of interest - calculate in flat space but export to base space
        % Now run cal_tonotopicMagnification to get data
        
        % overlays names to get for GLM and pRF analys
        glmoverlayNames = {'r2', [glmInfo.voxelPropertyNames{3} '_nERB'], [glmInfo.voxelPropertyNames{4} '_nERB']};
        pRFoverlayNames = {pRFInfo.pRFOverlayNames{1}, pRFInfo.pRFOverlayNames{2}, pRFInfo.pRFOverlayNames{3}};
        % names used in saved data structure
        estimateNames = {'r2', 'pCF', 'pTW'};
        
        for iSide = 1:length(Info.Sides)
            % non- fnirted flatmap names
            flatmapnames.gm = [subjectInfo.freeSurferName '_' Info.sides{iSide} '_GM.off'];
            flatmapnames.wm = [subjectInfo.freeSurferName '_' Info.sides{iSide} '_WM.off'];
            
            % set base
            
            if viewGet(thisView,'currentbase') ~= viewGet(thisView,'baseNum',[subjectInfo.flatmapNames{iSide}])
                baseNum = viewGet(thisView,'baseNum',[subjectInfo.flatmapNames{iSide}]);
                thisView = viewSet(thisView,'currentbase',baseNum);
            end
            
            % set path to match the current directory
            % its different for different computers..
            %             '/Volumes/data_PSY/data/Anatomy/freesurfer/subjects/11108_006/surfRelax'
            baseParams = viewGet(thisView,'basecoordmap');
            baseParams.path = fullfile(Info.dataDir,'Anatomy/freesurfer/subjects/',subjectInfo.freeSurferName,'surfRelax');
            thisView = viewSet(thisView,'basecoordmap',baseParams);
            
            for iGroup = 1:length(glmInfo.groupNames) % change change this if we want to calculate for both groups or not
                % select group
                groupName = glmInfo.groupNames{iGroup};
                if viewGet(thisView,'curgroup') ~= viewGet(thisView,'groupNum',groupName)
                    thisView = viewSet(thisView,'curgroup',groupName);
                end
                % Select analysis
                for iAnal = 1:length(analysisNames)
                    analysisName = analysisNames{iAnal};
                    
                    if viewGet(thisView,'curAnalysis') ~= viewGet(thisView,'analysisNum',analysisName)
                        thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',analysisName));
                    end
                    
                    for iOverlay = 1:3
                        % select overlay
                        clear overlayNum
                        switch analysisName
                            case 'glm_hrfDoubleGamma'
                                overlayNum = viewGet(thisView,'overlayNum',glmoverlayNames{iOverlay});
                                %                                 overlayNames = glmoverlayNames;
                            case pRFanalysisName
                                overlayNum = viewGet(thisView,'overlayNum',pRFoverlayNames{iOverlay});
                                %                                 overlayNames = pRFoverlayNames;
                        end
                        
                        % check if overlay has already been averaged
                        clear overlayCheck
                        clear overlayname
                        clear tempOverlayData tempOverlayRoiData overlay2getName
                        switch analysisName
                            case 'glm_hrfDoubleGamma'
                                overlay2check = ['averageDepthVol_' Info.Sides{iSide} '(', glmoverlayNames{iOverlay} ')'];
                                
                            case pRFanalysisName
                                overlay2check = ['averageDepthVol_' Info.Sides{iSide} '(', pRFoverlayNames{iOverlay}, ')'];
                        end
                        
                        overlayCheck = viewGet(thisView,'overlayNum',overlay2check);
                        
                        % average overlays over cortical depths
                        if isempty(overlayCheck)
                            [thisView,params] = combineTransformOverlays(thisView,[],'justGetParams=1','defaultParams=1',['overlayList=' mat2str(overlayNum)]);
                            params.combineFunction='averageDepthVol';
                            params.combinationMode = 'Apply function to each overlay';
                            params.baseSpace = 1; % perfrom analysis in base space (flatmap)
                            params.exportToNewGroup = 0; % don't export to volume in a new group
                            params.outputName = ['averageDepthVol_', Info.Sides{iSide}];
                            [thisView,params] = combineTransformOverlays(thisView,params);
                        end
                    end
                    
                    switch analysisName
                        case 'glm_hrfDoubleGamma'
                            pCFoverlayname = ['averageDepthVol_' Info.Sides{iSide} '(', glmoverlayNames{2} ')'];
                            
                        case pRFanalysisName
                            pCFoverlayname = ['averageDepthVol_' Info.Sides{iSide} '(', pRFoverlayNames{2}, ')'];
                    end
                    
                    AvOverlayNum = viewGet(thisView,'overlayNum',pCFoverlayname);
                    thisView = viewSet(thisView,'curOverlay',AvOverlayNum);
                    
                    for iAP = 1:length(AP)
                        % make ROI names and get numbers then put in a list
                        clear roiList
                        HighROIname = [Info.Sides{iSide}, 'High' AP{iAP} '_' analName{iAnal} '_ex_Vol'];
                        roiList(1) = viewGet(thisView,'roinum',HighROIname);
                        
                        LowROIname = [Info.Sides{iSide}, 'Low', '_' analName{iAnal} '_ex_Vol'];
                        roiList(2) = viewGet(thisView,'roinum',LowROIname);
                        
                        GRROIname = [Info.Sides{iSide}, 'GR' AP{iAP} '_' analName{iAnal} '_ex_Vol'];
                        roiList(3) = viewGet(thisView,'roinum',GRROIname);
                        
                        % calculate Tonotopic Magnification
                        clear relativeDistances overlayRoiData pathDistances
                        [relativeDistances, overlayRoiData, pathDistances, verticesScanCoords] = cal_tonotopicMagnification(thisView,roiList,flatmapnames,[]);
                        
                        % save to structure
                        roiSaveName = [Info.Sides{iSide}, 'GR' AP{iAP} '_' analName{iAnal}];
                        eval(['data.' roiSaveName '.' groupName '.' analysisName '.tonotopicMagnificaion.relativeDistances = relativeDistances;']);
                        eval(['data.' roiSaveName '.' groupName '.' analysisName '.tonotopicMagnificaion.pCFcheck = overlayRoiData;']);
                        eval(['data.' roiSaveName '.' groupName '.' analysisName '.tonotopicMagnificaion.pathDistances = pathDistances;']);
                        
                        
                        % get overlay data from voxels within the gradient - tonotopic estimates averaged over cortical depths
                        interpMethod = 'linear';
                        for iOverlay = 1:3
                            clear tempOverlayData tempOverlayRoiData overlay2getName
                            switch analysisName
                                case 'glm_hrfDoubleGamma'
                                    overlay2getName = ['averageDepthVol_' Info.Sides{iSide} '(', glmoverlayNames{iOverlay} ')'];
                                    
                                case pRFanalysisName
                                    overlay2getName = ['averageDepthVol_' Info.Sides{iSide} '(', pRFoverlayNames{iOverlay}, ')'];
                            end
                            overlay2get = viewGet(thisView,'overlayNum',overlay2getName);
                            tempOverlayData  = viewGet(thisView,'overlaydata',1,overlay2get);
                            tempOverlayRoiData = interpn(tempOverlayData,verticesScanCoords{3}(1,:),verticesScanCoords{3}(2,:),verticesScanCoords{3}(3,:),interpMethod);
                            eval(['data.' roiSaveName '.' groupName '.' analysisName '.tonotopicMagnificaion.' estimateNames{iOverlay} '= tempOverlayRoiData;']);
                            
                        end
                        
                    end
                end
            end
        end
        
        % save data
        save(saveName,'data','-v7.3');
        
        % save view
        mrSaveView(thisView);
        
    end
    
    %% Quit current mrLoadRet view
    mrQuit()
    
end