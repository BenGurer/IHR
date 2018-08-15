function CM_groupAnalysisScript
%% CM_subjectAnalysisScript
% Scripted group analysis for Comparisons at 7T and Cortical Magnification (CM) studies
% Study dates: 2015 - 2018

% Info:
% gets data from subject analysis
% save structure
% data.roi.dataset.analysis.voxelParameter
% ie. data.Left.scan1.glm_boxcar_nCons8.R2
% tidys data in tidyverse
% saves as .csv

% Aims:
% Measure cortical magnification
% Compare acquisition protocols
% Compare analysis methods
% plot scanner noise

%% Tasks:
%% ROI analysis
%% Comparing analysis methods
% GLM vs pRF

%% Comparison 1: GLM methods
% Index max (IM), Centroid (C), debiased Weighted Mean (dWM)
% Data = Sparse concatenated

%% Comparison 2: pRF vs GLM
% GLM dWM , pRF pCF
% Data = Sparse concatenated

%% Comparison 3: pTW estimation
% GLM Beta weights
% Split half tuning width
% Data = Scans
% pRF
% Split half tuning widths
% Data = Scans
% Notes: Calculate subject level with matlab and plot group level with R

%% Data vis and Stats
% 	pCF distribution
%   pCF Reliability for methods
% 	Between run correlation:
% 		GLM IM, GLM C, GLM dWN, pRF pCF
% 	pCF Correlation between methods
%       prove they are valid - correlation matrix
% 	pCF maps
%       Group level or representative subject

%% Comparing acquisition protocols
% Sparse vs Continuous

%% Comparison 1: BOLD activity
% GLM Beta weights
% Split half tuning width
% 	Data = Scans
% ROI average beta weights
% 	Data = concatenated
% Ratio between average beta weights
% 	Data = concatenated

%% Comparison 2: pCF estimation
% pRF
% Data = concatenated
% pCF distribution
% pCF Correlation between
% Runs
% Which is better?
% Sparse and Continuous
% How similar are they

% GLM dWeightMean
% Data = concatenated
% pCF distribution
% pCF Correlation between
% Runs
% Which is better?
% Sparse and Continuous
% How similar are they?

%% Comparison 3: pTW estimation
% GLM Beta weights
% Split half tuning width
% Data = Scans
% pRF
% Split half tuning widths
% Data = Scans
% Notes: Calculate subject level with matlab and plot group level with R

%% Other Analysis
% HRF estimation
% 	Average deconvolution estimate
% 	Average fitted params ( check if median is more appropriate than mean (in case of outliers))
% 	Plot averaged fitted params or average of the curves
% Notes: Calculate subject level with matlab and plot group level with R


%% Print mrView figures

%% save data for R analysis/plotting

%% Begin

%% close and clear everything
clear all; close all; clc

disp('Running CM_groupAnalysisScript...')
disp('Group analysis for Cortical Magnificaiton and Comparisions at 7T studies...')
disp('This function gets subject data and saves it in tidyVerse format...')

%% get study info
[stimInfo, glmInfo, pRFInfo, Info, plotInfo] = CM_setupStudyParams;

% use ispc to set data directory
if ispc
    Info.dataDir = 'E:\OneDrive - The University of Nottingham\data';
else
    Info.dataDir = '/Volumes/data_PSY/OneDrive - The University of Nottingham/data';
end

% set variable to be ' to use with eval
q = char(39);

%% set logicals to what we want to do
doHRF = 1;
doComparisions = 0;
doCorticalMagnification = 0;
doStudyPlots = 1;

%% define subjects
% nSubjects = 8;
iSubs2Run = [1,2,3,4,5,6,7,8];

% iSubs2Run = [1,2,5,6,7,8];
%% HRF params
hrf_name = [];
hrf_params = [];
hrf_subject = [];

%% Varibles to define
% Voxel estimates
ve_subjectID = [];
ve_voxelID = [];
ve_roi = [];
ve_r2 = [];
ve_r2ATrue = [];
%  pCF estimate value
ve_frequency_nERB = []; % estimate values in number of ERB
ve_frequency_kHz = []; % estimate values in kHz
% pTW estimate value
ve_selectivity_nERB = []; % estimate values in number of ERB
ve_selectivity_kHz = []; % estimate values in kHz
ve_acquistion = []; % acquistion protocol
ve_concatenation = []; % concatenated or individual scan
ve_analysis = []; % GLM or pRF
ve_hrf = [];
ve_estimation = []; % for GLM {'Centriod'  'Spread'  'julien_pCF'  'julien_pTW'  'indexMax'} or pRF for pRF

% Split Tuning Curves
tc_subjectID = [];
tc_roi = [];
tc_beta_weight = [];
tc_beta_weight_norm  = [];
tc_beta_weight_ConATrue  = [];
tc_beta_weight_norm_ConATrue  = [];
tc_pRF_tuning_curve = [];
tc_pRF_tuning_curve_ConATrue = [];
tc_beta_bin_kHz = [];
tc_beta_freq_kHz = [];
tc_beta_bin_nERB = [];
tc_beta_freq_nERB = [];
tc_beta_bin_id = [];
tc_beta_freq_id = [];
tc_acquistion = []; % acquistion protocol
tc_analysis = []; % GLM or pRF
tc_hrf = [];
tc_beta_averaging = [];

% ROI Average Beta Weights
bw_subjectID = [];
bw_roi = [];
bw_beta_weight = [];
bw_beta_weight_norm = [];
bw_beta_freq_id = [];
bw_beta_freq_kHz = [ ];
bw_beta_freq_nERB = [ ];
bw_acquistion = []; % acquistion protocol
bw_analysis = []; % GLM or pRF
bw_hrf = [];
bw_beta_averaging = [];

% Coritcal Magnifiction
cm_CorticalDistance = [];
cm_Frequency_nERB = [];
cm_Analysis = [];
cm_ROI = [];
cm_r2 = [];
cm_TuningWidth = [];
cm_SubjectID = [];

% HRF
hrf_SubjectID = [];
hrf_Est = [];
hrf_EstNorm = [];
hrf_Time = [];

hrftw_SubjectID = [];
hrftw_Est = [];
hrftw_EstNorm = [];
hrftw_Time = [];
hrftw_freq = [];

hrfparams_name = [];
hrfparams_params = [];
hrfparams_subject = [];

hrfav_plot = [];

%% Varibles to define
% Voxel estimates
% roi
roiNames_data = {'LeftGR_GLM', 'LeftGRa_GLM', 'LeftGRp_GLM', 'RightGR_GLM', 'RightGRa_GLM', 'RightGRp_GLM'};
roiNames_table = {'Left', 'Left Anterior', 'Left Posterior', 'Right', 'Right Anterior', 'Right Posterior'};

% acquistion
acquistionNames = {'Sparse', 'Continuous'};

% concatenation
concatNames_data = {'ConcatenationSparse', 'scan_1', 'scan_3';...
    'ConcatenationCont',  'scan_2', 'scan_4'};
concatNames_table = {'Sparse', 'Sparse_1', 'Sparse_2';...
    'Continuous',  'Continuous_3', 'Continuous_4'};

% analysis
pRFanalysisName = ['pRF_', pRFInfo.pRFrestrictROI];
analysisNames_data = {'glm_hrfDoubleGamma', 'glm_hrfBoxcar', 'pRF'};
analysisNames_table= {'GLM', 'GLM', 'pRF'};
hrfNames_table = {'Double Gamma', 'Box Car', 'Difference of Gamma'};

% estimation method
estimateGLMFreqNames_data = {'Centriod', 'julien_pCF', 'indexMax'};
estimateGLMTuningName_data = {'Spread', 'julien_pTW', 'NA'};
estimateGLMFreqNames_table = {'Centriod', 'Debiased Centriod', 'Max'};
estimateGLMTuningName_table = {'Spread', 'Debiased Spread', 'NA'};

% Split Tuning Curves

% acquistion
acquistionNames_data = {'ConcatenationSparse','ConcatenationCont'};
acquistionNames_table = {'Sparse', 'Continuous'};

run_names = {'scan_1', 'scan_3';...
    'scan_2', 'scan_4'};

% analysis
beta_analysisNames_data = {'glm_hrfDoubleGamma'}; %  {'glm_hrfDoubleGamma', 'glm_hrfBoxcar'};
beta_analysisNames_table = {'GLM'}; % {'GLM', 'GLM'};
betaAveraging_data_tc = {'_nCons_8'}; % {'_nCons_32'}
betaAveragingNames_table_tc = {'Binned'}; % {'Moving Average'}
beta_hrfNames_table = {'Double Gamma'}; % {'Double Gamma', 'Box Car'};

% ROI Average Beta Weights
% used Split tuning curve variables plus:
betaAveragingNames_table_bw = {'None','Moving Average'}; %{'Binned'};


%% HRF
hrfSaveNames = {'x_Gamma', 'x_doubleGamma', 'x_dGamma',};
hrfNames = {'Gamma', 'doubleGamma', 'diffofGamma'};

%% get data and tidy
for iSub = 1:length(iSubs2Run)
    
    data = [];
    % Get subject info
    subjectInfo = get_SubjectInfo_CM(iSubs2Run(iSub));
    % Subject ID, flatmap names
    saveName = [subjectInfo.subjectID '_data.mat'];
    % move to subject folder
    cd(fullfile(Info.dataDir,Info.studyDir,subjectInfo.subjectID));
    % Load subject data
    load(saveName);
    
    subjectNumber = iSubs2Run(iSub);
    disp(['tidying subject ' num2str(subjectNumber) ' data...']);
    
    %% HRF
    if doHRF
        % average hrf params for 3T analysis
        % TidyVerse
        % Row: observation = subject
        % Column: Variable = fitted hrf param
                
        % HRF fitted estimates
        % save params to plot histograms
        % save HRF from fitted params to show shape
        % outside loop - calculate averaged params hrf and save that for plotting


        
        for iHRF = 1:length(hrfSaveNames)
            
%             temp_hrfparams_params = [];
            temp_hrfparams_params = eval(['data.hrf.',hrfSaveNames{iHRF} ';']);
            temp_hrfparams_name = repmat(hrfNames{iHRF},length(temp_hrfparams_params),1);
            temp_hrfparams_subject = repmat(iSub,length(temp_hrfparams_params),1);
            
            hrfparams_name = [hrf_name; string(temp_hrfparams_name)];
            hrfparams_params = [hrf_params; temp_hrfparams_params'];
            hrfparams_subject = [hrf_subject; temp_hrfparams_subject];
            
        end
        
        % HRF estimate
        temp_Est = data.hrf.deconv';
        temp_EstNorm = data.hrf.deconv' ./ max(data.hrf.deconv);
        temp_Time = data.hrf.estimate.time';
        nEstimates = length(temp_Est);
        temp_subjectID = repmat(subjectNumber,nEstimates,1);                                    
                                    
        hrf_SubjectID = [hrf_SubjectID; temp_subjectID];
        hrf_Est = [hrf_Est; temp_Est];
        hrf_EstNorm = [hrf_EstNorm; temp_EstNorm];
        hrf_Time = [hrf_Time; temp_Time];
        
        % HRF TW estimate
%         data.hrf.deconvTW
%         data.hrf.estimate.time
        frequencyBin = -7:7;
        temp_hrftw_Est = [];
        temp_hrftw_Est_norm = [];
        temp_hrftw_norm = data.hrf.deconvTW ./ max(max(data.hrf.deconvTW));
        temp_hrftw_freq = [];
        for iFreq = 1:length(frequencyBin)
            temp_hrftw_Est = [temp_hrftw_Est, data.hrf.deconvTW(iFreq,:)];
            temp_hrftw_Est_norm = [temp_hrftw_Est_norm, temp_hrftw_norm(iFreq,:)];            
            temp_hrftw_freq = [temp_hrftw_freq, repmat(frequencyBin(iFreq),1,length(data.hrf.deconvTW(iFreq,:)))];
        end
        
        temp_hrftw_Time = repmat(data.hrf.estimate.time,1,length(frequencyBin));
        nEstimates = length(temp_hrftw_Time);
        temp_hrftw_subjectID = repmat(subjectNumber,nEstimates,1); 
        
        hrftw_SubjectID = [hrftw_SubjectID; temp_hrftw_subjectID];
        hrftw_Est = [hrftw_Est; temp_hrftw_Est'];
        hrftw_EstNorm = [hrftw_EstNorm; temp_hrftw_Est_norm'];
        hrftw_Time = [hrftw_Time; temp_hrftw_Time'];
        hrftw_freq = [hrftw_freq; temp_hrftw_freq'];
        
        % Matlab
        % save hrf params in seperate variables, n by m, n = suject; m = parameters.
        
        for iHRF = 1:length(hrfSaveNames)
            temp_hrf_params = eval(['data.hrf.',hrfSaveNames{iHRF} ';']);
            eval(['hrf.',  hrfNames{iHRF},'_data(iSubs2Run(iSub),:) = temp_hrf_params;']);
        end
    end
    
    %% 7T comparisions %%
    if doComparisions
        
        %% save data in TidyVerse
        % dataset 1: pCF estiamtes
        % need one table with all pCF estiamtes
        
        %% Comparing analysis methods
        %   GLM vs pRF
        
        %   Data = Sparse concatenated and split
        
        %% Data vis and Stats
        % 	pCF distribution
        %   pCF Reliability for methods
        % 	Between run correlation:
        % 		GLM IM, GLM C, GLM dWN, pRF pCF
        % 	pCF Correlation between methods
        %       prove they are valid - correlation matrix
        % 	pCF maps
        %       Group level or representative subject
        %   Average split half tuning curves
        
        % Comparison 1: GLM methods
        %   Index max (IM), Centroid (C), debiased Weighted Mean (dWM)
        
        % Comparison 2: pRF vs GLM
        %   GLM dWM , pRF pCF
        
        % Comparison 3: pTW estimation
        %   GLM Beta weights  , pRF pTW
        %   Data = Scans / Split half tuning curves
        %   Notes: need to calcuate split half tuning curves
        
        
        %% Begin comparing: Voxel estimates
        % outcomes:
        % Distribution
        % Correlation
        
        % Transform data into tidyVerse
        % Row: observation = voxel
        % Column: Variable = everything else
        % subject, pCF estimate value, pTW estimate value, roi, acquistion, concatenation, analysis, hrf, estimation method
        
        disp(['Subject ' num2str(iSub) ': Tidying Voxel Estimates'])
        
        for iROI = 1:length(roiNames_data)
            roiName_data = roiNames_data{iROI};
            roiName_table = roiNames_table{iROI};
            for iAcq = 1:length(acquistionNames)
                acquistionName = acquistionNames{iAcq};
                
                for iConcat = 1:size(concatNames_data,2)
                    % get from group and scans
                    concatName_data = concatNames_data{iAcq,iConcat};
                    concatName_table = concatNames_table{iAcq,iConcat};
                    
                    for iAnal = 1:length(analysisNames_data)
                        analysisName_data = analysisNames_data{iAnal};
                        analysisName_table = analysisNames_table{iAnal};
                        
                        if ~strcmp(analysisName_data,'pRF')
                            estimateFreqNames_data = estimateGLMFreqNames_data;
                            estimateTuningName_data = estimateGLMTuningName_data;
                            estimateFreqNames_table = estimateGLMFreqNames_table;
                            estimateTuningName_table = estimateGLMTuningName_table;
                            if ~(iConcat == 1)
                                analysisName_data = [analysisNames_data{iAnal} ,'_nCons_32'];
                                analysisName_table = analysisNames_table{iAnal};
                            end
%                             constantOffset = 1; % correct for not converting to nERB correctly
                        else
                            estimateFreqNames_data = {'PrefCentreFreq'};
                            estimateTuningName_data = {'rfHalfWidth'};
                            estimateFreqNames_table = {'population Centre Frequency'};
                            estimateTuningName_table = {'population Tuning Width'};
                            
%                             constantOffset = 0;
                        end
                        
                        for iEst = 1:length(estimateFreqNames_data)
                            
                            % to save: subject, pCF estimate value, pTW estimate value, roi, acquistion, concatenation, analysis, hrf, estimation method
                            
                            % pCF estimate value
                            % get pCF estimate values from data struct
                            % repeat subject, roi, acquistion, concatenation, analysis, estimation method - nVoxel times
                            % save in table after subject for loop
                            eval(['temp_frequency = data.' roiName_data '.' concatName_data '.' analysisName_data '.' estimateFreqNames_data{iEst} ';']);
                            if iscell(temp_frequency)
                                temp_frequency_nERB = temp_frequency{:}; % estimate values in number of ERB
                            else
                                temp_frequency_nERB = temp_frequency;
                            end
%                             temp_frequency_nERB = temp_frequency_nERB - constantOffset;
                            temp_frequency_kHz = funInvNErb(temp_frequency_nERB); % estimate values in kHz
                            
                            % number of observations
                            nVoxels = length(temp_frequency_nERB);

                            if iAcq == 1 && iConcat == 1 && iAnal == 1 && iEst == 1
                                if iROI == 1
                                    totalVoxels = 1;
                                end
                                startVoxels = totalVoxels;
                                totalVoxels = totalVoxels + nVoxels - 1;
                            end
                            
                            % pTW estimate value
                            if ~strcmp(estimateTuningName_data{iEst},'NA')
                                eval(['temp_selectivity = data.' roiName_data '.' concatName_data '.' analysisName_data '.' estimateTuningName_data{iEst} ';']);
                                if iscell(temp_frequency)
                                    temp_selectivity_nERB = temp_selectivity{:}; % estimate values in number of ERB
                                else
                                    temp_selectivity_nERB = temp_selectivity;
                                end
%                                 temp_selectivity_nERB = temp_selectivity_nERB - constantOffset;
                                temp_selectivity_kHz = funInvNErb(temp_selectivity_nERB); % estimate values in kHz
                            else
                                temp_selectivity_nERB = nan(1,nVoxels);
                                temp_selectivity_kHz = nan(1,nVoxels);
                                
                            end                            
                            
                            eval(['temp_r2 = data.' roiName_data '.' concatName_data '.' analysisName_data '.r2;']);
                            if iscell(temp_r2)
                                temp_r2 = temp_r2{:}; % estimate values in number of ERB
                            else
                                temp_r2 = temp_r2;
                            end
                            if strcmp(analysisName_data,'glm_hrfDoubleGamma') && strcmp(concatName_data,'ConcatenationSparse')
                                temp_r2ATrue = temp_r2;
                            end
                            
                            % subject
                            temp_subjectID = repmat(subjectNumber,nVoxels,1);
                            % voxel                            
                            temp_voxelID = startVoxels:totalVoxels;
                            % roi
                            temp_roi = repmat(roiName_table,nVoxels,1);
                            % acquistion
                            temp_acquistion = repmat(acquistionName,nVoxels,1);
                            % concatenation
                            temp_concatenation = repmat(concatName_table,nVoxels,1);
                            % analysis
                            temp_analysis = repmat(analysisName_table,nVoxels,1);
                            % hrf
                            temp_hrf = repmat(hrfNames_table{iAnal},nVoxels,1);
                            % estimation method
                            temp_estimation = repmat(estimateFreqNames_table{iEst},nVoxels,1);
                            
                            % concatenate data
                            ve_frequency_nERB = [ve_frequency_nERB; temp_frequency_nERB'];
                            ve_frequency_kHz = [ve_frequency_kHz; temp_frequency_kHz'];
                            ve_selectivity_nERB = [ve_selectivity_nERB; temp_selectivity_nERB'];
                            ve_selectivity_kHz = [ve_selectivity_kHz; temp_selectivity_kHz'];
                            ve_r2 = [ve_r2; temp_r2'];
                            ve_r2ATrue = [ve_r2ATrue; temp_r2ATrue'];
                            ve_subjectID = [ve_subjectID; temp_subjectID];
                            ve_voxelID = [ve_voxelID; temp_voxelID'];
                            ve_roi = [ve_roi; string(temp_roi)];
                            ve_acquistion = [ve_acquistion; string(temp_acquistion)];
                            ve_concatenation = [ve_concatenation; string(temp_concatenation)];
                            ve_analysis = [ve_analysis; string(temp_analysis)];
                            ve_hrf = [ve_hrf; string(temp_hrf)];
                            ve_estimation = [ve_estimation; string(temp_estimation)];
                            
                        end
                    end
                end
            end
            %     end
            
            %% fin comparing: Voxel estimates
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %% Begin comparing: ROI average GLM beta weights
            % this will require calucating ROI averages before exporting to tidyVerse
            % split tuning Curves
            % moving average
            
            % ROI average estimate
            % outcomes:
            % Split tuning Curves
            % data: binned
            % ROI Average Beta Weights
            % data: glm moving average of all betas; hrf = Double Gamma; concatenated data
            
            % Split Tuning Curves
            % data: glm 8 bins; hrf = Double Gamma; split runs
            
            %% Split Tuning Curves
            % data: split
            %       glm 8 bins
            
            % Transform data into tidyVerse
            % Row = observation: beta weight at x frequency point
            % Column = Variable:
            % subject, beta weight value, frequency point, frequency group, roi, acquistion
            
            disp(['Subject ' num2str(iSub) ': Tidying Split Tuning Curves'])
            
            for iAcq = 1:length(acquistionNames_table)
                
                acquistionName = acquistionNames{iAcq};
                runA_name = run_names{iAcq,1};
                runB_name = run_names{iAcq,2};
                
                for iAnal = 1:length(beta_analysisNames_data)
                    
                    for iCons = 1:length(betaAveraging_data_tc)
                        
                        analysisName_data = [beta_analysisNames_data{iAnal}, betaAveraging_data_tc{iCons}];
                        analysisName_table = beta_analysisNames_table{iAnal};
                        
                        % get data
                        eval(['temp_betas_runA_cell = data.' roiName_data '.' runA_name '.' analysisName_data '.betas;']);
                        eval(['temp_betas_runB_cell = data.' roiName_data '.' runB_name '.' analysisName_data '.betas;']);
                        
                        eval(['temp_pTW_runA = data.' roiName_data '.' runA_name '.pRF.rfHalfWidth;']);
                        eval(['temp_pTW_runB = data.' roiName_data '.' runB_name '.pRF.rfHalfWidth;']);
                        
                        eval(['temp_pCF_runA = data.' roiName_data '.' runA_name '.pRF.PrefCentreFreq;']);
                        eval(['temp_pCF_runB = data.' roiName_data '.' runB_name '.pRF.PrefCentreFreq;']);                       
                        
                        % unpack data
                        % glm
                        temp_betas_runA = nan(size(temp_betas_runA_cell{1},2),length(temp_betas_runA_cell));
                        temp_betas_runB = nan(size(temp_betas_runA_cell{1},2),length(temp_betas_runA_cell));
                        for iBeta = 1:length(temp_betas_runA_cell)

                            temp_betas_runA(:,iBeta) = temp_betas_runA_cell{iBeta};
                            temp_betas_runB(:,iBeta) = temp_betas_runB_cell{iBeta};
                        end
                        
                        % pRF
                        
%                         temp_pTW_runA = temp_pTW_runA_cell{:};
%                         temp_pTW_runB = temp_pTW_runB_cell{:};
%                         
%                         temp_pCF_runA = temp_pCF_runA_cell{:};
%                         temp_pCF_runB = temp_pCF_runB_cell{:};
                        
                        [~, VoxelIndex_A] = max(temp_betas_runA,[],2);
                        [~, VoxelIndex_B] = max(temp_betas_runB,[],2);
                        
                        if iAcq ~= 1
                            VoxelIndex_A_ConATrue = round((VoxelIndex_A + VoxelIndex_B)./2);
                            VoxelIndex_B_ConATrue = round((VoxelIndex_A + VoxelIndex_B)./2);
                        else
                            VoxelIndex_A_ConATrue = VoxelIndex_A;
                            VoxelIndex_B_ConATrue = VoxelIndex_B;
                        end
                        
                        % compute split mean
                        splitA = nan(size(temp_betas_runA,2),size(temp_betas_runA,2));
                        splitB = nan(size(temp_betas_runA,2),size(temp_betas_runA,2));
                        splitA_ConATrue = splitA;
                        splitB_ConATrue = splitB;
                        
                        splitA_pRF_pTW = nan(size(temp_betas_runA,2),1);
                        splitB_pRF_pTW = splitA_pRF_pTW;
                        splitA_pCF_pTW = splitA_pRF_pTW;
                        splitB_pCF_pTW = splitA_pRF_pTW;
                        
                        for n = 1:size(temp_betas_runA,2)
                            splitA(n,:) = sum(temp_betas_runA(VoxelIndex_B==n,:))/sum(VoxelIndex_B==n);
                            splitB(n,:) = sum(temp_betas_runB(VoxelIndex_A==n,:))/sum(VoxelIndex_A==n);
                            
                            splitA_ConATrue(n,:) = sum(temp_betas_runA(VoxelIndex_B_ConATrue==n,:))/sum(VoxelIndex_B_ConATrue==n);
                            splitB_ConATrue(n,:) = sum(temp_betas_runB(VoxelIndex_A_ConATrue==n,:))/sum(VoxelIndex_A_ConATrue==n);
                            
                            % if we want pRF split mean - get data and average pTW here
                            % just get the pRF (pCf and pTW) params for the voxels at each bin and average
                            
                            splitA_pRF_pTW(n) = nanmean(temp_pTW_runA(VoxelIndex_B==n));
                            splitB_pRF_pTW(n) = nanmean(temp_pTW_runB(VoxelIndex_A==n));
                            
                            splitA_pRF_pCF(n) = nanmean(temp_pCF_runA(VoxelIndex_B==n));
                            splitB_pRF_pCF(n) = nanmean(temp_pCF_runB(VoxelIndex_A==n));
                            
                            splitA_pRF_pTW_ConATrue(n) = nanmean(temp_pTW_runA(VoxelIndex_B_ConATrue==n));
                            splitB_pRF_pTW_ConATrue(n) = nanmean(temp_pTW_runB(VoxelIndex_A_ConATrue==n));
                            
                            splitA_pRF_pCF_ConATrue(n) = nanmean(temp_pCF_runA(VoxelIndex_B_ConATrue==n));
                            splitB_pRF_pCF_ConATrue(n) = nanmean(temp_pCF_runB(VoxelIndex_A_ConATrue==n));
                            
                        end
                        splitMean = (splitA + splitB) / 2;
                        splitMeanNorm = splitMean./ max(max(splitMean));
                        
                        splitMean_pRF_pTW = (splitA_pRF_pTW + splitB_pRF_pTW) / 2;
                        splitMean_pRF_pCF = (splitA_pRF_pCF + splitB_pRF_pCF) / 2;
                        
                        splitMean_ConATrue = (splitA_ConATrue + splitB_ConATrue) / 2;
                        splitMeanNorm_ConATrue = splitMean_ConATrue./ max(max(splitMean_ConATrue));
                        
                        splitMean_pRF_pTW_ConATrue = (splitA_pRF_pTW_ConATrue + splitB_pRF_pTW_ConATrue) / 2;
                        splitMean_pRF_pCF_ConATrue = (splitA_pRF_pCF_ConATrue + splitB_pRF_pCF_ConATrue) / 2;                       
                        
                        % tidy data
                        c = 1;
                        pRF_tuning_curves = nan(size(temp_betas_runA,2),size(temp_betas_runA,2));
                        for iBeta = 1:size(splitMean,1)
                        % pRF                        
                        pRF_tuning_curves(iBeta,:) = def_Gaussian([splitMean_pRF_pCF(iBeta), splitMean_pRF_pTW(iBeta)],stimInfo.stimNERBs_bin);                        
                        pRF_tuning_curves_ConATrue(iBeta,:) = def_Gaussian([splitMean_pRF_pCF_ConATrue(iBeta), splitMean_pRF_pTW_ConATrue(iBeta)],stimInfo.stimNERBs_bin); 
                            for iFreq = 1:size(splitMean,2)
                                eval(['temp_beta_bin' num2str(iBeta) '_freq_' num2str(iFreq) ' = splitMean(iBeta,iFreq);'])
                                
                                temp_beta_weight(c) = splitMean(iBeta,iFreq);
                                temp_beta_weight_norm(c) = splitMeanNorm(iBeta,iFreq);
                                temp_pRF_tuning_curve(c) = pRF_tuning_curves(iBeta,iFreq);
                                temp_beta_weight_ConATrue(c) = splitMean_ConATrue(iBeta,iFreq);
                                temp_beta_weight_norm_ConATrue(c) = splitMeanNorm_ConATrue(iBeta,iFreq);                                
                                temp_pRF_tuning_curve_ConATrue(c) = pRF_tuning_curves_ConATrue(iBeta,iFreq);
                                temp_beta_bin_id(c) = iBeta;
                                temp_beta_freq_id(c) = iFreq;
                                temp_beta_bin_kHz(c) = round(stimInfo.stimFreqs_bin(iBeta),3);
                                temp_beta_freq_kHz(c) = round(stimInfo.stimFreqs_bin(iFreq),3);
                                temp_beta_bin_nERB(c) = round(stimInfo.stimNERBs_bin(iBeta),3);
                                temp_beta_freq_nERB(c) = round(stimInfo.stimNERBs_bin(iFreq),3);
                                c = c + 1;
                            end
                        end
                        
                        nVoxels = length(temp_beta_weight);
                        
                        % subject
                        temp_subjectID = repmat(subjectNumber,nVoxels,1);
                        % roi
                        temp_roi = repmat(roiName_table,nVoxels,1);
                        % acquistion
                        temp_acquistion = repmat(acquistionName,nVoxels,1);
                        % analysis
                        temp_analysis = repmat(analysisName_table,nVoxels,1);
                        % number of conditions
                        temp_beta_averaging = repmat(betaAveragingNames_table_tc{iCons},nVoxels,1);
                        % hrf
                        temp_hrf = repmat(beta_hrfNames_table{iAnal},nVoxels,1);
                        
                        % concatenate data
                        tc_beta_weight = [tc_beta_weight; temp_beta_weight'];
                        tc_beta_weight_norm = [tc_beta_weight_norm; temp_beta_weight_norm'];
                        tc_beta_weight_ConATrue = [tc_beta_weight_ConATrue; temp_beta_weight_ConATrue'];
                        tc_beta_weight_norm_ConATrue = [tc_beta_weight_norm_ConATrue; temp_beta_weight_norm_ConATrue'];                        
                        tc_pRF_tuning_curve = [tc_pRF_tuning_curve; temp_pRF_tuning_curve'];
                        tc_pRF_tuning_curve_ConATrue = [tc_pRF_tuning_curve_ConATrue; temp_pRF_tuning_curve_ConATrue'];
                        tc_beta_bin_id = [tc_beta_bin_id; temp_beta_bin_id'];
                        tc_beta_freq_id = [tc_beta_freq_id; temp_beta_freq_id'];
                        tc_beta_bin_kHz = [tc_beta_bin_kHz; temp_beta_bin_kHz'];
                        tc_beta_freq_kHz = [tc_beta_freq_kHz; temp_beta_freq_kHz'];
                        tc_beta_freq_nERB = [tc_beta_freq_nERB; temp_beta_freq_nERB'];
                        tc_beta_bin_nERB = [tc_beta_bin_nERB; temp_beta_bin_nERB'];                        
                        tc_subjectID = [tc_subjectID; temp_subjectID];
                        tc_roi = [tc_roi; string(temp_roi)];
                        tc_acquistion = [tc_acquistion; string(temp_acquistion)];
                        tc_analysis = [tc_analysis; string(temp_analysis)];
                        tc_beta_averaging = [tc_beta_averaging; string(temp_beta_averaging)];
                        tc_hrf = [tc_hrf; string(temp_hrf)];
                    end
                    
                    %% ROI Average Beta Weights
                    % data: concatenated
                    %       glm all cons moving average
                    
                    % Transform data into tidyVerse
                    % Row = observation: beta weight at x frequency point
                    % Column = Variable:
                    % subject, beta weight value, frequency point, roi, acquistion
                    
                    disp(['Subject ' num2str(iSub) ': Tidying Average Beta Weights'])
                    
                    acquistionName_data = acquistionNames_data{iAcq};
                    acquistionName_table = acquistionNames_table{iAcq};
                    
                    temp_beta_freq = [];
                    
                    for iAv = 1:length(betaAveragingNames_table_bw)
                        
                        analysisName_data = beta_analysisNames_data{iAnal};
                        analysisName_table = beta_analysisNames_table{iAnal};
                        
                        % get data
                        eval(['temp_betas_cell = data.' roiName_data '.' acquistionName_data '.' analysisName_data '.betas;']);
                        
                        temp_betas = nan(size(temp_betas_cell{1},2),length(temp_betas_cell));
                        
                        for iBeta = 1:length(temp_betas_cell)
                            % unpack data
                            temp_betas(:,iBeta) = temp_betas_cell{iBeta};
                        end
                        
                        if strcmp(betaAveragingNames_table_bw{iAv},betaAveragingNames_table_bw{1})
                            % mean
                            temp_beta_weight = mean(temp_betas);
                            
                            temp_beta_freq_kHz = round(stimInfo.stimFreqs,3);
                            temp_beta_freq_nERB = round(stimInfo.stimNERBs,3);
                        else
                            % moving average
                            nBins = 8;
                            windowAvSize = size(temp_betas,2)/nBins;
                            
                            loopLength = size(temp_betas,2) - windowAvSize + 1;
                            temp_betas_mv = nan(size(temp_betas,1),loopLength);
                            temp_beta_freq_kHz = nan(1,loopLength);
                            temp_beta_freq_nERB = nan(1,loopLength);
                            
                            for iBeta = 1:loopLength
                                temp_betas_mv(:,iBeta) = nanmean(temp_betas(:,iBeta:iBeta+windowAvSize-1),2);
                                temp_beta_freq_kHz(iBeta) = round(mean(stimInfo.stimFreqs(iBeta:iBeta+windowAvSize-1)),3);                                
                                temp_beta_freq_nERB(iBeta) = round(mean(stimInfo.stimNERBs(iBeta:iBeta+windowAvSize-1)),3);
                            end
                            
                            temp_beta_weight = mean(temp_betas_mv);                            

                        end
                        
                        temp_beta_weight_norm = temp_beta_weight ./ max(temp_beta_weight);
                        
                        nFrequencies = length(temp_beta_weight);
                        temp_beta_freq_id = 1:nFrequencies;
                        
                        % subject
                        temp_subjectID = repmat(subjectNumber,nFrequencies,1);
                        % roi
                        temp_roi = repmat(roiName_table,nFrequencies,1);
                        % acquistion
                        temp_acquistion = repmat(acquistionName_table,nFrequencies,1);
                        % analysis
                        temp_analysis = repmat(analysisName_table,nFrequencies,1);
                        % number of conditions
                        temp_beta_averaging = repmat(betaAveragingNames_table_bw{iAv},nFrequencies,1);
                        % hrf
                        temp_hrf = repmat(beta_hrfNames_table{iAnal},nFrequencies,1);                        
                        
                        % concatenate data
                        bw_beta_weight = [bw_beta_weight; temp_beta_weight'];
                        bw_beta_weight_norm = [bw_beta_weight_norm; temp_beta_weight_norm'];
                        bw_beta_freq_kHz = [bw_beta_freq_kHz; temp_beta_freq_kHz'];
                        bw_beta_freq_nERB = [bw_beta_freq_nERB; temp_beta_freq_nERB'];
                        bw_beta_freq_id = [bw_beta_freq_id; temp_beta_freq_id'];
                        bw_subjectID = [bw_subjectID; temp_subjectID];
                        bw_roi = [bw_roi; string(temp_roi)];
                        bw_acquistion = [bw_acquistion; string(temp_acquistion)];
                        bw_analysis = [bw_analysis; string(temp_analysis)];
                        bw_beta_averaging = [bw_beta_averaging; string(temp_beta_averaging)];
                        bw_hrf = [bw_hrf; string(temp_hrf)];
                    end
                    
                end
            end
            
        end
        
        %% fin comparing: ROI average GLM beta weights
    end
    
    %% Cortical Magnification %%
    if doCorticalMagnification
        pRFanalysisName = ['pRF_', pRFInfo.pRFrestrictROI];
        analysisNames = {'glm_hrfDoubleGamma',pRFanalysisName};
        analysisSaveName = {'GLM','pRF'};
        AP = {'a','p'};
        analName = {'GLM', 'pRF'};
        roiAnalName = {'GLM'};
        roiAnalNum = 1;

        for iSide = 1:length(Info.Sides)
            for iGroup = 1:length(glmInfo.groupNames)
                groupName = glmInfo.groupNames{iGroup};
                for iAnal = 1:length(analysisNames)
                    
                    analysisName_data = analysisNames{iAnal};
                    for iAP = 1:length(AP)
                        
                        roiSaveName = [Info.Sides{iSide}, 'GR' AP{iAP} '_' roiAnalName{roiAnalNum}];
                        roiName_data = [Info.Sides{iSide}, 'GR' AP{iAP}];
                        
                        eval(['tempCorticalDistance = data.' roiSaveName '.' groupName '.' analysisName_data '.tonotopicMagnificaion.relativeDistances(2,:);']);
                        eval(['tempFrequency = data.' roiSaveName '.' groupName '.' analysisName_data '.tonotopicMagnificaion.pCF;']);
                        eval(['tempFrequencycheck = data.' roiSaveName '.' groupName '.' analysisName_data '.tonotopicMagnificaion.pCFcheck{3};']);
                        eval(['tempTuningWidth = data.' roiSaveName '.' groupName '.' analysisName_data '.tonotopicMagnificaion.pTW;']);
                        eval(['tempR2 = data.' roiSaveName '.' groupName '.' analysisName_data '.tonotopicMagnificaion.r2;']);
                        
                        nVoxels = length(tempFrequency);
                        tempAnalysis = repmat(analysisSaveName{iAnal},nVoxels,1);
                        tempROI = repmat(roiName_data,nVoxels,1);                        
                              
                        % subject
                        temp_subjectID = repmat(subjectNumber,nVoxels,1);
                        
                        cm_r2 = [cm_r2; tempR2'];
                        cm_TuningWidth = [cm_TuningWidth; tempTuningWidth'];
                        cm_CorticalDistance = [cm_CorticalDistance; tempCorticalDistance'];
                        cm_Frequency_nERB = [cm_Frequency_nERB; tempFrequency'];
                        if isempty(cm_Analysis)
                            cm_Analysis = tempAnalysis;
                            cm_ROI = tempROI;
                        else
                            cm_Analysis = char(cm_Analysis,tempAnalysis);
                            cm_ROI = char(cm_ROI,tempROI);
                        end
                        
                        cm_SubjectID = [cm_SubjectID; temp_subjectID];
                        
                    end
                    
                end
            end
        end
        
        
    end
    
end % subject loop end

%% save data to file
disp('saving data...')
% move to group data save location
cd(fullfile(Info.dataDir,Info.studyDir,'groupAnalysis'));

%% Save Comparisions
if doComparisions
    % Voxel estimates
    voxel_estimates_table = [];
    voxel_estimates_table = table(ve_subjectID,ve_voxelID,...
        ve_frequency_nERB, ve_frequency_kHz,...
        ve_selectivity_nERB, ve_selectivity_kHz,...
        ve_roi, ve_r2, ve_r2ATrue,...
        ve_acquistion, ve_concatenation,...
        ve_analysis, ve_hrf, ve_estimation,...
        'VariableNames',{'subjectID','voxelID',...
        'frequency_nERB', 'frequency_kHz', ...
        'selectivity_nERB', 'selectivity_kHz',...
        'roi','r2','r2ATrue',...
        'acquistion', 'concatenation',...
        'analysis', 'hrf',...
        'estimation'});
    
    writetable(voxel_estimates_table, 'Comparisions_voxel_estimates.csv')
    % Tuning curves
    tuning_curves_table = [];
    tuning_curves_table = table(tc_subjectID,...
        tc_beta_weight, tc_beta_weight_norm,...
        tc_beta_weight_ConATrue, tc_beta_weight_norm_ConATrue,...
        tc_pRF_tuning_curve, tc_pRF_tuning_curve_ConATrue,...
        tc_beta_bin_nERB,...
        tc_beta_freq_nERB,...
        tc_beta_bin_kHz,...
        tc_beta_freq_kHz,...
        tc_beta_bin_id,...
        tc_beta_freq_id,...
        tc_roi,...
        tc_acquistion, tc_beta_averaging,...
        tc_analysis, tc_hrf,...
        'VariableNames',{'subjectID',...
        'beta_weight', 'beta_weight_normalised',...
        'beta_weight_A_True', 'beta_weight_A_True_normalised',...
        'pRF_tuning_curve', 'pRF_tuning_curve_ConATrue',...
        'beta_bin_NERB', ...
        'beta_freq_NERB',...
        'beta_bin_kHz', ...
        'beta_freq_kHz',...
        'beta_bin_ID', ...
        'beta_freq_ID',...
        'roi',...
        'acquistion', 'beta_averaging',...
        'analysis', 'hrf'});
    
    
    writetable(tuning_curves_table, 'Comparisions_tuning_curves.csv')
    
    % average beta weights
    av_beta_weights_table = [];
    av_beta_weights_table = table(bw_subjectID,...
        bw_beta_weight, bw_beta_weight_norm,... 
        bw_beta_freq_nERB, bw_beta_freq_kHz, bw_beta_freq_id, bw_roi,...
        bw_acquistion, bw_beta_averaging,...
        bw_analysis, bw_hrf,...
        'VariableNames',{'subjectID',...
        'beta_weight', 'beta_weight_normalised',...
        'beta_freq_NERB', 'beta_freq_kHz', 'beta_freq_ID',...
        'roi',...
        'acquistion', 'beta_averaging',...
        'analysis', 'hrf'});
    
    writetable(av_beta_weights_table, 'Comparisions_beta_weights.csv')
end

%% Save HRF
if doHRF
    
    % subject loop would end here
    % move out of subject loop to save all subjects
    HRFparams_table = [];
    HRFparams_table = table(hrfparams_name,hrfparams_params,...
        hrfparams_subject,...
        'VariableNames',{'HRFparams_name', 'HRFparams_params', 'HRFparams_subject'});
        
    writetable(HRFparams_table, 'HRFparams_est.csv')
    
    hrf_table = [];
    hrf_table = table(hrf_SubjectID,hrf_Est,...
        hrf_EstNorm,hrf_Time,...
        'VariableNames',{'Subject', 'hrf_Est', 'hrf_EstNorm', 'hrf_Time'});    
    
    writetable(hrf_table, 'HRF_est.csv')
    
    hrftw_table = [];
    hrftw_table = table(hrftw_SubjectID, hrftw_Est, hrftw_EstNorm,...
        hrftw_freq, hrftw_Time,...
        'VariableNames',{'Subject', 'hrf_Est', 'hrf_EstNorm', 'hrf_freq', 'hrf_Time'});    
    
    writetable(hrftw_table, 'HRFTW_est.csv')
            
    % check param distribution
    % take median - minise the effects of outliers

    %% hrf params plotted
    hrfFunctions = {@get_HRFGamma;
        @get_HRFDoubleGamma;
        @get_HRFDiffOfGamma;};

    t = 0:20;
    hrfav_plot = [];
    for iHRF = 1:length(hrfNames)
        
        % get data
        hrf2av = [];
        eval(['hrf2av = hrf.',  hrfNames{iHRF},'_data;']);
        
        % check by plotting functions
        % average using median
        figure
        for iParam = 1:size(hrf2av,2)
            subplot(3,3,iParam)
            hist(hrf2av(:,iParam))
            hrf_av(iParam) = median(hrf2av(:,iParam));
        end
        
        eval(['hrf.',  hrfNames{iHRF},'_av = hrf_av;']);
        
        % check average functions
        subplot(3,3,9)
        plot(t,hrfFunctions{iHRF}(hrf_av,t))
        disp(hrf_av)
        
        temp_hrf_av_plot = [];
        temp_hrf_av_plot = hrfFunctions{iHRF}(hrf_av,t);
        
        hrfparamsav_plot = [hrfparamsav_plot, temp_hrf_av_plot'];
        hrfparamsav_time = [hrfparamsav_time, t']
        hrfparamsav_name = [hrfparamsav_time,  string(repmat(hrfNames{iHRF},length(temp_hrf_av_plot),1))]              
        
        % make param name list and save each param with its name
    end
    
     hrfparamsav_table = [];
    hrfparamsav_table = table(hrfparamsav_plot, hrfparamsav_time, hrfparamsav_name,...
        'VariableNames',{'Data', 'Time', 'Function name'});    
    
    % save hrf
    save('hrf.mat','hrf')
      %%%% to here £££££$$$£"T
    
end

%% Save Cortical Magnification 
if doCorticalMagnification
    
    T = table(cm_CorticalDistance,cm_Frequency_nERB,...
        cm_r2, cm_TuningWidth,...
        cm_Analysis,cm_ROI,...
        cm_SubjectID,...
        'VariableNames',{'CorticalDistance' 'Frequency' 'r2' 'TuningWidth' 'Analysis' 'ROI','SubjectID'});
    
    writetable(T, 'CorticalMagnification.csv')
    
end

%% Study stuff to plot
% could just add to the end of data as a different acquistion?
if doStudyPlots
    noise_temp_level = CM_plotScannerNoise;    
    noise_temp_freq = stimInfo.stimFreqs;
    
    % moving average
    nBins = 8;
    windowAvSize = size(noise_temp_freq,2)/nBins;
    loopLength = size(noise_temp_freq,2) - windowAvSize + 1;
    noise_temp_level_mv = nan(size(noise_temp_freq,1),loopLength);
    noise_temp_freq_mv = nan(1,loopLength);
    
    for iBeta = 1:loopLength
        noise_temp_level_mv(iBeta) = round(mean(noise_temp_level(iBeta:iBeta+windowAvSize-1)),3);
        noise_temp_freq_mv(iBeta) = round(mean(noise_temp_freq(iBeta:iBeta+windowAvSize-1)),3);
    end
    
    noise_freq = [noise_temp_level'; noise_temp_level_mv'];
    noise_level = [noise_temp_freq'; noise_temp_freq_mv'];
    
    hrf_table = [];
    hrf_table = table(noise_freq,noise_level,...
        'VariableNames',{'noise_freq', 'noise_level'});
    
    writetable(hrf_table, 'CM_scanner_noise.csv')    
end

%% Fin.
disp('Done!')

end
