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
doHRF = 0;
doComparisions = 1;
doCorticalMagnification = 0;

%% define subjects
% nSubjects = 8;
% iSubs2Run = [1,2,3,4,5,6,7,8];

iSubs2Run = [1,2,5,6,7,8];
%% HRF params
hrf_name = [];
hrf_params = [];
hrf_subject = [];

%% Varibles to define
% Voxel estimates
ve_subjectID = [];
ve_roi = [];
ve_r2 = [];
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
tc_beta_bin = [];
tc_beta_freq = [];
tc_acquistion = []; % acquistion protocol
tc_analysis = []; % GLM or pRF
tc_hrf = [];
tc_beta_averaging = [];

% ROI Average Beta Weights
bw_subjectID = [];
bw_roi = [];
bw_beta_weight = [];
bw_beta_weight_norm = [];
bw_beta_bin = [];
bw_beta_freq = [];
bw_acquistion = []; % acquistion protocol
bw_analysis = []; % GLM or pRF
bw_hrf = [];
bw_beta_averaging = [];

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

betaAveragingNames_table_bw = {'None','Moving Average'}; %{'Binned'};

% acquistion
%             acquistionNames_data = {'ConcatenationSparse','ConcatenationCont'};
%             acquistionNames_table = {'Sparse', 'Continuous'};

% analysis
% beta_analysisNames_data = {'glm_hrfDoubleGamma'}; %  {'glm_hrfDoubleGamma', 'glm_hrfBoxcar'};
% beta_analysisNames_table = {'GLM'}; % {'GLM', 'GLM'};
%     betaAveraging_data = {'_nCons_8'}; % {'_nCons_32'}

% could bin at this stage?
%     betaAveragingNames_table = {'None'}; %{'Binned'};
% hrfNames_table = {'Double Gamma'}; % {'Double Gamma', 'Box Car'};


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
        
        hrfSaveNames = {'x_Gamma', 'x_doubleGamma', 'x_dGamma',};
        hrfNames = {'Gamma', 'doubleGamma', 'diffofGamma'};
        hrfFunctions = {@get_HRFGamma;
            @get_HRFDoubleGamma;
            @get_HRFDiffOfGamma;};
        t = 0:15;
        
        for iHRF = 1:length(hrfSaveNames)
            
            temp_hrf_params = [];
            temp_hrf_params = eval(['data.hrf.',hrfSaveNames{iHRF} ';']);
            temp_hrf_name = repmat(hrfNames{iHRF},length(temp_hrf_params),1);
            temp_hrf_subject = repmat(iSub,length(temp_hrf_params),1);
            
            hrf_name = [hrf_name; string(temp_hrf_name)];
            hrf_params = [hrf_params; temp_hrf_params'];
            hrf_subject = [hrf_subject; temp_hrf_subject];
            
        end
        
        % subject loop would end here
        % move out of subject loop to save all subjects
        hrf_table = [];
        hrf_table = table(hrf_name,hrf_params,...
            hrf_subject,...
            'VariableNames',{'hrf_name', 'hrf_params', 'hrf_subject'});
        
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
                        else
                            estimateFreqNames_data = {'PrefCentreFreq'};
                            estimateTuningName_data = {'rfHalfWidth'};
                            estimateFreqNames_table = {'population Centre Frequency'};
                            estimateTuningName_table = {'population Tuning Width'};
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
                            temp_frequency_kHz = funInvNErb(temp_frequency_nERB); % estimate values in kHz
                            
                            % number of observations
                            nVoxels = length(temp_frequency_nERB);
                            
                            % pTW estimate value
                            if ~strcmp(estimateTuningName_data{iEst},'NA')
                                eval(['temp_selectivity = data.' roiName_data '.' concatName_data '.' analysisName_data '.' estimateTuningName_data{iEst} ';']);
                                if iscell(temp_frequency)
                                    temp_selectivity_nERB = temp_selectivity{:}; % estimate values in number of ERB
                                else
                                    temp_selectivity_nERB = temp_selectivity;
                                end
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
                            
                            % subject
                            temp_subjectID = repmat(subjectNumber,nVoxels,1);
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
                            ve_subjectID = [ve_subjectID; temp_subjectID];
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
            
            
            for iAcq = 1:length(acquistionNames_table)
                
                acquistionName = acquistionNames{iAcq};
                runA_name = run_names{iAcq,1};
                runB_name = run_names{iAcq,2};
                
                for iAnal = 1:length(beta_analysisNames_data)
                    
                    for iCons = 1:length(betaAveraging_data_tc)
                        
                        analysisName_data = [beta_analysisNames_data{iAnal}, betaAveraging_data_tc{iCons}];
                        analysisName_table = beta_analysisNames_table{iAnal};
                        
                        eval(['temp_betas_runA_cell = data.' roiName_data '.' runA_name '.' analysisName_data '.betas;']);
                        eval(['temp_betas_runB_cell = data.' roiName_data '.' runB_name '.' analysisName_data '.betas;']);
                        
                        temp_betas_runA = nan(size(temp_betas_runA_cell{1},2),length(temp_betas_runA_cell));
                        temp_betas_runB = nan(size(temp_betas_runA_cell{1},2),length(temp_betas_runA_cell));
                        for iBeta = 1:length(temp_betas_runA_cell)
                            % unpack data
                            temp_betas_runA(:,iBeta) = temp_betas_runA_cell{iBeta};
                            temp_betas_runB(:,iBeta) = temp_betas_runB_cell{iBeta};
                        end
                        
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
                        
                        for n = 1:size(temp_betas_runA,2)
                            splitA(n,:) = sum(temp_betas_runA(VoxelIndex_B==n,:))/sum(VoxelIndex_B==n);
                            splitB(n,:) = sum(temp_betas_runB(VoxelIndex_A==n,:))/sum(VoxelIndex_A==n);
                            
                            splitA_ConATrue(n,:) = sum(temp_betas_runA(VoxelIndex_B_ConATrue==n,:))/sum(VoxelIndex_B_ConATrue==n);
                            splitB_ConATrue(n,:) = sum(temp_betas_runB(VoxelIndex_A_ConATrue==n,:))/sum(VoxelIndex_A_ConATrue==n);
                            
                            % if we want pRF split mean - get data and average pTW here
                            % just get the pRF (pCf and pTW) params for the voxels at each bin and average
                        end
                        splitMean = (splitA + splitB) / 2;
                        splitMeanNorm = splitMean./ max(max(splitMean));
                        
                        splitMean_ConATrue = (splitA_ConATrue + splitB_ConATrue) / 2;
                        splitMeanNorm_ConATrue = splitMean_ConATrue./ max(max(splitMean_ConATrue));
                        
                        % tidy data
                        c = 1;
                        for iBeta = 1:size(splitMean,1)
                            for iFreq = 1:size(splitMean,2)
                                eval(['temp_beta_bin' num2str(iBeta) '_freq_' num2str(iFreq) ' = splitMean(iBeta,iFreq);'])
                                
                                temp_beta_weight(c) = splitMean(iBeta,iFreq);
                                temp_beta_weight_norm(c) = splitMeanNorm(iBeta,iFreq);
                                temp_beta_weight_ConATrue(c) = splitMean_ConATrue(iBeta,iFreq);
                                temp_beta_weight_norm_ConATrue(c) = splitMeanNorm_ConATrue(iBeta,iFreq);
                                temp_beta_bin(c) = iBeta;
                                temp_beta_freq(c) = iFreq;
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
                        tc_beta_bin = [tc_beta_bin; temp_beta_bin'];
                        tc_beta_freq = [tc_beta_freq; temp_beta_freq'];
                        tc_subjectID = [tc_subjectID; temp_subjectID];
                        tc_roi = [tc_roi; string(temp_roi)];
                        tc_acquistion = [tc_acquistion; string(temp_acquistion)];
                        tc_analysis = [tc_analysis; string(temp_analysis)];
                        tc_beta_averaging = [tc_beta_averaging; string(temp_beta_averaging)];
                        tc_hrf = [tc_hrf; string(temp_hrf)];
                    end
                    %                 end
                    
                    %             end
                    %     end
                    
                    %% ROI Average Beta Weights
                    % data: concatenated
                    %       glm all cons moving average
                    
                    % Transform data into tidyVerse
                    % Row = observation: beta weight at x frequency point
                    % Column = Variable:
                    % subject, beta weight value, frequency point, roi, acquistion
                    
                    
                    %     for iROI = 1:length(roiNames_data)
                    %         roiName_data = roiNames_data{iROI};
                    %         roiName_table = roiNames_table{iROI};
                    %             for iAcq = 1:length(acquistionNames_table)
                    %                 acquistionName_data = acquistionNames_data{iAcq};
                    acquistionName_table = acquistionNames_table{iAcq};
                    
                    %                 for iAnal = 1:length(beta_analysisNames_data)
                    
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
                        else
                            % moving average
                            nBins = 8;
                            windowAvSize = size(temp_betas,2)/nBins;
                            
                            loopLength = size(temp_betas,2) - windowAvSize + 1;
                            temp_betas_mv = nan(size(temp_betas,1),loopLength);
                            
                            for iBeta = 1:loopLength
                                temp_betas_mv(:,iBeta) = nanmean(temp_betas(:,iBeta:iBeta+windowAvSize-1),2);
                            end
                            temp_beta_weight = mean(temp_betas_mv);
                        end
                        
                        temp_beta_weight_norm = temp_beta_weight ./ max(temp_beta_weight);
                        
                        %
                        nFrequencies = length(temp_beta_weight);
                        
                        temp_beta_freq = 1:nFrequencies;
                        
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
                        bw_beta_freq = [bw_beta_freq; temp_beta_freq'];
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
        CorticalDistance = [];
        Frequency = [];
        Analysis = [];
        ROI = [];
        ve_r2 = [];
        TuningWidth = [];
        for iSide = 1:length(Info.Sides)
            for iGroup = 1:length(glmInfo.groupNames)
                groupName = glmInfo.groupNames{iGroup};
                for iAnal = 1:length(analysisNames)
                    
                    analysisName_data = analysisNames{iAnal};
                    for iAP = 1:length(AP)
                        
                        roiSaveName = [Info.Sides{iSide}, 'GR' AP{iAP} '_' analName{iAnal}];
                        roiName_data = [Info.Sides{iSide}, 'GR' AP{iAP}];
                        
                        eval(['tempCorticalDistance = data.' roiSaveName '.' groupName '.' analysisName_data '.tonotopicMagnificaion.relativeDistances(2,:);']);
                        eval(['tempFrequency = data.' roiSaveName '.' groupName '.' analysisName_data '.tonotopicMagnificaion.pCF;']);
                        eval(['tempFrequencycheck = data.' roiSaveName '.' groupName '.' analysisName_data '.tonotopicMagnificaion.pCFcheck{3};']);
                        eval(['tempTuningWidth = data.' roiSaveName '.' groupName '.' analysisName_data '.tonotopicMagnificaion.pTW;']);
                        eval(['tempR2 = data.' roiSaveName '.' groupName '.' analysisName_data '.tonotopicMagnificaion.r2;']);
                        
                        nVoxels = length(tempFrequency);
                        tempAnalysis = repmat(analysisSaveName{iAnal},nVoxels,1);
                        tempROI = repmat(roiName_data,nVoxels,1);
                        
                        ve_r2 = [ve_r2; tempR2'];
                        TuningWidth = [TuningWidth; tempTuningWidth'];
                        CorticalDistance = [CorticalDistance; tempCorticalDistance'];
                        Frequency = [Frequency; tempFrequency'];
                        if isempty(Analysis)
                            Analysis = tempAnalysis;
                            ROI = tempROI;
                        else
                            Analysis = char(Analysis,tempAnalysis);
                            ROI = char(ROI,tempROI);
                        end
                        
                        
                    end
                    
                end
            end
        end
        
        subject = repmat(iSub,length(CorticalDistance),1);
        
    end
    
end % subject loop end

%% save data to file

% move to group data save location
cd(fullfile(Info.dataDir,Info.studyDir,'groupAnalysis'));

if doComparisions
    % Voxel estimates
    voxel_estimates_table = [];
    voxel_estimates_table = table(ve_subjectID,...
        ve_frequency_nERB, ve_frequency_kHz,...
        ve_selectivity_nERB, ve_selectivity_kHz,...
        ve_roi, ve_r2,...
        ve_acquistion, ve_concatenation,...
        ve_analysis, ve_hrf, ve_estimation,...
        'VariableNames',{'subjectID',...
        'frequency_nERB', 'frequency_kHz', ...
        'selectivity_nERB', 'selectivity_kHz',...
        'roi','r2',...
        'acquistion', 'concatenation',...
        'analysis', 'hrf',...
        'estimation'});
    
    writetable(voxel_estimates_table, 'Comparisions_voxel_estimates.csv')
    % Tuning curves
    tuning_curves_table = [];
    tuning_curves_table = table(tc_subjectID,...
        tc_beta_weight, tc_beta_weight_norm,...
        tc_beta_weight_ConATrue,tc_beta_weight_norm_ConATrue,...
        tc_beta_bin,...
        tc_beta_freq,...
        tc_roi,...
        tc_acquistion, tc_beta_averaging,...
        tc_analysis, tc_hrf,...
        'VariableNames',{'subjectID',...
        'beta_weight', 'beta_weight_normalised',...
        'beta_weight_A_True', 'beta_weight_A_True_normalised',...
        'beta_bin', ...
        'beta_freq',...
        'roi',...
        'acquistion', 'beta_averaging',...
        'analysis', 'hrf'});
    
    
    writetable(tuning_curves_table, 'Comparisions_tuning_curves.csv')
    
    % average beta weights
    av_beta_weights_table = [];
    av_beta_weights_table = table(bw_subjectID,...
        bw_beta_weight, bw_beta_weight_norm,... 
        bw_beta_freq,bw_roi,...
        bw_acquistion, bw_beta_averaging,...
        bw_analysis, bw_hrf,...
        'VariableNames',{'subjectID',...
        'beta_weight', 'beta_weight_normalised',...
        'beta_freq',...
        'roi',...
        'acquistion', 'beta_averaging',...
        'analysis', 'hrf'});
    
    writetable(av_beta_weights_table, 'Comparisions_beta_weights.csv')
end

% HRF
if doHRF
    % check param distribution
    % take median - minise the effects of outliers
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
    end
    
    % save hrf
    save('hrf.mat','hrf')
end

%% Cortical Magnification %%
if doCorticalMagnification
    
    T = table(CorticalDistance,Frequency,...
        ve_r2, TuningWidth,...
        Analysis,ROI,...
        subject,...
        'VariableNames',{'CorticalDistance' 'Frequency' 'r2' 'TuningWidth' 'Analysis' 'ROI','Subject'});
    
    writetable(T, [saveName, '_CM.csv'])
    
end

%% Fin.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%     TimePoints = [];
%     hrfEstimate = [];
%     hrfEstimateName = [];
%
%     % get_HRFDoubleGamma(x_doubleGamma,t)
%     % get_HRFGamma(x_Gamma,t)
%     % get_HRFDiffOfGamma(x_dGamma,t)
%
%
%     % Two data frames
%     % one to average curves
%     % one to average fits?
%     % OR
%     % fit to average curves
%     % average and fit in matlab??
%     % or export fitted curves as functions and then fit to them?
%
%     hrfFunctions = {@get_HRFGamma;
%         @get_HRFDoubleGamma;
%         @get_HRFDiffOfGamma;
%         @get_HRFBoxCar};
%
%     for iHRF = 1:length(hrfNames)
%
%         tempTimePoints = data.hrf.estimate.time;
%
%         if ~strcmpi(hrfSaveNames{iHRF},hrfSaveNames{end})
%             eval(['temphrfFit = data.hrf.' hrfSaveNames{iHRF} ';']);
%             temphrfEstimate = hrfFunctions{iHRF}(temphrfFit,tempTimePoints);
%         elseif hrfSaveNames{iHRF} == 'BoxCar'
%             temphrfEstimate = hrfFunctions{iHRF}([2.5, 2.5],tempTimePoints);
%
%         else
%             eval(['temphrfEstimate = data.hrf.' hrfSaveNames{iHRF} ';']);
%         end
%
%         temphrfEstimateName = repmat(hrfNames{iHRF},nTimepoints,1);
%
%         TimePoints = [TimePoints, tempTimePoints'];
%         hrfEstimate = [hrfEstimate, temphrfEstimate'];
%
%         if isempty(hrfEstimateName)
%             hrfEstimateName = temphrfEstimateName;
%         else
%             hrfEstimateName = char(hrfEstimateName,temphrfEstimateName);
%         end
%
%     end
%     subject = repmat(iSub,length(TimePoints),1);
%
%     T_hrf = table(hrfEstimate,TimePoints,...
%         hrfEstimateName,...
%         subject,...
%         'VariableNames',{'hrfEstimate', 'Gamma', 'doubleGamma', 'differenceOfGamma', 'Time(s)', 'Subject'});
%
%     writetable(T, 'HRF_data.csv')
%
%     % just HRF estimate
%     nTimepoints = length(data.hrf.estimate.time);
%
%     TimePoints = [];
%     hrfEstimate = [];
%     hrfEstimateName = [];
%     tempTimePoints = data.hrf.estimate.time;
%     temphrfEstimate = data.hrf.deconv;
%
%     TimePoints = [TimePoints, tempTimePoints'];
%     hrfEstimate = [hrfEstimate, temphrfEstimate'];
%
%
%     subject = repmat(iSub,length(TimePoints),1);
%
%     T = table(hrfEstimate,TimePoints,...
%         subject,...
%         'VariableNames',{'hrfBetaEstimate', 'Time_sec', 'Subject'});
%
%     writetable(T, [saveName, '_HRF.csv'])
%
%    % average estimate for science
%        %     data.hrf.x_doubleGamma, data.hrf.x_Gamma, data.hrf.x_dGamma, data.hrf.estimate, data.hrf.deconv, data.hrf.deconvTW
%        hrfNames = {'Gamma', 'Double Gamma', 'Diff of Gamma', 'BoxCar', 'Deconvolution'};
%     hrfSaveNames = {'x_Gamma', 'x_doubleGamma', 'x_dGamma', 'deconv'};
%     nTimepoints = length(data.hrf.estimate.time);

%% HRF params
% need to:
% get data for subject
% tidy
% save to group data table in tidy format - this is then loaded by R or do group analysis

% average params using median - plot histogram to check

%% save average recentred Deconvoluction



%% 7T comparisions %%
%
% % compare:
% %   pCF estimates
%
% %% Concatenated data
% Frequency = [];
% Frequency_kHz = [];
% TuningWidth = [];
% TuningWidth_kHz = [];
% Analysis = [];
% ROI = [];
% ve_r2 = [];
% % voxelPropertyNames: {'Centriod'  'Spread'  'julien_pCF'  'julien_pTW'  'indexMax'}
% % pRFOverlayNames: {'r2'  'PrefCentreFreq'  'rfHalfWidth'}
% estimateFreqNames = {'Centriod', 'julien_pCF', 'indexMax'};
% estimateTuningName = {'Spread', 'julien_pTW', 'NA'};
%
% analysisType = {'GLM', 'pRF'};
% hrfType = {'boxcar', 'Gamma'};
%
%
% pRFanalysisName = ['pRF_', pRFInfo.pRFrestrictROI];
% analysisNames = {'glm_hrfDoubleGamma', 'glm_hrfBoxCar', pRFanalysisName};
% nCons = [32, 8];
% nScans = 4;
%
% for iSide = 1:length(Info.Sides)
%     for iGroup = 1:length(glmInfo.groupNames)
%         groupName = glmInfo.groupNames{iGroup};
%         for iAnal = 1:length(analysisType)
%
%             analysisName_data = analysisNames{iAnal};
%             if analysisName_data ~= analysisNames{3}
%
%                 estimateFreqNames = {'Centriod', 'julien_pCF', 'indexMax'};
%                 estimateTuningName = {'Spread', 'julien_pTW', 'NA'};
%
%             else
%
%                 %%%%%%% FIND pRF save names %%%%%%%%%%%%%%%%
%                 estimateFreqNames = {'pCF'};
%                 estimateTuningName = {'pTW'};
%
%             end
%
%             for iHRF = 1:length(hrfType)
%                 %                 roiSaveName = [Info.Sides{iSide}, 'GR_' analName{iAnal}];
%                 roiSaveName = [Info.Sides{iSide}, 'GR_GLM']; % compare using ROI data from derived GLM analysis
%                 roiName_data = [Info.Sides{iSide}, 'GR'];
%
%                 eval(['tempFrequency = data.' roiSaveName '.' groupName '.' analysisName_data '.' estimateFreqNames{iEst}]);
%                 tempFrequency_kHz = funInvNErb(tempFrequency);
%
%                 if ~strcmp('NA',estimateTuningName{iEst})
%                     eval(['tempTuningWidth = data.' roiSaveName '.' groupName '.' analysisName_data '.' estimateTuningName{iEst}]);
%                 else
%                     tempTuningWidth = nan(length(tempFrequency),1);
%                 end
%                 tempTuningWidth_kHz = funInvNErb(tempTuningWidth);
%
%
%                 eval(['tempR2 = data.' roiSaveName '.' groupName '.' analysisName_data '.r2;']);
%
%
%                 % seperate loop for scan data
%                 if analysisName_data == analysisNames{1} || analysisNames{2}
%
%                     for iBeta = 1:nCons(1)
%                         eval(['beta' mat2str(nCons(1)) ' = data.' roiSaveName '.' groupName '.' analysisName_data '.betas{iBeta,:}']);
%                     end
%
%                 else
%
%                     for iBeta = 1:nCons(1)
%                         eval(['betas' mat2str(nCons(1)) ' = nans(length(tempFrequency),1)']);
%                     end
%                 end
%
%                 nVoxels = length(tempFrequency);
%
%                 tempEstimateFreqName = repmat(estimateFreqNames{iEst},nVoxels,1);
%                 tempEstimateTuningName = repmat(estimateTuningName{iEst},nVoxels,1);
%                 tempAnalysis = repmat(analysisSaveName{iAnal},nVoxels,1);
%                 tempROI = repmat(roiName_data,nVoxels,1);
%
%                 ve_r2 = [ve_r2; tempR2'];
%
%                 Frequency = [Frequency; tempFrequency'];
%                 TuningWidth = [TuningWidth; tempTuningWidth'];
%                 Frequency_kHz = [Frequency; tempFrequency_kHz'];
%                 TuningWidth_kHz = [TuningWidth; tempTuningWidth_kHz'];
%
%                 if isempty(Analysis)
%                     Analysis = tempAnalysis;
%                     ROI = tempROI;
%                     estimateFreqName = tempEstimateFreqName;
%                     estimateTuningName = tempEstimateTuningName;
%                 else
%                     Analysis = char(Analysis,tempAnalysis);
%                     ROI = char(ROI,tempROI);
%                     estimateFreqName = char(estimateFreqName,tempEstimateFreqName);
%                     estimateTuningName = char(estimateTuningName,tempEstimateTuningName);
%                 end
%
%
%             end
%
%         end
%     end
% end
%
%
% subject = repmat(iSub,length(Frequency),1);
%
% T = table(Frequency,Frequency_kHz,...
%     TuningWidth,TuningWidth_kHz,...
%     ve_r2, Betas,...
%     Analysis,ROI,...
%     subject,...
%     'VariableNames',{'Frequency' 'Frequency_kHz' 'TuningWidth' 'TuningWidth_kHz' 'r2' 'Betas' 'Analysis' 'ROI','Subject'});
%
% writetable(T, [saveName, '_Comparisions.csv'])


%% Scan data







% end

%% get GLM data
% group

% scans


%% get pRF data

%% Estimating HRF
% average HRF estiamte


%% Compare Analysis



%% Compare acquistion
% Sparses Vs Continuous
% use best analysis
% Calculate
% Visualize


%% get hrf

%% notes from CM_subjectAnalysisScript
%% Tidy data
% cortical magnificaiton
%
%
% pRFrestrictROI = 'ARexp';
% pRFanalysisName = ['pRF_', pRFrestrictROI];
% analysisNames = {'glm_hrfDoubleGamma',pRFanalysisName};
% analysisSaveName = {'GLM','pRF'};
% AP = {'a','p'};
% analName = {'GLM', 'pRF'};
% CorticalDistance = [];
% Frequency = [];
% Analysis = [];
% ROI = [];
% ve_r2 = [];
% TuningWidth = [];
% for iSide = 1:length(Info.Sides)
%     for iGroup = 1:length(glmInfo.groupNames)
%         groupName = glmInfo.groupNames{iGroup};
%         for iAnal = 1:length(analysisNames)
%
%             analysisName_data = analysisNames{iAnal};
%             for iAP = 1:length(AP)
%
%                 roiSaveName = [Info.Sides{iSide}, 'GR' AP{iAP} '_' analName{iAnal}];
%                 roiName_data = [Info.Sides{iSide}, 'GR' AP{iAP}];
%
%                 eval(['tempCorticalDistance = data.' roiSaveName '.' groupName '.' analysisName_data '.tonotopicMagnificaion.relativeDistances(2,:);']);
%                 eval(['tempFrequency = data.' roiSaveName '.' groupName '.' analysisName_data '.tonotopicMagnificaion.pCF;']);
%                 eval(['tempFrequencycheck = data.' roiSaveName '.' groupName '.' analysisName_data '.tonotopicMagnificaion.pCFcheck{3};']);
%                 eval(['tempTuningWidth = data.' roiSaveName '.' groupName '.' analysisName_data '.tonotopicMagnificaion.pTW;']);
%                 eval(['tempR2 = data.' roiSaveName '.' groupName '.' analysisName_data '.tonotopicMagnificaion.r2;']);
%
%                 %             CorticalDistance =
%                 %             Frequency =
%
%                 %                 if tempFrequency == tempFrequencycheck
%                 nVoxels = length(tempFrequency);
%                 tempAnalysis = repmat(analysisSaveName{iAnal},nVoxels,1);
%                 tempROI = repmat(roiName_data,nVoxels,1);
%
%                 ve_r2 = [ve_r2; tempR2'];
%                 TuningWidth = [TuningWidth; tempTuningWidth'];
%                 CorticalDistance = [CorticalDistance; tempCorticalDistance'];
%                 Frequency = [Frequency; tempFrequency'];
%                 if isempty(Analysis)
%                     Analysis = tempAnalysis;
%                     ROI = tempROI;
%                 else
%                     Analysis = char(Analysis,tempAnalysis);
%                     ROI = char(ROI,tempROI);
%                 end
%
%
%             end
%
%         end
%     end
% end
% % Analysis = Analysis(2:end,:);
% % ROI = ROI(2:end,:);
%
% T = table(CorticalDistance,Frequency,...
%     ve_r2, TuningWidth,...
%     Analysis,ROI,...
%     'VariableNames',{'CorticalDistance' 'Frequency' 'r2' 'TuningWidth' 'Analysis' 'ROI'});
%
% writetable(T, [subjectInfo.subjectID, '_CM.csv'])

%% Comparisions
% what to compare?
% what extra things do I need from subjects

%% Comparisions: analysis
% convert units to match - convert overlays so figures are compareable
% convert stimulus to ERB for glm - pRF already convert - steal that code

%% Comparisions: aquistion
% use best analysis


% %% difference map
% % use averaged over depth overlays just created to make difference maps
%
% % make difference between groups maps - Sparse vs Continuous
% % make difference between analysis maps - GLM vs pRF
% % go to flatmap groups > take overlay from each group > subtrack them from
% % each other > install as new overlay
%
% for iSide = 1:length(Info.Sides)
%     overlay = cell(size(pRFInfo.analysisNames_Groups));
%     for iGroup = 1:length(glmInfo.groupNames)
%
%         thisView = viewSet(thisView,'curGroup',glmInfo.groupNames{iGroup});
%
%         % 'overlay'
%         %    overlay = viewGet(view,'overlay',[overlayNum],[analysisNum])
%         %    overlay = viewGet(view,'overlay',overlayNum,[])
%         %    overlay = viewGet(view,'overlay',[],analysisNum)
%         %    overlay = viewGet(view,'overlay',[],[])
%         %    overlay = viewGet(view,'overlay',overlayNum)
%         %    overlay = viewGet(view,'overlay')
%
%         %% loop over analysis
%         % add this in pRF analysis loop?
%         for iAnal = 1:length(pRFInfo.analysisNames_Groups{iGroup})
%             %         analysisName = pRFInfo.analysisNames_Groups{iGroup}{iAnal};
%             analysisName = [pRFInfo.analysisNames_Groups{iGroup}{iAnal}, '_', pRFInfo.pRFrois{iSide}, 'Vol' ];
%
%             thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',analysisName));
%             overlayNum = viewGet(thisView,'overlayNum','PrefCentreFreq');
%             overlay{iGroup}{iAnal} = viewGet(thisView,'overlay',overlayNum);
%         end
%     end
%     for iAnal = 1:length(pRFInfo.analysisNames_Groups{2})
%         %     analysisName = pRFInfo.analysisNames_Groups{iGroup}{iAnal};
%         analysisName = [pRFInfo.analysisNames_Groups{iGroup}{iAnal}, '_', pRFInfo.pRFrois{iSide}, 'Vol' ];
%
%         thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',analysisName));
%         [ thisView , differenceData ] = script_createDifferenceMaps(thisView,overlay{1}{1},overlay{2}{iAnal});
%     end
% end


% %% Save/export data for group average
% % all overlays = averaged across depths
% % pCF
% % pTW
% % Difference maps
% % PSIR - need to do the extra regression stuff - thickness and curvature
% % Gradient reversals - will need to perform analysis again and export to volume space
%
% %load MNI single subject flat sampling subject space
% ssMNI = {'Colin27','Colin27_flipX'};
% flatName = {'_34_122_117_Rad60','_141_153_99';'_138_126_114','_34_122_117'};
% mniRotation = [230,330;60,230];
% flatWarp = {'',['_invFNIRT_' subjects{iSubj}]};
% params.anatFileName = fullfile(dataDir,'Anatomy/freesurfer/subjects/', freeSurferName{iSubj}, 'surfRelax',[freeSurferName{iSubj} '_mprage_pp.nii']);
% params.flatRes=3;
% for iMNI =1:2
%   params.path = fullfile(dataDir,'Anatomy/freesurfer/subjects/',ssMNI{iMNI},'surfRelax');
%   for iWarp=1:2
%     for iSide=1:2
%         params.flatFileName = [ssMNI{iMNI} '_' sides{iSide} '_Flat' flatName{iMNI,iSide} '.off'];
%         params.outerCoordsFileName = [ssMNI{iMNI} '_' sides{iSide} '_GM_' freeSurferName{iSubj} flatWarp{iWarp} '.off'];
%         params.innerCoordsFileName = [ssMNI{iMNI} '_' sides{iSide} '_WM_' freeSurferName{iSubj} flatWarp{iWarp} '.off'];
%         params.curvFileName = [ssMNI{iMNI} '_' sides{iSide} '_Curv.vff'];
%         base = importFlatOFF(params);
%         base.name = [ssMNI{iMNI} '_' sides{iSide} '_Flat' flatWarp{iWarp}];
%         thisView = viewSet(thisView, 'newbase', base);
%         thisView = viewSet(thisView,'rotate',mniRotation(iMNI,iSide));
%     end
%   end
% end
%
% %%% anatomy
% thisView = viewSet(thisView,'curgroup',psirGroup);
% thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',psirAnalysis));
% thisView = viewSet(thisView,'curOverlay',1);
% for iSide = 1:2
%   %export data to Colin27
%   thisView = viewSet(thisView,'curbase',viewGet(thisView,'basenum',['Colin27_' sides{iSide} '_Flat']));
%   mrExport2SR(thisView.viewNum,fullfile(dataDir,studyDir,'flatExport',[subjects{iSubj} Sides{iSide} 'PSIR.nii']));
%   thisView = viewSet(thisView,'curbase',viewGet(thisView,'basenum',['Colin27_flipX_' sides{iSide} '_Flat']));
%   mrExport2SR(thisView.viewNum,fullfile(dataDir,studyDir,'flatExport',[subjects{iSubj} Sides{iSide} 'PSIRflipX.nii']));
% end
%
% %%% Tuning
% %set group
% thisView = viewSet(thisView,'curgroup',concatenationGroup);
% thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',functionalAnalysis));
% [thisView,params] = combineTransformOverlays(thisView,[],'justGetParams=1','defaultParams=1','overlayList=[2 3 4 5 6 7 8]');
% params.combineFunction='searchlightTuningWidth';
% params.additionalArgs = '[3 3 7],1';
% params.baseSpace=1;
% params.nOutputOverlays=7;
% for iSide = 1:2
%   thisView = viewSet(thisView,'curbase',viewGet(thisView,'basenum',[freeSurferName{iSubj} '_' sides{iSide} '_Flat_invFNIRT_' subjects{iSubj} '_lowres']));
%   params.outputName=['searchlightTuningWidth ' sides{iSide} ' flat lowres'];
%   [thisView,params] = combineTransformOverlays(thisView,params);
%   thisView = viewSet(thisView,'curOverlay',mainOverlays(iSubj,[1:10 10+iSide]));
%   thisView = viewSet(thisView,'curbase',viewGet(thisView,'basenum',['Colin27_' sides{iSide} '_Flat_invFNIRT_' subjects{iSubj}]));
%   mrExport2SR(thisView.viewNum,fullfile(dataDir,studyDir,'flatExport',[subjects{iSubj} Sides{iSide} 'Tonotopy.nii']));
%   thisView = viewSet(thisView,'curbase',viewGet(thisView,'basenum',['Colin27_flipX_' sides{iSide} '_Flat_invFNIRT_' subjects{iSubj}]));
%   mrExport2SR(thisView.viewNum,fullfile(dataDir,studyDir,'flatExport',[subjects{iSubj} Sides{iSide} 'TonotopyFlipX.nii']));
% end
%
%
% %% plot study information
% [ data ] = plot_studyInfo(stimInfo, glmInfo, pRFInfo, Info, plotInfo);


%% Comparisions
% what to compare?
% what extra things do I need from subjects

%% Comparisions: analysis
% convert units to match - convert overlays so figures are compareable
% convert stimulus to ERB for glm - pRF already convert - steal that code

%% Comparisions: aquistion
% use best analysis

%% Group analysis function will:
% import data
% tidy data
% statiscal analysis: average
% plot

%% perform ROI analysis
% NOTE: selecting data should happen outside of function
% roiAnalysis = script_ROIAnalysis(roiData,glmInfo.analysisBaseNames_Scans,Info,stimInfo,plotInfo,Info.conditionRunIndex,glmInfo.analysisScanNum,'GLM');
% for iSide = 1:length(Info.Sides)
%     eval(['ROInames = Info.' Info.Sides{iSide} 'ROInames;']);
%     for iROI = 1:length(ROInames)
%         eval(['roidata = data.' Info.Sides{iSide} '.' ROInames{iROI} ';']);
%         eval(['data.' Info.Sides{iSide} '.' ROInames{iROI} ' = script_ROIAnalysis(roidata,Info,glmInfo,stimInfo,plotInfo,subjectInfo,glmInfo.analysisScanNum,' q 'overlays' q ',ROInames{iROI});']);
%     end
% end


%     % delete view now we are done with it
%     deleteView(thisView);




%
%
%
%
%     % for iSub = 1:8
%
%     % tasks:
%     % get nVoxel data values
%     % repeat the properties nVoxel times
%     % save all in a table
%
%     subject = [];
%     tc_roi = [];
%     %  pCF estimate value
%     frequency_nERB = []; % estimate values in number of ERB
%     frequency_kHz = []; % estimate values in kHz
%     % pTW estimate value
%     ve_selectivity_nERB = []; % estimate values in number of ERB
%     ve_selectivity_kHz = []; % estimate values in kHz
%     tc_acquistion = []; % acquistion protocol
%     concatenation = []; % concatenated or individual scan
%     tc_analysis = []; % GLM or pRF
%     method = []; % for GLM {'Centriod'  'Spread'  'julien_pCF'  'julien_pTW'  'indexMax'} or pRF for pRF
%
%     % Concatenated data
%     Frequency = [];
%     Frequency_kHz = [];
%     TuningWidth = [];
%     TuningWidth_kHz = [];
%     Analysis = [];
%     ROI = [];
%     ve_r2 = [];
%     % voxelPropertyNames: {'Centriod'  'Spread'  'julien_pCF'  'julien_pTW'  'indexMax'}
%     % pRFOverlayNames: {'r2'  'PrefCentreFreq'  'rfHalfWidth'}
%     estimateFreqNames = {'Centriod', 'julien_pCF', 'indexMax'};
%     estimateTuningName = {'Spread', 'julien_pTW', 'NA'};
%
%     analysisType = {'GLM', 'pRF'};
%     hrfType = {'boxcar', 'Gamma'};
%
%
%     pRFanalysisName = ['pRF_', pRFInfo.pRFrestrictROI];
%     analysisNames = {'glm_hrfDoubleGamma', 'glm_hrfBoxCar', pRFanalysisName};
%     nCons = [32, 8];
%     nScans = 4;
%
%     for iSide = 1:length(Info.Sides)
%         for iGroup = 1:length(glmInfo.groupNames)
%             groupName = glmInfo.groupNames{iGroup};
%             for iAnal = 1:length(analysisType)
%
%                 analysisName_data = analysisNames{iAnal};
%                 if analysisName_data ~= analysisNames{3}
%                     estimateFreqNames = {'Centriod', 'julien_pCF', 'indexMax'};
%                     estimateTuningName = {'Spread', 'julien_pTW', 'NA'};
%                 else
%                     %%%%%%% FIND pRF save names %%%%%%%%%%%%%%%%
%                     %                     pRFInfo.pRFOverlayNames = {'r2','PrefCentreFreq','rfHalfWidth'};
%                     estimateFreqNames = {'pCF'};
%                     estimateTuningName = {'pTW'};
%                 end
%
%                 for iHRF = 1:length(hrfType)
%                     %                 roiSaveName = [Info.Sides{iSide}, 'GR_' analName{iAnal}];
%                     roiSaveName = [Info.Sides{iSide}, 'GR_GLM']; % compare using ROI data from derived GLM analysis
%                     roiName_data = [Info.Sides{iSide}, 'GR'];
%                     for iEst = 1:length(estimateFreqNames)
%                         eval(['tempFrequency = data.' roiSaveName '.' groupName '.' analysisName_data '.' estimateFreqNames{iEst}]);
%                         tempFrequency_kHz = funInvNErb(tempFrequency);
%
%                         if ~strcmp('NA',estimateTuningName{iEst})
%                             eval(['tempTuningWidth = data.' roiSaveName '.' groupName '.' analysisName_data '.' estimateTuningName{iEst}]);
%                         else
%                             tempTuningWidth = nan(length(tempFrequency),1);
%                         end
%                         tempTuningWidth_kHz = funInvNErb(tempTuningWidth);
%
%
%                         eval(['tempR2 = data.' roiSaveName '.' groupName '.' analysisName_data '.r2;']);
%
%
%                         % seperate loop for scan data
%                         if analysisName_data == analysisNames{1} || analysisNames{2}
%
%                             for iBeta = 1:nCons(1)
%                                 eval(['beta' mat2str(nCons(1)) ' = data.' roiSaveName '.' groupName '.' analysisName_data '.betas{iBeta,:}']);
%                             end
%
%                         else
%
%                             for iBeta = 1:nCons(1)
%                                 eval(['betas' mat2str(nCons(1)) ' = nans(length(tempFrequency),1)']);
%                             end
%                         end
%
%                         nVoxels = length(tempFrequency);
%
%                         tempEstimateFreqName = repmat(estimateFreqNames{iEst},nVoxels,1);
%                         tempEstimateTuningName = repmat(estimateTuningName{iEst},nVoxels,1);
%                         tempAnalysis = repmat(analysisSaveName{iAnal},nVoxels,1);
%                         tempROI = repmat(roiName_data,nVoxels,1);
%
%                         ve_r2 = [ve_r2; tempR2'];
%
%                         Frequency = [Frequency; tempFrequency'];
%                         TuningWidth = [TuningWidth; tempTuningWidth'];
%                         Frequency_kHz = [Frequency; tempFrequency_kHz'];
%                         TuningWidth_kHz = [TuningWidth; tempTuningWidth_kHz'];
%
%                         if isempty(Analysis)
%                             Analysis = tempAnalysis;
%                             ROI = tempROI;
%                             estimateFreqName = tempEstimateFreqName;
%                             estimateTuningName = tempEstimateTuningName;
%                         else
%                             Analysis = char(Analysis,tempAnalysis);
%                             ROI = char(ROI,tempROI);
%                             estimateFreqName = char(estimateFreqName,tempEstimateFreqName);
%                             estimateTuningName = char(estimateTuningName,tempEstimateTuningName);
%                         end
%
%                     end
%                 end
%
%             end
%         end
%     end
%
%
%     subject = repmat(iSub,length(Frequency),1);
%
%     T = table(Frequency,Frequency_kHz,...
%         TuningWidth,TuningWidth_kHz,...
%         ve_r2, Betas,...
%         Analysis,ROI,...
%         subject,...
%         'VariableNames',{'Frequency' 'Frequency_kHz' 'TuningWidth' 'TuningWidth_kHz' 'r2' 'Betas' 'Analysis' 'ROI','Subject'});
%
%     writetable(T, [saveName, '_Comparisions.csv'])
%
%
% Scan data

end
