%% pilot_subjectAnalysisScript
doLoadView = 1;
doDeconv = 1;

%% Setting up subject
% Get subject info
disp('Getting subject info...')
subjectInfo = get_SubjectInfo_sHL(iSubs2Run(iSub));
% Subject ID, flatmap names etc
pRFInfo.sHLscans = subjectInfo.conditionOrder{2};

% define file name for data
saveName = [subjectInfo.subjectID '_data.mat'];

% move to subject folder
cd(fullfile(Info.dataDir,Info.studyDir,subjectInfo.subjectID));

% load data if it exists
if exist(saveName, 'file')
    disp('Loading data...')
    load(saveName)
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

%% GLM Analysis with deconvolved HRF
% Estiamting hemodynamic response function
if doDeconv
    % estimate hrf using deconvolution
    disp('Performing GLM analysis (HRF=deconvolution)')
    % get view so we have the ROIs
    thisView = getMLRView;
    
    % hrf
    thisView = viewSet(thisView,'curGroup','ConcatenationContHDR');
    refreshMLRDisplay(thisView);
    
    [thisView, glmParams] = glmAnalysis(thisView,[],'justGetParams=1','defaultParams=1');
    glmParams.hrfModel = 'hrfDeconvolution';
    [thisView, glmParams] = glmAnalysis(thisView,glmParams,'justGetParams=1','defaultParams=1');
    glmParams.numberContrasts = 1;
    glmParams.hrfParams.hdrlenS = 16;
    glmParams.scanParams{1}.stimDurationMode = 'From file';
    glmParams.scanParams{1}.supersamplingMode =  'Set value';
    glmParams.scanParams{1}.designSupersampling = 3;
    glmParams.scanParams{1}.acquisitionDelay = .75;
    glmParams.contrasts = [1;1;1;1];
    glmParams.componentsToTest = [0 1 1 1 1 0 0 0];
    glmParams.numberEVs = 1;
    glmParams.computeTtests = 1;
    glmParams.alphaContrastOverlay = 'Uncorrected';
    glmParams.parametricTests = 1;
    glmParams.fweAdjustment = 0;
    glmParams.fdrAdjustment = 0;
    glmParams.outputStatistic = 0;
    [thisView, glmParams] = glmAnalysis(thisView,glmParams,'justGetParams=1','defaultParams=1');
    glmParams.saveName = 'GLM_Deconvolution';
    glmParams.hrfParams.description = 'GLM Deconvolution';
    
    [thisView, glmParams] = glmAnalysis(thisView,glmParams);
    
    thisView = getMLRView;
    % get data from analysis
    analysisData = get_analysisData(thisView,'GLM_Deconvolution');
    
    roi = viewGet(thisView,'roi','ARHG');
    
    % restrict data by ROI
    roi.scanCoords = getROICoordinates(thisView,roi);
    % get ROI estimates
    r2data = analysisData.overlays(1).data;
    volumeIndices = sub2ind(size(r2data{:}),roi.scanCoords(1,:),roi.scanCoords(2,:),roi.scanCoords(3,:));
    [estimate,volumeIndices] = getEstimates(analysisData.d{:}, analysisData.params ,volumeIndices');
    % nVoxels = length(volumeIndices);
    
    e_mean = mean(estimate.hdr,3);
    e_sum = nansum(estimate.hdr,3);
    e_sum_norm = e_sum./max(e_sum);
    e_std = std(estimate.hdr,0,3);
    %     eval(['data.hrf.' hrfGroupName{iGroup} '.estimate = e_sum_norm;']);
    data.hrf.noise = e_sum_norm';
    
    % tono
    
    thisView = viewSet(thisView,'curGroup','ConcatenationContpRF');
    thisView = getMLRView;
    refreshMLRDisplay(thisView);
    
    [thisView, glmParams] = glmAnalysis(thisView,[],'justGetParams=1','defaultParams=1');
    glmParams.hrfModel = 'hrfDeconvolution';
    glmParams.scanParams{1, 1}.preprocess  = 'binStim_pilot_hrf';
    glmParams.EVnames = [];
    [thisView, glmParams] = glmAnalysis(thisView,glmParams,'justGetParams=1','defaultParams=1');
    glmParams.hrfParams.hdrlenS = 16;
    glmParams.scanParams{1}.stimDurationMode = 'From file';
    glmParams.scanParams{1}.supersamplingMode =  'Set value';
    glmParams.scanParams{1}.designSupersampling = 3;
    glmParams.scanParams{1}.acquisitionDelay = .75;
    glmParams.numberContrasts = 10;
    glmParams.componentsToTest = [0 1 1 1 1 0 0 0];
    glmParams.numberEVs = 10;
    glmParams.computeTtests = 1;
    glmParams.alphaContrastOverlay = 'Uncorrected';
    glmParams.parametricTests = 1;
    glmParams.fweAdjustment = 0;
    glmParams.fdrAdjustment = 0;
    glmParams.outputStatistic = 0;
    [thisView, glmParams] = glmAnalysis(thisView,glmParams,'justGetParams=1','defaultParams=1');
    glmParams.saveName = 'GLM_Deconvolution';
    glmParams.hrfParams.description = 'GLM Deconvolution';
    
    [thisView, glmParams] = glmAnalysis(thisView,glmParams);
    
    thisView = getMLRView;
    % get data from analysis
    analysisData = get_analysisData(thisView,'GLM_Deconvolution');
    
    roi = viewGet(thisView,'roi','ARHG');
    
    % restrict data by ROI
    roi.scanCoords = getROICoordinates(thisView,roi);
    % get ROI estimates
    r2data = analysisData.overlays(1).data;
    volumeIndices = sub2ind(size(r2data{:}),roi.scanCoords(1,:),roi.scanCoords(2,:),roi.scanCoords(3,:));
    [estimate,volumeIndices] = getEstimates(analysisData.d{:}, analysisData.params ,volumeIndices');
    
    e = estimate.hdr;
    [ x_doubleGamma, x_Gamma, x_dGamma, hrf_deconv, hrf_deconvTW] = cal_hrfROIAverage(e,estimate.time,analysisData.d{:});
    
    e_norm = hrf_deconv./ max(hrf_deconv);
    data.hrf.tono = e_norm;
    
    data.hrf.time = estimate.time;
    
    disp('saving HRF estimate data...')
    save('hrf','data','-v7.3');
end

% tidy HRF estimates
tidy_hrf = [data.hrf.tono'; data.hrf.noise'];
tidy_time = [data.hrf.time'; data.hrf.time'];
tidy_names = [string(repmat('tono',length(data.hrf.tono),1));...
    string(repmat('noise',length(data.hrf.noise),1))];

hrf_table = [];
hrf_table = table(tidy_hrf,tidy_time,tidy_names,...
    'VariableNames',{'HRF', 'Time', 'Name'});

writetable(hrf_table, 'HRF_est.csv')

%% tonotopic estiamtes
%% GLM
thisView = viewSet(thisView,'curGroup','ConcatenationContpRF');

[thisView, glmParams] = glmAnalysis(thisView,[],'justGetParams=1','defaultParams=1');
glmParams.scanParams{1, 1}.preprocess  = 'binStim_pilot_tono';
glmParams.hrfModel = 'hrfDoubleGamma';
glmParams.EVnames = [];
[thisView, glmParams] = glmAnalysis(thisView,glmParams,'justGetParams=1','defaultParams=1');
glmParams.saveName = 'GLM_hrfDoubleGamma';
glmParams.hrfParams.description = 'hrfDoubleGamma';
glmParams.hrfParams.x =  4;
glmParams.hrfParams.y = 11;
glmParams.hrfParams.z =  4;
glmParams.scanParams{1}.stimDurationMode = 'From file';
glmParams.scanParams{1}.supersamplingMode =  'Set value';
glmParams.scanParams{1}.designSupersampling = 3;
glmParams.scanParams{1}.acquisitionDelay = .75;
glmParams.numberFtests  = 1;
glmParams.fTestNames{1, 1} = 'fTest - all conditions';
glmParams.restrictions{1, 1} = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
    0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
    0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
    0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
    0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
    0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
    0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0;...
    0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0;...
    0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0;...
    0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0;...
    0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0;...
    0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0;...
    0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0;...
    0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0;...
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0;...
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0;...
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0;...
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0;...
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0;...
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1;...
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
glmParams.alphaContrastOverlay = 'Uncorrected';
glmParams.parametricTests = 1;
glmParams.fweAdjustment = 0;
glmParams.fdrAdjustment = 1;
glmParams.outputStatistic = 1;
glmParams.numberContrasts = 0;
glmParams.outputEstimatesAsOverlays = 1;
[thisView, glmParams] = glmAnalysis(thisView,glmParams);

[thisView,params] = combineTransformOverlays(thisView,[],'justGetParams=1','defaultParams=1',['overlayList=' mat2str([2:21])]);
params.combineFunction='indexMax';
params.nOutputOverlays=2;
[thisView,params] = combineTransformOverlays(thisView,params);
curOverlay=viewGet(thisView,'curOverlay');
thisView = viewSet(thisView,'overlaycolorrange',[0 32],curOverlay-1);

% Weighted mean and corrected weighted mean
params.combineFunction='weightedMeanStd_CM';
params.nOutputOverlays=4;
[thisView,params] = combineTransformOverlays(thisView,params);
curOverlay=viewGet(thisView,'curOverlay');
thisView = viewSet(thisView,'overlaycolorrange',[0 20],curOverlay-3);
thisView = viewSet(thisView,'overlaycolorrange',[0 20],curOverlay-2);
thisView = viewSet(thisView,'overlaycolorrange',[0 28],curOverlay-1);
thisView = viewSet(thisView,'overlaycolorrange',[0 28],curOverlay);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thisView = viewSet(thisView,'curGroup','ConcatenationSparsepRF');

[thisView, glmParams] = glmAnalysis(thisView,[],'justGetParams=1','defaultParams=1');
glmParams.scanParams{1, 1}.preprocess  = 'binStim_pilot_tono';
glmParams.hrfModel = 'hrfDoubleGamma';
glmParams.EVnames = [];
[thisView, glmParams] = glmAnalysis(thisView,glmParams,'justGetParams=1','defaultParams=1');
glmParams.saveName = 'GLM_hrfDoubleGamma';
glmParams.hrfParams.description = 'hrfDoubleGamma';
glmParams.hrfParams.x =  4;
glmParams.hrfParams.y = 11;
glmParams.hrfParams.z =  4;
glmParams.scanParams{1}.stimDurationMode = 'From file';
glmParams.scanParams{1}.supersamplingMode =  'Set value';
glmParams.scanParams{1}.designSupersampling = 3;
glmParams.scanParams{1}.acquisitionDelay = .75;

glmParams.alphaContrastOverlay = 'Uncorrected';
glmParams.parametricTests = 1;
glmParams.fweAdjustment = 0;
glmParams.fdrAdjustment = 1;
glmParams.outputStatistic = 1;
glmParams.numberContrasts = 0;
glmParams.outputEstimatesAsOverlays = 1;
[thisView, glmParams] = glmAnalysis(thisView,glmParams);

[thisView,params] = combineTransformOverlays(thisView,[],'justGetParams=1','defaultParams=1',['overlayList=' mat2str([2:16])]);
params.combineFunction='indexMax';
params.nOutputOverlays=2;
[thisView,params] = combineTransformOverlays(thisView,params);
curOverlay=viewGet(thisView,'curOverlay');
thisView = viewSet(thisView,'overlaycolorrange',[0 32],curOverlay-1);

% Weighted mean and corrected weighted mean
params.combineFunction='weightedMeanStd_CM';
params.nOutputOverlays=4;
[thisView,params] = combineTransformOverlays(thisView,params);
curOverlay=viewGet(thisView,'curOverlay');
thisView = viewSet(thisView,'overlaycolorrange',[0 16],curOverlay-3);
thisView = viewSet(thisView,'overlaycolorrange',[0 16],curOverlay-2);
thisView = viewSet(thisView,'overlaycolorrange',[0 24],curOverlay-1);
thisView = viewSet(thisView,'overlaycolorrange',[0 24],curOverlay);

%% pRF
tonoGroupName = {'ConcatenationSparsepRF', 'ConcatenationContpRF'};
roiName = 'AR_exp';
for iGroup = 1:length(tonoGroupName)
    thisView = viewSet(thisView,'curGroup',tonoGroupName{iGroup});
    analysisSaveName = ['pRF' '_' roiName];
    [thisView, pRFParams] = pRF_auditory(thisView,[],'justGetParams=1','defaultParams=1');
    pRFParams.saveName = analysisSaveName;
    pRFParams.restrict = ['ROI: ' roiName];
    pRFParams.pRFFit.supersampling = 1;
    pRFParams.pRFFit.fitHDR = 0;
    pRFParams.pRFFit.fwhm = 0;
    pRFParams.voxelScale = 'lin'; %'Scaling domain of voxel function.
    pRFParams.betaEachScan = true; %'Compute a separate beta weight (scaling) for each scan in the concanetation. This may be useful if there is some reason to believe that different scans have different magnitude responses, this will allow the fit to scale the magnitude for each scan'};
    pRFParams.algorithm = 'Levenberg-marquardt'; %'Which algorithm to use for optimization. Levenberg-marquardt seems to get stuck in local minimum, so the default is nelder-mead. However, levenberg-marquardt can set bounds for parameters, so may be better for when you are trying to fit the hdr along with the rf, since the hdr parameters can fly off to strange values.'};
    pRFParams.defaultConstraints = 0;
    pRFParams.pRFFit.diffOfGamma = 0;
    [thisView, pRFParams] = pRF_auditory(thisView,pRFParams);
    
    thisView = viewSet(thisView,'overlaycolorrange',[0 40],2);
    thisView = viewSet(thisView,'overlaycolorrange',[0 40],3);
    
    thisView = viewSet(thisView,'overlayrange',[0 40],2);
    thisView = viewSet(thisView,'overlayrange',[0 40],3);
    
    thisView = viewSet(thisView,'clipacrossoverlays',0);
end

    if doPrintFigures
        
        thisView = viewSet(thisView,'showrois','hide');        
        thisView = viewSet(thisView,'basecorticaldepth', [0.3, 0.7]);
        thisView = viewSet(thisView,'alpha',1);
        
        glmoverlayNames = {'r2', [glmInfo.voxelPropertyNames{1} '_nERB'], [glmInfo.voxelPropertyNames{2} '_nERB'],...
            [glmInfo.voxelPropertyNames{3} '_nERB'], [glmInfo.voxelPropertyNames{4} '_nERB'], [glmInfo.voxelPropertyNames{5} '_nERB']};
        pRFoverlayNames = pRFInfo.pRFOverlayNames;
        % names used in saved data structure
        estimateNames_pRF = {'r2', 'pCF', 'pTW'};
        estimateNames_GLM = {'r2', 'pCFc', 'pTWs', 'pCFj', 'pTWj', 'pCFi'};
        
        pRFrestrictROI = 'ARexp';
        pRFanalysisName = ['pRF_', pRFrestrictROI];
        analysisNames = {'glm_hrfDoubleGamma',pRFanalysisName};
        analysisSave_Names = {'GLM','pRF'};
        
        roiNames = {'GR' 'GRa' 'GRp'};
        AP = {'a','p'};
        %         analName = {'GLM', 'pRF'};
        roiAnalName = {'GLM'};
        
        groupSaveNames = {'sparse','cont'};
        
        for iSide = 1:length(Info.Sides)
            if viewGet(thisView,'currentbase') ~= viewGet(thisView,'baseNum',[subjectInfo.flatmapNames{iSide}])
                baseNum = viewGet(thisView,'baseNum',[subjectInfo.flatmapNames{iSide}]);
                thisView = viewSet(thisView,'currentbase',baseNum);
            end
            
            for iGroup = 1:length(glmInfo.groupNames)
                if viewGet(thisView,'curgroup') ~= viewGet(thisView,'groupNum',glmInfo.groupNames{iGroup})
                    thisView = viewSet(thisView,'curgroup',glmInfo.groupNames{iGroup});
                end
                
                for iAnal = 1:length(analysisNames)
                    analysisName = analysisNames{iAnal};
                    
                    if viewGet(thisView,'curAnalysis') ~= viewGet(thisView,'analysisNum',analysisName)
                        thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',analysisName));
                    end                    
                    
                    switch analysisName
                        case 'glm_hrfDoubleGamma'
                            overlayNames = glmoverlayNames;
                            estimateNames = estimateNames_GLM;
                        case pRFanalysisName
                            overlayNames = pRFoverlayNames;
                            estimateNames = estimateNames_pRF;
                    end
                    
                    for iOverlay = 1:length(overlayNames)
                        
                        overlayNum = viewGet(thisView,'overlayNum',overlayNames{iOverlay});
                        thisView = viewSet(thisView,'curOverlay',overlayNum);
                        
                        saveName = [groupSaveNames{iGroup} '_' analysisSave_Names{iAnal} '_' estimateNames{iOverlay} '_' Info.sides{iSide}];
                        
                        refreshMLRDisplay(thisView);
                        
                        mrPrint(thisView,'useDefault=1','roiSmooth=0','roiLabels=0')
                        fh = findobj( 'Type', 'Figure', 'Name', 'Print figure' );
                        saveas(fh,saveName,'svg')
                        close(fh)
                        
                    end
                end
            end
        end
        
        thisView = viewSet(thisView,'alpha',0);
        thisView = viewSet(thisView,'currentroi',viewGet(thisView,'roiNum','AR'));
        
        for iSide = 1:length(Info.Sides)
            if viewGet(thisView,'currentbase') ~= viewGet(thisView,'baseNum',[subjectInfo.flatmapNames{iSide}])
                baseNum = viewGet(thisView,'baseNum',[subjectInfo.flatmapNames{iSide}]);
                thisView = viewSet(thisView,'currentbase',baseNum);
            end
            refreshMLRDisplay(thisView);
            
            % anatomy
            thisView = viewSet(thisView,'showrois','hide');
            saveName = ['anat_', Info.sides{iSide}];
            mrPrint(thisView,'useDefault=1','roiSmooth=0','roiLabels=0')
            fh = findobj( 'Type', 'Figure', 'Name', 'Print figure' );
            saveas(fh,saveName,'svg')
            close(fh)
            
            % AR roi
            thisView = viewSet(thisView,'showrois','selected perimeter');
            saveName = ['roi_AR_', Info.sides{iSide}];
            mrPrint(thisView,'useDefault=1','roiSmooth=0','roiLabels=0')
            fh = findobj( 'Type', 'Figure', 'Name', 'Print figure' );
            saveas(fh,saveName,'svg')
            close(fh)
            
        end
        
        % GR rois
        thisView = viewSet(thisView,'showrois','selected perimeter');
        for iSide = 1:length(Info.Sides)
            if viewGet(thisView,'currentbase') ~= viewGet(thisView,'baseNum',[subjectInfo.flatmapNames{iSide},'Volume'])
                baseNum = viewGet(thisView,'baseNum',[subjectInfo.flatmapNames{iSide},'Volume']);
                thisView = viewSet(thisView,'currentbase',baseNum);
            end
            
            for iROI = 1:length(roiNames)
                
                saveName = ['roi_', roiNames{iROI}, Info.sides{iSide}];
                
                thisView = viewSet(thisView,'currentroi',viewGet(thisView,'roiNum',[Info.Sides{iSide},roiNames{iROI}, '_GLM']));
                refreshMLRDisplay(thisView);
                
                mrPrint(thisView,'useDefault=1','roiSmooth=0','roiLabels=0')
                fh = findobj( 'Type', 'Figure', 'Name', 'Print figure' );
                saveas(fh,saveName,'svg')
                close(fh)
            end
            
        end
        
        thisView = viewSet(thisView,'alpha',1);
    end



