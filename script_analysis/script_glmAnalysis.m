function [thisView] = script_glmAnalysis(thisView,glmInfo,groupNames,hrfModel,runSplitHalf,runTonotopic)
%
%   usage: script_glmAnalysis(thisView,glmInfo)
%      by: Ben Gurer
%    date: 13/03/2017
% purpose: script glm analysis for tonotopy
%   input: thisView, glmInfo
%  output: thisView
%
% hrfModel = glmInfo.hrfModel;
for i = 1:length(groupNames)
    for iHRF = 1:length(hrfModel)
        analysisName = ['glm_' hrfModel{iHRF}];
        thisView = viewSet(thisView,'curGroup',groupNames{i});
        [thisView, glmParams] = glmAnalysis(thisView,[],'justGetParams=1','defaultParams=1');
        glmParams.hrfModel = hrfModel{iHRF};
        [thisView, glmParams] = glmAnalysis(thisView,glmParams,'justGetParams=1','defaultParams=1');
        glmParams.saveName = analysisName;
        glmParams.hrfParams.description = hrfModel{iHRF};
        switch hrfModel{iHRF}
            case 'hrfBoxcar'
                %                 switch glmInfo.studyNames
                %                     case 'CM'
                %                         switch glmInfo.groupNames{i}
                %                             case glmInfo.groupNames{1}
                %                                 glmParams.hrfParams.delayS =  2.5;
                %                                 glmParams.hrfParams.durationS = 2.5;
                %                             case glmInfo.groupNames{2}
                %                                 glmParams.hrfParams.delayS =  2.5;
                %                                 glmParams.hrfParams.durationS = 2.5;
                %                         end
                %                     case 'sHL'
                glmParams.hrfParams.delayS =  2.5;
                glmParams.hrfParams.durationS = 2.5;
                %                 end
            case 'hrfDoubleGamma'
                if isfield(glmInfo,'hrfParamsDoubleGamma')
                    glmParams.hrfParams.x = glmInfo.hrfParamsDoubleGamma(1);
                    glmParams.hrfParams.y = glmInfo.hrfParamsDoubleGamma(2);
                    glmParams.hrfParams.z = glmInfo.hrfParamsDoubleGamma(3);
                else
                    glmParams.hrfParams.x =  4;
                    glmParams.hrfParams.y = 11;
                    glmParams.hrfParams.z =  4;
                end
        end
        glmParams.scanParams{1}.stimDurationMode = 'From file';
        glmParams.scanParams{1}.supersamplingMode =  'Set value';
        glmParams.scanParams{1}.designSupersampling = 3;
        glmParams.scanParams{1}.acquisitionDelay = .75;
        glmParams.computeTtests = 1;
        glmParams.numberContrasts  = 1;
        glmParams.contrasts = ones(1,32);
        glmParams.numberFtests  = 1;
        glmParams.fTestNames{1, 1} = 'fTest - all conditions';
        glmParams.restrictions{1, 1} = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
            0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
            0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
            0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
            0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
            0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
            0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
            0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
            0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
            0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
            0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
            0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
            0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
            0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0;...
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0;...
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0;...
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0;...
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0;...
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0;...
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0;...
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0;...
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0;...
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0;...
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0;...
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0;...
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1];
        glmParams.alphaContrastOverlay = 'Uncorrected';
        glmParams.parametricTests = 1;
        glmParams.fweAdjustment = 0;
        glmParams.fdrAdjustment = 1;
        glmParams.outputStatistic = 1;
        glmParams.numberContrasts = 0;
        glmParams.outputEstimatesAsOverlays = 1;
        [thisView, glmParams] = glmAnalysis(thisView,glmParams);
        
        if runTonotopic
            %Tonotopy analysis
            % Index max
            [thisView,params] = combineTransformOverlays(thisView,[],'justGetParams=1','defaultParams=1',['overlayList=' mat2str([2:33])]);
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
            thisView = viewSet(thisView,'overlaycolorrange',[0 32],curOverlay-3);
            thisView = viewSet(thisView,'overlaycolorrange',[0 32],curOverlay-2);
            thisView = viewSet(thisView,'overlaycolorrange',[-2 38],curOverlay-1);
            thisView = viewSet(thisView,'overlaycolorrange',[-2 38],curOverlay);
            
            %% save analysis
            saveAnalysis(thisView,analysisName)
            
        end
        
    end
end

%% Split Run Analysis
%% GLM Analysis
% nStim = [8 32];
% Set hrf type based on acquisition type
if runSplitHalf
    % splitHRFmodel = hrfModel;
    thisView = viewSet(thisView,'curGroup',glmInfo.scanGroupName);
    for iScan = 1:glmInfo.nScans
        thisView = viewSet(thisView,'curScan', iScan);
        for iHRF = 1:length(hrfModel)
            for iStim = 1:length(glmInfo.nStim)
                analysisName_split{iScan} = ['glm_' hrfModel{iHRF} '_nCons_' mat2str(glmInfo.nStim(iStim)) '_Scan_' mat2str(iScan)];
                [thisView, glmParams] = glmAnalysis(thisView,[],'justGetParams=1','defaultParams=1',['scanList=' mat2str(iScan)]);
                glmParams.hrfModel = hrfModel{iHRF};
                [thisView, glmParams] = glmAnalysis(thisView,glmParams,'justGetParams=1','defaultParams=1',['scanList=' mat2str(iScan)]);
                %         glmParams.saveName = analysisName;
                glmParams.hrfParams.description = hrfModel{iHRF};
                switch hrfModel{iHRF}
                    case 'hrfBoxcar'
                        glmParams.hrfParams.delayS =  2.5;
                        glmParams.hrfParams.durationS = 2.5;
                    case 'hrfDoubleGamma'
                        if isfield(glmInfo,'hrfParamsDoubleGamma')
                            glmParams.hrfParams.x = glmInfo.hrfParamsDoubleGamma(1);
                            glmParams.hrfParams.y = glmInfo.hrfParamsDoubleGamma(2);
                            glmParams.hrfParams.z = glmInfo.hrfParamsDoubleGamma(3);
                        else
                            glmParams.hrfParams.x =  4;
                            glmParams.hrfParams.y = 11;
                            glmParams.hrfParams.z =  4;
                        end
                end
                
                if glmInfo.nStim(iStim) == 8
                    glmParams.scanParams{1, iScan}.preprocess  = 'binStimFreq';
                    glmParams.EVnames = {'1','2','3','4','5','6','7','8'};
                end
                glmParams.numberContrasts = glmInfo.nStim(iStim);
                glmParams.numberEVs = glmInfo.nStim(iStim);
                [thisView, glmParams] = glmAnalysis(thisView,glmParams,'justGetParams=1','defaultParams=1',['scanList=' mat2str(iScan)]);
                glmParams.scanParams{iScan}.stimDurationMode = 'From File';
                glmParams.scanParams{iScan}.supersamplingMode =  'Set value';
                glmParams.scanParams{iScan}.designSupersampling = 3;
                glmParams.scanParams{iScan}.acquisitionDelay = .75;
                glmParams.numberContrasts = 0;
                glmParams.parametricTests = 0;
                glmParams.outputEstimatesAsOverlays = 1;
                glmParams.saveName = analysisName_split{iScan};
                glmParams.scanNum = iScan;
                [thisView, glmParams] = glmAnalysis(thisView,glmParams,['scanList=' mat2str(iScan)]);
                
                %                 switch hrfModel{iHRF}
                %                     case 'hrfDoubleGamma'
                %                         if glmInfo.nStim(iStim) == 32
                if runTonotopic
                    %Tonotopy analysis
                    % Index max
                    if glmInfo.nStim(iStim) == 32
                        [thisView,params] = combineTransformOverlays(thisView,[],'justGetParams=1','defaultParams=1',['overlayList=' mat2str([2:33])]);
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
                        thisView = viewSet(thisView,'overlaycolorrange',[0 32],curOverlay-3);
                        thisView = viewSet(thisView,'overlaycolorrange',[0 32],curOverlay-2);
                        thisView = viewSet(thisView,'overlaycolorrange',[0 40],curOverlay-1);
                        thisView = viewSet(thisView,'overlaycolorrange',[0 40],curOverlay);
                    else
                        
                        [thisView,params] = combineTransformOverlays(thisView,[],'justGetParams=1','defaultParams=1',['overlayList=' mat2str([2:9])]);
                        params.combineFunction='indexMax';
                        params.nOutputOverlays=2;
                        [thisView,params] = combineTransformOverlays(thisView,params);
                        curOverlay=viewGet(thisView,'curOverlay');
                        thisView = viewSet(thisView,'overlaycolorrange',[0 8],curOverlay-1);
                        
                        % Weighted mean and corrected weighted mean
                        params.combineFunction='weightedMeanStd_CM';
                        params.nOutputOverlays=4;
                        [thisView,params] = combineTransformOverlays(thisView,params);
                        curOverlay=viewGet(thisView,'curOverlay');
                        thisView = viewSet(thisView,'overlaycolorrange',[0 8],curOverlay-3);
                        thisView = viewSet(thisView,'overlaycolorrange',[0 8],curOverlay-2);
                        thisView = viewSet(thisView,'overlaycolorrange',[0 10],curOverlay-1);
                        thisView = viewSet(thisView,'overlaycolorrange',[0 10],curOverlay);
                    end
                    %% save analysis
                        saveAnalysis(thisView,analysisName_split{iScan})
                end
                
                %                         end
                %                 end
                
            end
        end
    end
    
end
% save view and quit
mrSaveView(thisView);
