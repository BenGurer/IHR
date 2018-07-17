function [thisView, pRFParams] = script_pRFAnalysis(thisView,pRFInfo,glmInfo,roiName,runSplitHalf,weightStim)
%
%   usage: script_pRFAnalysis(thisView,pRFInfo)
%      by: Ben Gurer
%    date: 22/11/2017
% purpose: script pRF analysis of fMRI data
%   input: mrView, pRF analysis information
%  output: updated mrView, pRF data
%
% NEED: ROI, save name, fit hdr?, weight stim?
% this function will be run after glm, gradient reversals and roi creation

% start parallel processing
nProcessors = mlrNumWorkers(4);

if ~weightStim
    for iGroup = 1:length(glmInfo.groupNames)
        thisView = viewSet(thisView,'curGroup',glmInfo.groupNames{iGroup});
        analysisSaveName = ['pRF' '_' roiName];
        [thisView, pRFParams] = pRF_auditory(thisView,[],'justGetParams=1','defaultParams=1');
        pRFParams.saveName = [analysisSaveName];
        pRFParams.restrict = ['ROI: ' roiName];
        pRFParams.pRFFit.supersampling = 1;
        pRFParams.pRFFit.fitHDR = 0;
        pRFParams.pRFFit.fwhm = 0;
        pRFParams.voxelScale = 'lin'; %'Scaling domain of voxel function.
        pRFParams.betaEachScan = true; %'Compute a separate beta weight (scaling) for each scan in the concanetation. This may be useful if there is some reason to believe that different scans have different magnitude responses, this will allow the fit to scale the magnitude for each scan'};
        pRFParams.algorithm = 'Levenberg-marquardt'; %'Which algorithm to use for optimization. Levenberg-marquardt seems to get stuck in local minimum, so the default is nelder-mead. However, levenberg-marquardt can set bounds for parameters, so may be better for when you are trying to fit the hdr along with the rf, since the hdr parameters can fly off to strange values.'};
        pRFParams.defaultConstraints = 0;
        if isfield(pRFInfo,'hrfParamsGamma')            
            % timelag = x(1);
            % tau = x(2);
            % exponent = x(3);            
            pRFParams.pRFFit.timelag = pRFInfo.hrfParamsGamma(1);
            pRFParams.pRFFit.tau = pRFInfo.hrfParamsGamma(2);
            pRFParams.pRFFit.exponent = pRFInfo.hrfParamsGamma(3);
        end
        %         if isfield(pRFInfo,'hrfParamsDiffofGamma')
        %             pRFParams.pRFFit.diffOfGamma = 1;
        %             % pRFParams find params for HRF
        %             % timelag = x(1);
        %             % tau = x(2);
        %             % exponent = x(3);
        %             % timelag2 = x(4);
        %             % tau2 = x(5);
        %             % exponent2 = x(6);
        %             % amplitude2 = x(7);
        %             pRFParams.pRFFit.timelag = pRFInfo.hrfParamsDiffofGamma(1);
        %             pRFParams.pRFFit.tau = pRFInfo.hrfParamsDiffofGamma(2);
        %             pRFParams.pRFFit.exponent = pRFInfo.hrfParamsDiffofGamma(3);
        %             pRFParams.pRFFit.timelag2 = pRFInfo.hrfParamsDiffofGamma(4);
        %             pRFParams.pRFFit.tau2 = pRFInfo.hrfParamsDiffofGamma(5);
        %             pRFParams.pRFFit.exponent2 = pRFInfo.hrfParamsDiffofGamma(6);
        %             pRFParams.pRFFit.amplitudeRatio = pRFInfo.hrfParamsDiffofGamma(7);
        %         end
        [thisView, pRFParams] = pRF_auditory(thisView,pRFParams);
        
        thisView = viewSet(thisView,'overlaycolorrange',[0 40],2);
        thisView = viewSet(thisView,'overlaycolorrange',[0 40],3);
        
        thisView = viewSet(thisView,'overlayrange',[0 40],2);
        thisView = viewSet(thisView,'overlayrange',[0 40],3);
        
        
        thisView = viewSet(thisView,'clipacrossoverlays',0);
        
        % save view
        mrSaveView(thisView);
    end
    
    if runSplitHalf
        % splitHRFmodel = hrfModel;
        thisView = viewSet(thisView,'curGroup',glmInfo.scanGroupName);
        for iScan = 1:glmInfo.nScans
            thisView = viewSet(thisView,'curScan', iScan);
            analysisSaveName = ['pRF' '_' roiName '_Scan_' num2str(iScan)];
            [thisView, pRFParams] = pRF_auditory(thisView,[],'justGetParams=1','defaultParams=1');
            pRFParams.saveName = [analysisSaveName];
            pRFParams.restrict = ['ROI: ' roiName];
            pRFParams.pRFFit.supersampling = 1;
            pRFParams.pRFFit.fitHDR = 0;
            pRFParams.pRFFit.fwhm = 0;
            pRFParams.scanNum = iScan;
            
            if isfield(pRFInfo,'hrfParamsGamma')                
                % timelag = x(1);
                % tau = x(2);
                % exponent = x(3);                
                pRFParams.pRFFit.timelag = pRFInfo.hrfParamsGamma(1);
                pRFParams.pRFFit.tau = pRFInfo.hrfParamsGamma(2);
                pRFParams.pRFFit.exponent = pRFInfo.hrfParamsGamma(3);
            end
            
            [thisView, pRFParams] = pRF_auditory(thisView,pRFParams);
            
            thisView = viewSet(thisView,'overlaycolorrange',[0 40],2);
            thisView = viewSet(thisView,'overlaycolorrange',[0 40],3);
            
            thisView = viewSet(thisView,'overlayrange',[0 40],2);
            thisView = viewSet(thisView,'overlayrange',[0 40],3);
            
            
            thisView = viewSet(thisView,'clipacrossoverlays',0);
            
            % save view
            mrSaveView(thisView);
        end
    end
    
else
    
    % pRFInfo.stimulusWeighting{1} = {'None'};
    % pRFInfo.stimulusWeighting{2} = {'None','SL_level','BOLD','fit'};
    
    
    for iGroup = 1: length(glmInfo.groupNames)
        thisView = viewSet(thisView,'curGroup',glmInfo.groupNames{iGroup});
        stimulusWeighting = pRFInfo.stimulusWeighting{iGroup};
        for iWeight = 1:length(stimulusWeighting)
            analysisSaveName = [ pRFInfo.analysisNames_Groups{iGroup}{iWeight} '_' roiName];
            [thisView, pRFParams] = pRF_auditory(thisView,[],'justGetParams=1','defaultParams=1');
            pRFParams.saveName = analysisSaveName;
            pRFParams.restrict = ['ROI: ' roiName];
            pRFParams.pRFFit.supersampling = 1;
            pRFParams.pRFFit.fitHDR = 0;
            pRFParams.pRFFit.fwhm = 0;
            pRFParams.pRFFit.stimulusWeighting = stimulusWeighting{iWeight}; % {'None','SL_level','BOLD','fit'}
            if strcmpi(stimulusWeighting{iWeight},'BOLD')
                pRFParams.pRFFit.SWgradient = glmInfo.m;
                pRFParams.pRFFit.SWoffset = glmInfo.b;
            end
            if isfield(pRFInfo,'hrfParamsGamma')                
                % timelag = x(1);
                % tau = x(2);
                % exponent = x(3);                
                pRFParams.pRFFit.timelag = pRFInfo.hrfParamsGamma(1);
                pRFParams.pRFFit.tau = pRFInfo.hrfParamsGamma(2);
                pRFParams.pRFFit.exponent = pRFInfo.hrfParamsGamma(3);
            end
            
            [thisView, pRFParams] = pRF_auditory(thisView,pRFParams);
            
            thisView = viewSet(thisView,'overlaycolorrange',[0 40],2);
            thisView = viewSet(thisView,'overlaycolorrange',[0 40],3);
            
            thisView = viewSet(thisView,'overlayrange',[0 40],2);
            thisView = viewSet(thisView,'overlayrange',[0 40],3);
            
            
            thisView = viewSet(thisView,'clipacrossoverlays',0);
            
            % save view
            mrSaveView(thisView);
            
        end
        
    end
    
    % run on scans without weighting
    if runSplitHalf
        % splitHRFmodel = hrfModel;
        thisView = viewSet(thisView,'curGroup',glmInfo.scanGroupName);
        for iScan = 1:glmInfo.nScans
            thisView = viewSet(thisView,'curScan', iScan);
            analysisSaveName = ['pRF' '_' roiName '_Scan - ' num2str(iScan)];
            [thisView, pRFParams] = pRF_auditory(thisView,[],'justGetParams=1','defaultParams=1');
            pRFParams.saveName = [analysisSaveName];
            pRFParams.restrict = ['ROI: ' roiName];
            pRFParams.pRFFit.supersampling = 1;
            pRFParams.pRFFit.fitHDR = 0;
            pRFParams.pRFFit.fwhm = 0;
            pRFParams.scanNum = iScan;
            if isfield(pRFInfo,'hrfParamsGamma')
                % timelag = x(1);
                % tau = x(2);
                % exponent = x(3);
                pRFParams.pRFFit.timelag = pRFInfo.hrfParamsGamma(1);
                pRFParams.pRFFit.tau = pRFInfo.hrfParamsGamma(2);
                pRFParams.pRFFit.exponent = pRFInfo.hrfParamsGamma(3);
            end
            
            [thisView, pRFParams] = pRF_auditory(thisView,pRFParams);
            
            thisView = viewSet(thisView,'overlaycolorrange',[0 40],2);
            thisView = viewSet(thisView,'overlaycolorrange',[0 40],3);
            
            thisView = viewSet(thisView,'overlayrange',[0 40],2);
            thisView = viewSet(thisView,'overlayrange',[0 40],3);
            
            
            thisView = viewSet(thisView,'clipacrossoverlays',0);
            
            % save view
            mrSaveView(thisView);
        end
    end
    
end

end