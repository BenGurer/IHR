
function thisView = script_hrfAnalysis(thisView,group)
%
%   usage: script_hrfAnalysis(thisView,glmInfo)
%      by: Ben Gurer
%    date: 13/03/2017
% purpose: script hrf analysis. Uses glm deconvolution to estimate hrf.
%   input: thisView, glmInfo
%  output: thisView
%

% set to continouus conca group

% thisView = viewSet(thisView,'curGroup','ConcatenationCont');
thisView = viewSet(thisView,'curGroup',group);
% thisView = viewSet(thisView,'curGroup','Concatenation Cont');
refreshMLRDisplay(thisView);

[thisView, glmParams] = glmAnalysis(thisView,[],'justGetParams=1','defaultParams=1');
glmParams.hrfModel = 'hrfDeconvolution';
[thisView, glmParams] = glmAnalysis(thisView,glmParams,'justGetParams=1','defaultParams=1');
glmParams.scanParams{1, 1}.preprocess  = 'binStimFreq';
glmParams.hrfParams.hdrlenS = 15;
glmParams.numberContrasts = 8;
glmParams.componentsToTest = [0 1 1 1 1 0 0 0 0 0];
glmParams.numberEVs = 8;
glmParams.computeTtests = 1;
glmParams.numberFtests  = 1;
glmParams.fTestNames{1, 1} = 'fTest - all conditions';
glmParams.restrictions{1, 1} = [1,0,0,0,0,0,0,0;0,1,0,0,0,0,0,0;0,0,1,0,0,0,0,0;0,0,0,1,0,0,0,0;0,0,0,0,0,1,0,0;0,0,0,0,0,0,1,0;0,0,0,0,0,0,0,1];
glmParams.alphaContrastOverlay = 'Uncorrected';
glmParams.parametricTests = 1;
glmParams.fweAdjustment = 0;
glmParams.fdrAdjustment = 0;
glmParams.outputStatistic = 0;
[thisView, glmParams] = glmAnalysis(thisView,glmParams,'justGetParams=1','defaultParams=1');
glmParams.saveName = 'GLM_Deconv_8bins';
glmParams.hrfParams.description = 'GLM Deconvolution';

[thisView, glmParams] = glmAnalysis(thisView,glmParams);
