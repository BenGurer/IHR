% methods info

lowFreqkHz = 0.1;
highFreqkHz = 8;
nStim = 32;

%% Get stimulus properties
% get stimulus frequencies
stimInfo.lowFreqkHz = 0.1;
stimInfo.highFreqkHz = 8;
stimInfo.nStim = 32;

[stimInfo.stimFreqs, stimInfo.stimFreqs_bin, stimInfo.stimFreqs_mv, stimInfo.stimNERBs, stimInfo.stimNERBs_bin, stimInfo.stimNERBs_mv] = convertStimIDtoFrequency(stimInfo.lowFreqkHz,stimInfo.highFreqkHz,stimInfo.nStim);

stimOne = stimInfo.stimNERBs(1);
stimTwo = stimInfo.stimNERBs(2);
% 
% [ thisView , ~ ] = convertOverlay_GLMCF2NERB(thisView,overlayIN,stimOne,stimTwo,[glmInfo.voxelPropertyNames{i} '_nERB'],debaised);


minERB = -2.* (stimTwo-stimOne) + (stimOne - 1);

maxERB = 38* (stimTwo-stimOne) + (stimOne - 1);

maxkHZ = cal_nERB2kHz(maxERB)
minkHz = cal_nERB2kHz(minERB)