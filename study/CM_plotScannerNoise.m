function av_dB = CM_plotScannerNoise
%% output: scanner noise at each stimulus frequency
% data = struct;
% scanID = 10:28;
% nScans = length(scanID);

% for i = 1:nScans
%     tf.i = i;
%     n =  scanID(i);

%% get study info
[stimInfo, glmInfo, pRFInfo, Info, plotInfo] = CM_setupStudyParams;

% use ispc to set data directory
if ispc
    Info.dataDir = 'E:\OneDrive - The University of Nottingham\data';
else
    Info.dataDir = '/Volumes/data_PSY/OneDrive - The University of Nottingham/data';
end

cd(fullfile(Info.dataDir, 'Pilot\11651_002 (Aud_Noise_Test)\AcousticMeasurements\' ));
filename = 'scan24FFT CH1.csv';
[tf.path,tf.file,tf.extension] = fileparts(filename);
%     filestring = [path '\' file extension];
transfer_function = lcfloadTransferFunction(tf);
%     tf.filename = filename;
% figure; semilogx(tf.frequencies,tf.fft);
plotSpectrum(transfer_function)

% convert frequency values from kHz to nERB
frequencies_nERB = cal_nERB(transfer_function.frequencies);
% interp to scale linearly spaced in nERBs
lin_nERB = linspace(cal_nERB(0.05),cal_nERB(16),1000);
dB_values = interp1(frequencies_nERB,transfer_function.fft,lin_nERB,'spline');
% average dB values at stimulus frequencies
% get stimulus values - could be an input..

stimulus = linspace(cal_nERB(0.1),cal_nERB(8),32);

[~, c_index] = min(abs(lin_nERB-stimulus(1)))

for i = 1:length(stimulus)
[~, u_index] = min(abs(lin_nERB-(stimulus(i)+0.5)))
[~, l_index] = min(abs(lin_nERB-(stimulus(i)-0.5)))

av_dB(i) = mean(dB_values(l_index:u_index))
end

figure; plot(stimulus,av_dB)
figure; plot(lin_nERB,dB_values)
    
    
end
function plotSpectrum(tf)
x = tf.frequencies;
y = tf.fft;
% Plot sound
% Create figure
f = figure ('Visible','off',...
    'Menu', 'None', ...
    'Color', [1 1 1]);

set(0,'Units','pixels')
p = get(0, 'MonitorPositions');
H = size (p);

switch H(1)
    case 1
        scnsize = p(1,:);
    case 2
        %                 scnsize = p(2,:);
        scnsize = p(1,:);
end

position = get(f,'Position');
outerpos = get(f,'OuterPosition');
borders = outerpos - position;
edge = -borders(1)/2;

pos =   [scnsize(1),...
    scnsize(4)./2,...
    (scnsize(3)./2) - edge,...
    scnsize(4)./2];


set(f,'OuterPosition',pos)

% Create axes
axes1 = axes('Parent',f,'YGrid','on',...
    'XScale','log',...
    'XMinorTick','on',...
    'XMinorGrid','on',...
    'XGrid','on');
xlim(axes1,[min(x) max(x)]);
ylim(axes1,[25 120]);
hold(axes1,'all');
plot(x,y)
title(tf.file)
xlabel('Frequency (kHz)');
ylabel('Amplitude (RMS) dB SPL');
text(max(xlim)-10,max(ylim)-10,sprintf('dBz RMS  = %.2f \ndBa RMS  = %.2f',tf.dB(1),tf.dB(2)))

% text(100,100,sprintf('dBz RMS  = %.2f \ndBa RMS  = %.2f',tf.dB(1),tf.dB(2)))
text(tf.spectrumPeak,tf.spectrumPeakValue,sprintf('Spectrum Peak = %.2f kHz\nPeak RMS  = %.2f dB SPL',tf.spectrumPeak,tf.spectrumPeakValue))
set (f, 'Visible', 'on');

% save plot
saveas(f,[tf.path '\' tf.file],'png')
end

function tf = lcfloadTransferFunction(tf)
% transfer functions measured using IHR Brüel&Kjær microphone and 2cc coupler

% [path,file,extension] = fileparts(filename);
ffts = csvread([tf.path '\' tf.file tf.extension],29,1);
%         ffts = csvread(filename,29,1);
        threshold = 100;
        index = ffts(:,end-2)<threshold;
        ffts(index,:) = nan;
        fftSum = ffts(:,end-2:end-1);
        tf.dB = nanmean(fftSum);
        fft = ffts(:,1:end-3); %remove last 3 columns corresponding to two different weighted averages and an empty column
        tf.frequencies = (0:size(fft,2)-1)*(20000./(size(fft,2)-1)); %assuming frequency resolution of 3.125 Hz, but should find a way to read from file (readcsv encounters an error)
        tf.frequencies = tf.frequencies/1000; %convert to kHz
        tf.fft = nanmean(fft);
        [peakValue, peakIndex] = max(tf.fft);
        tf.spectrumPeak = tf.frequencies(peakIndex);
        tf.spectrumPeakValue = peakValue;
        %     tf.fft = conv(tf.fft,ones(1,20)/20,'same');
end