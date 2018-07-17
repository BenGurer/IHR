function hrf = get_HRFBoxCar(x,t)

% x(1) = boxcar delay in seconds
% x(2) = boxcar duration in seconds
% t = time points of hrf

% delayS = 2.5;
% durationS = 2.5;
% sampleDuration = 1.5;
delayS = x(1);
durationS =x(2);
sampleDuration = t(2)-t(1);

totalDurationS = delayS+durationS;
totalDuration = round(totalDurationS/sampleDuration);  %total duration in samples
delay = round(delayS/sampleDuration);  %duration of delay in samples
duration = totalDuration - delay;
totalT = length(t);
hrf = [zeros(1,delay),ones(1,duration),zeros(1,totalT-totalDuration)];
% subplot(3,2,6)
% plot(t,hrfboxcar)