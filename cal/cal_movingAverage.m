function betas_mv = cal_movingAverage(betas)
nBins = 8;
windowAvSize = size(betas,1)/nBins;

if isreal(windowAvSize) && rem(windowAvSize,1)==0
    
loopLength = size(betas,1) - windowAvSize + 1;
betas_mv = zeros(loopLength,size(betas,2));

for i = 1:loopLength
betas_mv(i,:) = nanmean(betas(i:i+windowAvSize-1,:));
end

else
    error('Moving average window not an integer')
end

end