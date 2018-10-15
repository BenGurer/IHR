function d = binFreq_pilot_hrf(d)

% keyboard
groupSize = 20;
loopLength = (length(d.stimNames))/groupSize;
c = 1;
% index = {'a' 'b' 'c' 'd' 'e' 'f' 'g' 'h' 'i' 'j' 'k' 'l' 'm' 'n' 'o' 'p' 'r' 's' 't' 'u' 'v' 'w' 'x' 'y' 'z'};
for i = 1:loopLength
% %     c = c + 1;
% %     for ii = 1:groupSize
% %     d.stimfile{j}.stimNames{c+ii} = i;
% %     end
stimNames{i} = [d.stimNames{c} ' to ' d.stimNames{c+groupSize-1}];
% stimNames{i} = num2str(i);
% stimNames{i} = index(i);
stimdurations{i} = [d.stimDurations{c:c+groupSize-1}];
stimvol{i} = [d.stimvol{c:c+groupSize-1}];
    c = c + groupSize;
end

d.stimNames = stimNames;
% d.stimNames = index(1:loopLength);
d.stimDurations = stimdurations;
d.stimvol = stimvol;

end