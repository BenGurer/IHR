function ERB = cal_ERB(f)
%% calculate ERB at a given frequency in kHz
% ERBs as per Glasberg and Moore (1990);
A = 24.7/1000; B = 4.37;
ERB = A*(B*f+1);
end