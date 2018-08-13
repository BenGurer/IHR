function nERB = cal_nERB(kHz)
%% converts kHz to nERB (ERB numbers)
A = 24.7/1000; B = 4.37;
nERB = 1/(A*B)*log(B*kHz+1);
end