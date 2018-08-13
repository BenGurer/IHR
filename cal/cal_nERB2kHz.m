function kHz = cal_nERB2kHz(nERB)
%% converts nERB (ERB numbers) to  kHz
A = 24.7/1000; B = 4.37;
kHz = 1/B*(exp(A*B*nERB)-1);
end