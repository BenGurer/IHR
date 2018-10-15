function DLF= cal_DLF(f)

%% log DLF = a'sqrt(F) + k' + m'(SL^-1) % Nelson et al (1983)
% Nelson et al k' = -0.15, m' = 5.056, a' = 0.0214

%% Weir (1977) = log DLF = a sqrt(F) + b % depended on dB SL

a = 0.0214; k = -0.15; m = 5.056; SL = 40;
f = f.*1000; % convert from kHz to Hz

DLF = log10(a.*(sqrt(f)) + k + m.*(SL^-1));

end