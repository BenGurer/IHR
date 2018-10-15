function f = cal_nDLF2kHz(nDLF)
%% define frequency resolution
fs = 0:0.01:20;
% convert from kHz to Hz
fs = fs.*1000;
% constants for frequency discrimination
% SL is dB SL
a = 0.0214; k = -0.15; m = 5.056; SL = 40;
% y = exp(a .* sqrt(f) + k + m.*(SL.^-1));
% equation from Nelson
% DLM = exp(a .* sqrt(f) + k + m.*(SL.^-1));
% integral for 1/exp(a .* sqrt(f) + k + m.*(SL.^-1))
% describes how distance changes as a function of characteristic frequency
% Fc
d = -((2.*(a.*sqrt(fs) + 1) .* exp (-(a.* SL .* sqrt(fs) + k .*SL + m)./SL))./a.^2);

% log df = exp(a*sqrt(f)+B)) Weir 1977

 % for 40 dB SL
% a = 0.026;
% b = -0.533;
% fd = exp(a*sqrt(fs)+b);
% df = 1./(exp(a*sqrt(fs)+b));
% integral of 1/e^(a*sqrt(f)+b)
% d = -((2.*(a.*sqrt(fs) +1).*exp(a.*(-sqrt(fs)-b))./a.^2));
dq = 0:0.001:30; % distance of cortex
% x = cortical distance = d
% v = characterisic frequency = f

f = interp1(d,fs,nDLF,'spline');
f = f./1000

((c k + m)^2 (1 + 
   ProductLog[-((E^(-1 - m/(c k) - (a^2 y)/(2 k)))^(((c k)/(
     c k + m)))/(c k + m))])^2)/(a^2 c^2)

end