function nDLF = cal_nDLF(f)
%% define frequency resolution
% fs = 0:0.1:20;
% convert from kHz to Hz
f = f.*1000;
% constants for frequency discrimination
% SL is dB SL
a = 0.0214; k = -0.15; m = 5.056; % SL = 40;
SL = 45;
% y = exp(a .* sqrt(f) + k + m.*(SL.^-1));
% equation from Nelson
% DLM = exp(a .* sqrt(f) + k + m.*(SL.^-1));
% integral for 1/exp(a .* sqrt(f) + k + m.*(SL.^-1))
% describes how distance changes as a function of characteristic frequency
% Fc
nDLF = -((2.*(a.*sqrt(f) + 1) .* exp (-(a.* SL .* sqrt(f) + k .*SL + m)./SL))./a.^2);

% log df = exp(a*sqrt(f)+B)) Weir 1977

 % for 40 dB SL
% a = 0.026;
% b = -0.533;
% fd = exp(a*sqrt(fs)+b);
% df = 1./(exp(a*sqrt(fs)+b));
% integral of 1/e^(a*sqrt(f)+b)
% nDLF = -((2.*(a.*sqrt(f) +1).*exp(a.*(-sqrt(f)-b))./a.^2));

% integral of:
% 1/log(exp(a * sqrt(f) + k + m*(SL^-1))) = 
% nDLF = (2*(a.*sqrt(f) + (a.*sqrt(f) - log(exp(k + m/SL + a .* sqrt(f)))) .* log(SL .* log(exp((k + (m/SL) + a .* sqrt(f))))))) ./ a.^2;

end


% y = f;
% c = SL;
% 
% A = ((c * k + m)^2) +(log(a * c * sqrt(y) + c * k + m))/(a^2 * c^2);
% B = (sqrt(y) * (c*k +m)) / a*c;
% C = y * log(a * sqrt(y + (m/c) + k) - 0.5*y);
% 
% nDLF = 0.434294 * (-A + B + C);
% 
% 1/(a Sqrt[y] + m/c + k)
% 
% (2 a c Sqrt[y]-2 (c k+m) Log[c k+m+a c Sqrt[y]])/(a^2 c)

% (((c k+m) Sqrt[y])/(a c)-y/2+y Log[k+m/c+a Sqrt[y]]-((c k+m)^2 Log[c k+m+a c Sqrt[y]])/(a^2 c^2))/Log[10]
 
%  (2 (a Sqrt[y]+(a Sqrt[y]-Log[E^(k+m/c+a Sqrt[y])]) Log[c Log[E^(k+m/c+a Sqrt[y])]]))/a^2