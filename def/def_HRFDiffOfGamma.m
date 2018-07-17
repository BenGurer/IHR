function gammafun = get_HRFDiffOfGamma(x,xdata) 

time = xdata;
timelag = x(1);
tau = x(2);
exponent = round(x(3));
timelag2 = x(4);
tau2 = x(5);
exponent2 = round(x(6));
amplitude2 = x(7);
amplitude = 1; % x(4);
offset = 0; % x(5);

% exponent = round(exponent);
% gamma function
gammafun = (((time-timelag)/tau).^(exponent-1).*exp(-(time-timelag)/tau))./(tau*factorial(exponent-1));

% negative values of time are set to zero,
% so that the function always starts at zero
gammafun(find((time-timelag) < 0)) = 0;

% normalize the amplitude
% if (max(gammafun)-min(gammafun))~=0
%     gammafun = (gammafun-min(gammafun)) ./ (max(gammafun)-min(gammafun));
% end

gammafun = gammafun ./max(gammafun);
gammafun2 = amplitude2 .* (((time-timelag2)/tau2).^(exponent2-1).*exp(-(time-timelag2)/tau2))./(tau2*factorial(exponent2-1));

gammafun = gammafun - gammafun2;
gammafun = (amplitude*gammafun+offset);

end