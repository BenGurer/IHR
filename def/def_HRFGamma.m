
function gammafun = get_HRFGamma(x,xdata) 
time = xdata;
timelag = x(1);
tau = x(2);
exponent = round(x(3));
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

gammafun = (amplitude*gammafun+offset);

end