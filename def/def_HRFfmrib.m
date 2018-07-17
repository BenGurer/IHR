function y = get_HRFfmrib(x,xdata) 

t = xdata;

if length(x) == 5
    x(6) = 1;
elseif length(x) > 6
    error('oops! Your p is too long')
end

y = fmribHRF(t, x(1), x(2), x(3), x(4), x(5), x(6));

y = y ./ max(y);


end