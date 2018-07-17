function [ gaussianFun ] = get_Gaussian(x,xdata)
mu = x(1);
sigma = x(2);

gaussianFun = exp(-(((xdata-mu).^2)/(2*(sigma^2))));

end