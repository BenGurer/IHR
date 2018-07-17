function rfModel = get_ROEX(x,xdata)
mu = x(1);
sigma = x(2);

fun = @(xdata,mu,sigma) 1 * exp(-(xdata - mu).^2/2/sigma^2);
pTW = integral(@(xdata)fun(xdata,mu,sigma),-100,100);

P = 4*mu/pTW;
g = abs(xdata-mu)/mu;
rfModel = (1+P*g).*exp(-P*g);

end