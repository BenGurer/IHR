
function doublegammafun = get_HRFDoubleGamma(x,xdata) 
t = xdata;
amplitude = 1; %x(4);
offset = 0; %x(5);
% modelHrf = gampdf(t, params.x, 1) - gampdf(t, params.y, 1)/params.z;
modelHrf = gampdf(t, x(1), 1) - gampdf(t, x(2), 1)/x(3);
% compare normalising intergral with normalising max

% create wrapper function to allow normalisation

%normalise so that integral of sum = 1
modelHrf = modelHrf./sum(modelHrf(:));

% normalize the amplitude
if (max(modelHrf)-min(modelHrf))~=0
    modelHrf = (modelHrf-min(modelHrf)) ./ (max(modelHrf)-min(modelHrf));
end

doublegammafun= (amplitude*modelHrf+offset);
end