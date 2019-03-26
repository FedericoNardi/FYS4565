% SET ENERGY SPREAD
function [] = setEnergySpread(sigma_E)
    
    % set only energy offset param
    params = getParams();
    params(5) = sigma_E;
    setParams(params);

end
