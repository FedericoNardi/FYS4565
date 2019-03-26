% GET PARTICLE NUMBER FROM PARAMETER FILE
function [ sigma_E ] = getEnergySpread()

    % reads energy spread from param file
    params = getParams();
    sigma_E = params(5);
    
end


