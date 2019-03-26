% GET ENERGY OFFSET (delta) FROM PARAMETER FILE
function [ delta ] = getEnergyOffset()

    % reads energy offset from param file
    params = getParams();
    delta = params(11);
    
end

