% GET ENERGY FROM PARAMETER FILE
function [ energy ] = getEnergy()

    % reads energy from param file
    params = getParams();
    energy = params(2);
    
end

