% GET PARTICLE NUMBER FROM PARAMETER FILE
function [ particleNumber ] = getParticleNumber()

    % reads energy from param file
    params = getParams();
    particleNumber = params(1);
    
end
