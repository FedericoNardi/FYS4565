% GET QUADRUPOLE STRENGTH (k) FROM PARAMETER FILE
function [ kQuad ] = getQuadStrength()

    % reads quad strength from param file
    params = getParams();
    kQuad = params(10);
    
end


