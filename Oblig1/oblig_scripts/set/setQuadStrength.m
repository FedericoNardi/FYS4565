% SET QUADRUPOLE STRENGTH (k)
function [] = setQuadStrength(kQuad)
    
    % set only quadrupole strength param
    params = getParams();
    params(10) = kQuad;
    setParams(params);

end
