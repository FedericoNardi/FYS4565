% SETS ENERGY OFFSET (delta) PARAMETER
function [] = setEnergyOffset( delta )
    
    % set only energy offset param
    params = getParams();
    params(11) = delta;
    setParams(params);
    
end

