% SETS ENERGY PARAMETER
function []= setEnergy( energy )
    
    % set only energy param
    params = getParams();
    params(2) = energy;
    setParams(params);
    
end
