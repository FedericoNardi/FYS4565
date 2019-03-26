% SETS NORMALIZED EMITTANCE PARAMETERS
function [] = setEmittances( emit_nx,  emit_ny )
    
    % set only (normalized) emittance params
    params = getParams();
    params(3) = emit_nx;
    params(4) = emit_ny;
    setParams(params);
    
end

