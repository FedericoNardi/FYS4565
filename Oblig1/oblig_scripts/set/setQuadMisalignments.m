% SET ST.DEV. OF NORMALLY DISTRIBUTED MISALIGNMENTS of QUADS
function [] = setQuadMisalignments(quad_dx, quad_dy)
    
    % set only energy offset param
    params = getParams();
    params(6) = quad_dx;
    params(7) = quad_dy;
    setParams(params);

end

