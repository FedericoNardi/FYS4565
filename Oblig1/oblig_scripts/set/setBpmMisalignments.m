% SET ST.DEV. OF NORMALLY DISTRIBUTED MISALIGNMENTS of BPMS
function [] = setBpmMisalignments(bpm_dx, bpm_dy)
    
    % set only energy offset param
    params = getParams();
    params(8) = bpm_dx;
    params(9) = bpm_dy;
    setParams(params);

end

