% GET NORMALIZED EMITTANCES FROM PARAMETER FILE
function [ emit_Nx, emit_Ny ] = getEmittances()

    % reads energy from param file
    params = getParams();
    emit_Nx = params(3);
    emit_Ny = params(4);
    
end

