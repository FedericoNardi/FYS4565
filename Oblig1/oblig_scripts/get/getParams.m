% GET PARAMETERS FROM FILE "params.tfs"
% Format: [ nParticles, energy, emit_Nx, emit_Ny, sigma_E, quad_dx, quad_dy, bpm_dx, bpm_dy, kQuad, deltaE ]
function [ params ] = getParams()

    % declare parameter file
    paramsFile = 'params.tfs';
    
    % extract quadrupole misalignments in x and y
    headerLinesIn = 4;
    paramsData = importdata(paramsFile, ' ', headerLinesIn);
    
    % return parameter list
    params = paramsData.data;
    
end

