% EXTRACT QUAD MISALIGNMENT SIGMAS and ACTUAL NUMBERS
function [ sigma_dx, sigma_dy, quad_dxs, quad_dys ] = getQuadMisalignments()

    % misalignments from parameters (standard devs)
    params = getParams();
    sigma_dx = params(6);
    sigma_dy = params(7);
    
    % declare misalignment file
    misalignmentFile = 'quadMisalignments.tfs';

    % extract quadrupole misalignments in x and y
    headerLinesIn = 8;
    misalignmentData = importdata(misalignmentFile, ' ', headerLinesIn);
    
    % read out misalignments of quads (DX, DY)
    quad_dxs = misalignmentData.data(:, 43);
    quad_dys = misalignmentData.data(:, 44);
    
end

