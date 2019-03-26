% EXTRACT BPM MISALIGNMENT SIGMAS and ACTUAL NUMBERS
function [ sigma_dx, sigma_dy, bpm_dxs, bpm_dys ] = getBpmMisalignments()

    % misalignments from parameters (standard devs)
    params = getParams();
    sigma_dx = params(8);
    sigma_dy = params(9);
    
    % declare misalignment file
    misalignmentFile = 'bpmMisalignments.tfs';

    % extract quadrupole misalignments in x and y
    headerLinesIn = 8;
    misalignmentData = importdata(misalignmentFile, ' ', headerLinesIn);
    
    % read out misalignments of BPMs (MREX, MREY: "monitor readout error)
    bpm_dxs = misalignmentData.data(:, 43);
    bpm_dys = misalignmentData.data(:, 44);
    
end

