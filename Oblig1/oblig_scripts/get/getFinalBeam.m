% EXTRACT FINAL BEAM (AFTER LATTICE)
% Format: (x, xp, y, yp, z, dE/E)
function [ beam ] = getFinalBeam()
    
    % define tracking file
    trackFile = 'particles.one';

    % extract beam 6D phase space (after last BPM)
    bpmNumber = 8;
    headerLinesIn = 54 + bpmNumber*(getParticleNumber()+1);
    track_data = importdata(trackFile, ' ', headerLinesIn);
    beam = track_data.data(:,3:8);
    
end

    
