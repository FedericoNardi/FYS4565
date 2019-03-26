% EXTRACT INITIAL BEAM (BEFORE LATTICE)
% Format: (x, xp, y, yp, z, dE/E)
%         [m, rad, m, rad, m, %/100]
function [ beam ] = getInitialBeam()
    
    % declare tracking file
    trackFile = 'particles.one';

    % extract beam 6D phase space (after last BPM)
    headerLinesIn = 54;
    track_data = importdata(trackFile, ' ', headerLinesIn);
    beam = track_data.data(:,3:8);
    
end
