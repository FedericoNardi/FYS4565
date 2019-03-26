function [] = plotBPMreadings()
    
    % get all bpm readings
    [xs, ys, ss] = getBPMreadings();
    
    sigma_E = getEnergySpread();
    [emit_Nx, emit_Ny] = getEmittances();
    [quad_dx, quad_dy] = getQuadMisalignments();
    [bpm_dx, bpm_dy] = getBpmMisalignments();
    
    % create plot with appropriate labels
    plot(ss, xs, ss, ys);
    title(['\sigma_E = ' num2str(sigma_E) ',  ' ...
           '\epsilon_N = (' num2str(emit_Nx) ', ' num2str(emit_Ny) ') m rad,  ' ...
           '\sigma_{QD} = (' num2str(quad_dx) ', ' num2str(quad_dy) ') m,  ' ...
           '\sigma_{BPM} = (' num2str(bpm_dx) ', ' num2str(bpm_dy) ') m']);
    xlabel('s / m');
    ylabel('BPM readings / m');
    legend('x','y');
    xlim([min(ss), max(ss)]);
    
end

