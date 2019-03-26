% SET PARAMETERS USED BY MADX
% Format: [ nParticles, energy, emit_Nx, emit_Ny, sigma_E, quad_dx, quad_dy, bpm_dx, bpm_dy, kQuad, deltaE ]
function [] = setParams(params)
    
    % define parameter file used by MADX
    paramFile = 'params.tfs';
        
    % save parameters
    fileID = fopen(paramFile,'w');
    fprintf(fileID, '@ NAME %%06s "PARAMS"\n');
    fprintf(fileID, '@ TYPE %%04s "USER"\n');
    fprintf(fileID, '* N_PARTICLES ENERGY EMIT_NX EMIT_NY SIGMA_E QUAD_DX QUAD_DY BPM_DX BPM_DY K_QUAD DELTA_E\n');
    fprintf(fileID, '$ %%le %%le %%le %%le %%le %%le %%le %%le %%le %%le %%le\n');
    fprintf(fileID, ' %4.2e %4.2e %4.2e %4.2e %4.2e %4.2e %4.2e %4.2e %4.2e %4.2e %4.2e\n', params);
    fclose(fileID);

end
