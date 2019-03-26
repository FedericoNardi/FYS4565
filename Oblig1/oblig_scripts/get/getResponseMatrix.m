% CALCULATE RESPONSE MATRIX (BETWEEN KICKERS AND BPM READINGS)
function [Rx, Ry] = getResponseMatrix(delta)
    
    % Note: to save time in calculating response matrix, we use a
    %  single zero-emittance particle (reference orbit).
    
    % time-saver: remember beam params and optimise for tracking (reset later)
    particleNumber = getParticleNumber();
    sigma_E = getEnergySpread();
    [emit_Nx, emit_Ny] = getEmittances();
    setParticleNumber(1);
    setEnergySpread(0);
    setEmittances(0,0);
    
    % set beam offset energy (default 0 offset)
    if (nargin < 1); delta = 0; end
    setEnergyOffset(delta);
    
    % define number of kickers and bpms
    nBPMs = 8;
    nKickers = 8;
    
    % declare response matrices in x,y
    Rx = zeros(nBPMs, nKickers);
    Ry = zeros(nBPMs, nKickers);
    
    % track reference orbit
    resetKickers();
    fprintf('Tracking reference orbit... ');
    evalc('runMADX()');
    [xs0, ys0] = getBPMreadings();
    disp('Done');
    
    % amplitude of test kicks to be applied
    testKick = 1e-5;
    
    % do test kicks sequentially to calculate matrix in x,y
    fprintf('Calculating response from kickers... ');
    for j = 1:nKickers
        
        % set all kickers to zero but one
        xkicks = zeros(1,nKickers);
        ykicks = zeros(1,nKickers);
        xkicks(j) = testKick;
        ykicks(j) = testKick;
        setKickers(xkicks, ykicks);
        
        % track beam through
        evalc('runMADX()');
        fprintf([num2str(j) ' ']);
        
        % calculate response matrix row for active kicker
        [xs, ys] = getBPMreadings();
        Rx(:, j) = (xs - xs0)/testKick;
        Ry(:, j) = (ys - ys0)/testKick;
    
    end
    disp('Done.');
    
    % reset params to no energy offset
    setEnergyOffset(0);
    
    % time-saver: restore original beam
    setParticleNumber(particleNumber);
    setEnergySpread(sigma_E);
    setEmittances(emit_Nx, emit_Ny);
    
end

