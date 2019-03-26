% BPM READINGS (MEAN BEAM POSITION) IN X, Y with optional BPM misalignments
function [ xs, ys, ss ] = getBPMreadings()
    
    % tracking file
    trackFile = 'particles.one';
    
    % find actual bpm offsets (not sigmas)
    [~, ~, bpm_dx, bpm_dy] = getBpmMisalignments();
    
    % number of BPMs
    nBPMs = 8;
    
    % declare x, y, s
    xs = zeros(nBPMs,1);
    ys = zeros(nBPMs,1);
    ss = zeros(nBPMs,1);
    
    % extract position data for each BPM
    fileID = fopen(trackFile);
    data = textscan(fileID, '%d %d %f %f %f %f %f %f %f %f', 'HeaderLines', 54, 'Delimiter', '\t', 'CommentStyle','#');
    fclose(fileID);
    particleNumber = getParticleNumber();
    for i = 1:nBPMs
        xs(i) = mean(data{3}((1:particleNumber) + i*particleNumber)) - bpm_dx(i); % mean x with BPM misalignment
        ys(i) = mean(data{5}((1:particleNumber) + i*particleNumber)) - bpm_dy(i); % mean y with BPM misalignment
        ss(i) = data{9}(i*particleNumber+1); % BPM s-position
    end
    
end

