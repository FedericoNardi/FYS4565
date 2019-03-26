% RUNS SCRIPT "main.madx" IN MADX
function [] = runMADX()
    
    % main file
    mainFile = 'main.madx';
    addpath("/home/federico/Desktop/FYS4565/Oblig1/oblig_scripts/")
    % run MADX and print progress
    fprintf('Running MADX... ');
    syscmd = ['./madx ' mainFile];
    %evalc('system(syscmd)'); % quiet execution
    system(syscmd); % verbose execution
    disp('Done.');

end

