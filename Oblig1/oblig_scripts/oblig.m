%% FYS4565 - OBLIG1
clear all 
%clc

% Set working folder
cd ~/Desktop/FYS4565/Oblig1/oblig_scripts/

% Add to path getters and setters
addpath('get','set');

n = 2e4;
setParticleNumber(n);

% Run MAD-X
resetKickers;
energy = 1.0;
setEnergy(energy);
setEnergySpread(0.01);
setBpmMisalignments(0.0,0.0)
setQuadMisalignments(0,0);

runMADX;

%% Importing beam in workspace
% get initial beam 
Beam = getInitialBeam();
mass = 0.511e-3; %GeV

energy = getEnergy();
gamma = energy/mass;


%% calculating RMS emittance
% geometric
[Ex, Ey] = emittance(Beam);

fprintf('\n============================\n')
fprintf('Geometric emittance: \n')
fprintf('Ex = %d\nEy = %d\n',Ex,Ey)

% normalized

fprintf('----------------------------\n')
fprintf('Normalized emittance: \n')
fprintf('Ex_norm = %d\nEy_norm = %d\n',gamma*Ex,gamma*Ey)
fprintf('============================\n')


%% TWISS parameters
[Betax, Betay, ax, ay] = twiss(Beam);

fprintf('\n============================\n')
fprintf('TWISS parameters: \n')
fprintf('alpha_x = %d, beta_x=%d\nalpha_y = %d, beta_y=%d\n',ax,Betax,ay,Betay)

%% Plot beam orbit 
[x,y,s] = getBPMreadings;

fig1 = figure();
plot(s,x,'Linestyle','--','Marker','o','linewidth',1);
hold on;
grid on;
plot(s,y,'Linestyle','-.','Marker','x','linewidth',1);
axs = gca;
set(gca,'TickLabelInterpreter','latex','fontsize',11);
ylabel('[m]','Interpreter','Latex');
legend(axs,'$x$ orbit','$y$ orbit','Interpreter','Latex','fontsize',15,'Location','best');
title('BPM orbit','Fontsize',18,'Interpreter','Latex');
saveas(gcf,'../figures/orbit','jpg');


%% Introduce misalignments
setEnergy(energy)
setQuadMisalignments(0.001,0.001)
runMADX

% plot new orbit
[x,y,s] = getBPMreadings;

fig2 = figure();
plot(s,x,'Linestyle','--','Marker','o','linewidth',1);
hold on
grid on
plot(s,y,'Linestyle','-.','Marker','x','linewidth',1);
axs = gca;
set(gca,'TickLabelInterpreter','latex','fontsize',11);
ylabel('[m]','Interpreter','Latex');
legend(axs,'$x$ orbit','$y$ orbit','Interpreter','Latex','fontsize',15,'Location','best');
title('BPM orbit - misalignment','Fontsize',18,'Interpreter','Latex');
saveas(gcf,'../figures/orbit_mis','jpg');


%% Dispersion function
% Generate beam with different energies

energy = 1.0;
setEnergy(energy);
setQuadMisalignments(0.001,0.001);
setEnergyOffset(0.0);

runMADX;
[x0,y0,s0] = getBPMreadings;

newEnergy = 1.05;
newP = sqrt(newEnergy^2 + mass^2);
setEnergyOffset(newEnergy-energy);

runMADX;
[x1,y1,s] = getBPMreadings;

p0 = sqrt(energy^2 + mass^2);
dp = newP - p0;
Disp_x = (x1-x0)*p0/dp;
Disp_y = (y1-y0)*p0/dp;

% plot dispersion function
fig3 = figure();
plot(s,Disp_x,'Linestyle','--','Marker','none','linewidth',1.5)
hold on
plot(s,Disp_y,'Linestyle','-.','Marker','none','linewidth',1.5);
grid on
axs = gca;
set(gca,'TickLabelInterpreter','latex','fontsize',11);
ylabel('[m]','Interpreter','Latex');
legend(axs,'$D_x$','$D_y$','Interpreter','Latex','fontsize',15,'Location','best');
title('Dispersion function','Fontsize',18,'Interpreter','Latex');
saveas(gcf,'dispersion','jpg');


%% Emittance growth

n_points = 10;
setEnergy(1.0);
setEnergySpread(0);
setEnergyOffset(0.01);
misalign = linspace(0,2e-3,n_points);

[EmGrowthX, EmGrowthY] = EmittanceGrowth(misalign);

fig4 = figure();
plot(misalign,EmGrowthX,'Linestyle','--','Marker','o','linewidth',1)
hold on
plot(misalign,EmGrowthY,'Linestyle','-.','Marker','x','linewidth',1);
grid on
axs = gca;
set(gca,'TickLabelInterpreter','latex','fontsize',11);
legend(axs,'$\Delta\varepsilon_x/\varepsilon_x$','$\Delta\varepsilon_y/\varepsilon_y$','Interpreter','Latex','fontsize',15,'Location','best');
title('Emittance growth','Fontsize',18,'Interpreter','Latex');
saveas(gcf,'../figures/growth','jpg');

%% Emittance growth with energy spread

n_points = 10;
setEnergy(1.0);
setEnergySpread(0.01);
misalign = linspace(0,2e-3,n_points);

[EmGrowthX, EmGrowthY] = EmittanceGrowth(misalign);

fig4 = figure();
plot(misalign,EmGrowthX,'Linestyle','--','Marker','none','linewidth',1.5)
hold on
plot(misalign,EmGrowthY,'Linestyle','-.','Marker','none','linewidth',1.5);
grid on
axs = gca;
set(gca,'TickLabelInterpreter','latex','fontsize',11);
legend(axs,'$\Delta\varepsilon_x/\varepsilon_x$','$\Delta\varepsilon_y/\varepsilon_y$','Interpreter','Latex','fontsize',15,'Location','best');
title('Emittance growth - $\Delta E/E=0.01$\%','Fontsize',18,'Interpreter','Latex');
saveas(gcf,'../figures/growth_spread','jpg');


%% Change quadrupole strength
strength = linspace(1,1.8,5);

fig5 = figure();

for (i=1:length(strength))
    setQuadStrength(strength(i));

    [EmGrowthX, EmGrowthY] = EmittanceGrowth(misalign);

    plot(misalign,EmGrowthX,'Linestyle','--','Marker','x','linewidth',1.5)
    hold on
%plot(misalign,EmGrowthY,'Linestyle','-.','Marker','none','linewidth',1.5);
end

grid on
axs = gca;
set(gca,'TickLabelInterpreter','latex','fontsize',11);
legend(axs,'$k=$'+string(strength(1)),'$k=$'+string(strength(2)),'$k=$'+string(strength(3)),'$k=$'+string(strength(4)),'$k=$'+string(strength(5)),'Interpreter','Latex','fontsize',15,'Location','best');
title('Emittance growth - $\Delta E/E=0.01$\%','Fontsize',18,'Interpreter','Latex');
saveas(gcf,'../figures/growth_spread_k','jpg');


%% PART 2
% good cross check: compare \sigma = \beta\epsilon + D^2\sigma^2p/p0;

clear all 
clc

n = 2e4;
setParticleNumber(n)

energy = 1.0;
mass = 0.511e-3;

% running original orbit
setEnergy(1.0);
setEnergySpread(0.01);
setEnergyOffset(0.0);
setQuadMisalignments(0.001,0.001);
setBpmMisalignments(0,0);
resetKickers;

runMADX;
[x0,y0,s] = getBPMreadings;

% Adding the 1-1 steering corrections
x1 = OneToOne();

% plotting the orbits
fig1 = figure();
plot(s,x0,'Linestyle','--','Marker','o','linewidth',1.5);
hold on;
plot(s,x1,'Linestyle','-.','Marker','x','linewidth',1.5);
grid on;
axs = gca;
set(gca,'TickLabelInterpreter','latex','fontsize',11);
legend(axs,'Original','Corrected','Interpreter','Latex','fontsize',15,'Location','best');
title('BPM orbit - kickers','Fontsize',18,'Interpreter','Latex');

%% Plot emittance growth
setBpmMisalignments(0.0,0.0);
setEnergySpread(0.01);

n_points = 10;
misalign = linspace(0,2e-3,n_points);

resetKickers;
EmGrowthX0 = EmittanceGrowth(misalign,'none');
resetKickers;
EmGrowthX_c = EmittanceGrowth(misalign,'1to1');

fig2 = figure();
plot(misalign,EmGrowthX_c,'Linestyle','--','Marker','o','linewidth',1.5);
hold on;
plot(misalign,EmGrowthX0,'Linestyle','-.','Marker','x','linewidth',1.5);
grid on;
axs = gca;
set(gca,'TickLabelInterpreter','latex','fontsize',11);
legend(axs,'1-to-1 steer','Uncorrected','Interpreter','Latex','fontsize',15,'Location','best');
title('Emittance growth - 1to1','Fontsize',18,'Interpreter','Latex');


%% Introduce BPM misalignments
setBpmMisalignments(0.001,0.001);
setQuadMisalignments(0.001,0.001);
resetKickers;
runMADX;

% 1-to-1 steering

[x0,y0,s] = getBPMreadings;
x1 = OneToOne();

% plotting the orbits
fig3=figure();
plot(s,x0,'Linestyle','--','Marker','o','linewidth',1.5);
hold on;
plot(s,x1,'Linestyle','-.','Marker','x','linewidth',1.5);
grid on;
axs = gca;
set(gca,'TickLabelInterpreter','latex','fontsize',11);
legend(axs,'Original','Corrected','Interpreter','Latex','fontsize',15,'Location','best');
title('BPM orbit - BPM misalignment','Fontsize',18,'Interpreter','Latex');



%% Plot emittance growth
setBpmMisalignments(0.001,0.001);
setEnergySpread

n_points = 10;
misalign = linspace(0,2e-3,n_points);

resetKickers;
[EmGrowthX0, EmGrowthY0] = EmittanceGrowth(misalign,'none');
[EmGrowthX_c, EmGrowthY_c] = EmittanceGrowth(misalign,'1to1');

fig4 = figure();
plot(misalign,EmGrowthX_c,'Linestyle','--','Marker','none','linewidth',1.5);
hold on;
plot(misalign,EmGrowthX0,'Linestyle','-.','Marker','none','linewidth',1.5);
grid on;
axs = gca;
set(gca,'TickLabelInterpreter','latex','fontsize',11);
legend(axs,'1-to-1 steer','Uncorrected','Interpreter','Latex','fontsize',15,'Location','best');
title('Emittance growth - BPM misalignments','Fontsize',18,'Interpreter','Latex');



%% Dispersion-free steering

% Get trajectories and response matrix for reference beam

QuadMis = 1e-3;

% Get uncorrected trajectory
resetKickers;
setEnergy(1.0);
setQuadStrength(1.4142);
setEnergySpread(0.01);
setEnergyOffset(0.0);
setBpmMisalignments(0.001,0.001);
setQuadMisalignments(QuadMis,QuadMis);
runMADX;

[x0, y0, s] = getBPMreadings;

[x_c, y_c, s] = DFS(QuadMis);

% plot
plot(s,y0,'Linestyle','--','Marker','none','linewidth',1.5)
hold on
plot(s,y_c,'Linestyle','-.','Marker','none','linewidth',1.5)
grid on
axs = gca;
set(gca,'TickLabelInterpreter','latex','fontsize',11);
legend(axs,'Free','DFS','Interpreter','Latex','fontsize',15,'Location','best');
title('BPM orbit - DFS steering','Fontsize',18,'Interpreter','Latex');


%% Plot emittance growth
setEnergy(energy);
setBpmMisalignments(0.001,0.001);
setQuadMisalignments(0.001,0.001);

[EmGrowthX, EmGrowthY] = EmittanceGrowth_DFS(misalign)

fig5 = figure();
plot(misalign,EmGrowthX_c,'Linestyle','--','Marker','none','linewidth',1.5);
%hold on;
%plot(misalign,EmGrowthX0,'Linestyle','-.','Marker','none','linewidth',1.5);
grid on;
axs = gca;
set(gca,'TickLabelInterpreter','latex','fontsize',11);
legend(axs,'DFS steer','Uncorrected','Interpreter','Latex','fontsize',15,'Location','best');
title('Emittance growth - BPM misalignments','Fontsize',18,'Interpreter','Latex');


%% Functions

% RMS emittance
function [Ex , Ey] = emittance(Beam)

    Covx = cov(Beam(:,1),Beam(:,2));
    Covy = cov(Beam(:,3),Beam(:,4));
    
    Ex = sqrt(det(Covx));
    Ey = sqrt(det(Covy));
    
end

% Twiss parameters
function [Betax, Betay, ax, ay] = twiss(Beam)

    Covx = cov(Beam(:,1),Beam(:,2));
    Covy = cov(Beam(:,3),Beam(:,4));
    Ex = sqrt(det(Covx));
    Ey = sqrt(det(Covy));
    
    matrix = Covx/Ex;
    matriy = Covy/Ey;
    
    Betax = matrix(1,1);
    ax = -matrix(1,2);

    Betay = matriy(1,1);
    ay = -matriy(1,2);
    
end

% Emittance growth
function [EmGrowthX, EmGrowthY] = EmittanceGrowth(misalign,varargin)
    EmGrowthX = zeros(size(misalign));
    EmGrowthY = zeros(size(misalign));
    
    for (i=1:length(misalign))
        
        kick = char(varargin{1})
        
        switch kick
            case 'dfs'
                resetKickers;
                DFS(misalign(i));
        
            case '1to1'
                setQuadMisalignments(misalign(i),misalign(i));
                resetKickers;
                OneToOne();
                
            case 'none'
                setQuadMisalignments(misalign(i),misalign(i));
                resetKickers;
                runMADX;
        end     
   
        Beam0 = getInitialBeam;
        Beam1 = getFinalBeam;
   
        [ex0, ey0] = emittance(Beam0);
        [ex1, ey1] = emittance(Beam1);
   
        EmGrowthX(i) = (ex1-ex0)/ex0;
        EmGrowthY(i) = (ey1-ey0)/ey0; 
   
        clear ex0 ey0 ex1 ey1
    end
end


% 1-to-1 steering
function [x,y,s] = OneToOne()

    [x0,y0,s] = getBPMreadings;
    [Rx, Ry] = getResponseMatrix;
    setEnergy(1.0);
    setEnergyOffset(0);
    DthetaX = -pinv(Rx)*x0;
    DthetaY = -pinv(Ry)*y0;
    setKickers(DthetaX,DthetaY); 
    runMADX;

    [x,y] = getBPMreadings;
end


% DFS steering
function [x,y,s] = DFS(QuadMis)
    % get reference orbit 
    resetKickers;
    setEnergy(1.0);
    setEnergySpread(0.0);
    setEnergyOffset(0.0);
    setBpmMisalignments(0.001,0.001);
    setQuadMisalignments(QuadMis,QuadMis);
    runMADX;

    [x0, y0] = getBPMreadings;
    [Rx0, Ry0] = getResponseMatrix;
    setEnergy(1.0);
    setEnergyOffset(0);

    % Get trajectories and response matrix for test beam
    setEnergySpread(0.01);
    runMADX;

    [x1, y1] = getBPMreadings;
    [Rx1, Ry1] = getResponseMatrix;
    setEnergy(1.0);
    setEnergyOffset(0);

    % set up kickers
    DthetaY = -pinv(Ry1-Ry0)*(y1-y0);
    DthetaX = -pinv(Rx1-Rx0)*(x1-x0);

    % correct orbit
    setKickers(DthetaX,DthetaY);
    
    runMADX;

    [x, y, s] = getBPMreadings;
end



