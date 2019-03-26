%% FYS4565 - OBLIG1
clear all 
%clc

% Set working folder
cd ~/Desktop/FYS4565/Oblig1/oblig_scripts/

% Add to path getters and setters
addpath('get','set');

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

fig1 = figure()
plot(s,x,'Linestyle','--','Marker','none','linewidth',2)
hold on
grid on
plot(s,y,'Linestyle','-.','Marker','none','linewidth',2)
axs = gca
set(gca,'TickLabelInterpreter','latex','fontsize',11)
legend(axs,'$x$ orbit','$y$ orbit','Interpreter','Latex','fontsize',15,'Location','best')
title('BPM orbit','Fontsize',18,'Interpreter','Latex')


%% Introduce misalignments
setEnergy(energy)
setQuadMisalignments(0.001,0.001)
runMADX

% plot new orbit
[x,y,s] = getBPMreadings;

fig2 = figure();
plot(s,x,'Linestyle','--','Marker','none','linewidth',2);
hold on
grid on
plot(s,y,'Linestyle','-.','Marker','none','linewidth',2);
axs = gca;
set(gca,'TickLabelInterpreter','latex','fontsize',11);
legend(axs,'$x$ orbit','$y$ orbit','Interpreter','Latex','fontsize',15,'Location','best');
title('BPM orbit - misalignment','Fontsize',18,'Interpreter','Latex');


%% Dispersion function
% Generate beam with different energies
energy = 1.0;
setEnergy(energy);
setQuadMisalignments(0.001,0.001);

runMADX;
[x0,y0,s0] = getBPMreadings;

newEnergy = 1.05;
newP = sqrt(newEnergy^2 + mass^2);
setEnergyOffset(newEnergy-energy);

runMADX;
[x1,y1,s1] = getBPMreadings;

p0 = sqrt(energy^2 + mass^2);
dp = newP - p0;
Disp_x = (x1-x0)*p0/dp;
Disp_y = (y1-y0)*p0/dp;

% plot dispersion function
fig3 = figure()
plot(s1,Disp_x,'Linestyle','--','Marker','none','linewidth',1.5)
hold on
plot(s1,Disp_y,'Linestyle','-.','Marker','none','linewidth',1.5);
grid on
axs = gca;
set(gca,'TickLabelInterpreter','latex','fontsize',11);
legend(axs,'$D_x$','$D_y$','Interpreter','Latex','fontsize',15,'Location','best');
title('Dispersion function','Fontsize',18,'Interpreter','Latex');


%% Emittance growth
% ASK: Very sensitive to  number of particles, it decreases for n=1e5, but
% not for n = 2e4
n = 2e4;
setParticleNumber(n);

n_points = 10
setEnergy(1.0)
setEnergySpread(0);
misalign = linspace(0,2e-3,n_points);
EmGrowthX = zeros(size(misalign));
EmGrowthY = zeros(size(misalign));

%Phasespaces = figure()

for (i=1:length(misalign))
   setQuadMisalignments(misalign(i),misalign(i))
   
   runMADX;
   
   Beam0 = getInitialBeam;
   Beam1 = getFinalBeam;
   
   [ex0, ey0] = emittance(Beam0);
   [ex1, ey1] = emittance(Beam1);
   
   EmGrowthX(i) = (ex1-ex0)/ex0;
   EmGrowthY(i) = (ey1-ey0)/ey0; 
   
   clear ex0 ey0 ex1 ey1
   
   
   %figure(Phasespaces)
   %clf;
   %subplot(2,2,1), plot(Beam0(:,1),Beam0(:,2),'r.')
   %subplot(2,2,2), plot(Beam0(:,3),Beam0(:,4),'b.')
   %subplot(2,2,3), plot(Beam1(:,1),Beam1(:,2),'r.')
   %subplot(2,2,4), plot(Beam1(:,3),Beam1(:,4),'b.')
   %hold on
   %pause(0.01)
   
end

fig4 = figure()
plot(misalign,EmGrowthX,'Linestyle','--','Marker','none','linewidth',1.5)
hold on
plot(misalign,EmGrowthY,'Linestyle','-.','Marker','none','linewidth',1.5);
grid on
axs = gca;
set(gca,'TickLabelInterpreter','latex','fontsize',11);
legend(axs,'$D_x$','$D_y$','Interpreter','Latex','fontsize',15,'Location','best');
title('Emittance growth - n= '+string(n),'Fontsize',18,'Interpreter','Latex');





%% PART 2
%% Correct orbit with kickers 
clear all 
clc

% running original orbit
setEnergy(1.0);
setQuadMisalignments(0.001,0.001);
setBpmMisalignments(0,0);
resetKickers;

runMADX;
[x0,y0,s0] = getBPMreadings;

% Adding the 1-1 steering corrections
[Rx, Ry] = getResponseMatrix;
DthetaX = -pinv(Rx)*x0;
DthetaY = -pinv(Ry)*y0;
setKickers(DthetaX,DthetaY);
fprintf('\n-----kickers are set!-----\n')

runMADX;
[x1,y1,s1] = getBPMreadings;

% plotting the orbits
plot(s0,x0,'Linestyle','--','Marker','none','linewidth',1.5)
hold on
plot(s1,x1,'Linestyle','-.','Marker','none','linewidth',1.5)
grid on
axs = gca;
set(gca,'TickLabelInterpreter','latex','fontsize',11);
legend(axs,'Original','Corrected','Interpreter','Latex','fontsize',15,'Location','best');
title('BPM orbit - kickers','Fontsize',18,'Interpreter','Latex');

%% Plot emittance growth




%% Introduce BPM misalignments
setBpmMisalignments(0.001,0.001);
runMADX;

[x1,y1,s1] = getBPMreadings;

% plotting the orbits
plot(s0,x0,'Linestyle','--','Marker','none','linewidth',1.5)
hold on
plot(s1,x1,'Linestyle','-.','Marker','none','linewidth',1.5)
grid on
axs = gca;
set(gca,'TickLabelInterpreter','latex','fontsize',11);
legend(axs,'Original','Corrected','Interpreter','Latex','fontsize',15,'Location','best');
title('BPM orbit - BPM misalignment','Fontsize',18,'Interpreter','Latex');

%% Dispersion-free steering
% ASK! 

clear all 
clc
% Get trajectories and response matrix for reference beam
resetKickers;
energy = 1.0;
setEnergy(energy);
setEnergySpread(0.01);
setBpmMisalignments(0.0,0.0)
setQuadMisalignments(0,0);
runMADX;

[x0, y0, s] = getBPMreadings;
[Rx0, R0] = getResponseMatrix;

% Get trajectories and response matrix for test beam
setBpmMisalignments(0.001,0.001)
setQuadMisalignments(0.001,0.001)
runMADX;

[x1, y1, s] = getBPMreadings;
[Rx1, R1] = getResponseMatrix;

% set up kickers
DthetaY = -pinv(R1-R0)*(y1-y0);
DthetaX = -pinv(Rx1-Rx0)*(x1-x0);

% correct orbit
setKickers(DthetaX,DthetaY);
setEnergy(1.0);
setEmittances(1e-5,1e-5);
setQuadStrength(1.4142);
runMADX;

[x_c, y_c, s] = getBPMreadings;

% plot
plot(s,y1,'Linestyle','--','Marker','none','linewidth',1.5)
hold on
plot(s,y_c,'Linestyle','-.','Marker','none','linewidth',1.5)
grid on
axs = gca;
set(gca,'TickLabelInterpreter','latex','fontsize',11);
legend(axs,'Free','DFS','Interpreter','Latex','fontsize',15,'Location','best');
title('BPM orbit - DFS steering','Fontsize',18,'Interpreter','Latex');



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





