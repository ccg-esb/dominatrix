
clc
clear all
close all

run('lib/addpath_recurse');
addpath_recurse('lib/');
addpath_recurse('src/');
 %%
 

params=defineParams();

params.numDivisions=10;
params.T=1;
numLevels=6;
nH=10;
Hs=linspace(1/nH,1,nH);  %Dominances to test


G1=10;
R12=1;

params.iniPlasmids=[G1, R12];

%% PLOT TREE

[freqsT, this_frac_survivors, num_mutations_tree]=plotTree(params, numLevels, Hs);
 


%% PLOT MANY REALIZATIONS


nsims=10;

cmap=cbrewer('qual', 'Paired', 6);

figure('Position', [500 500 1000 200])
clf('reset');set(gcf,'DefaultLineLineWidth',1); set(gcf, 'color', 'white');
for i=1:nsims
    [itime, ipG1, ipR12]=simulatePlasmidDynamics(params);

    stairs(itime,ipG1,'-','Color',cmap(1,:),'LineWidth',1); hold on;
    stairs(itime,ipR12,'-','Color',cmap(5,:),'LineWidth',1);
    axis([0 itime(end) -.5 params.pcn+.5]);
end
set(gca,'FontSize',18)
xlabel('Time (divisions)','FontSize',20);
ylabel('PCN','FontSize',20);

legend('WT','MUT','Location','NorthEastOutside');

export_fig 'figures/_model_manyRealizations.pdf';

%% PLOT ONE REALIZATION
%close all
[t, A]=simulateReplication(params, [G1; R12], 1);
[time, pG1, pR12]=simulatePlasmidDynamics(params);

figure('Position', [500 500 1000 200])
clf('reset');
set(gcf,'DefaultLineLineWidth',1); set(gcf, 'color', 'white');

    stairs(time,pG1,'-','Color',cmap(2,:),'LineWidth',2); hold on;
    stairs(time,pR12,'-','Color',cmap(6,:),'LineWidth',2);
axis([0 time(end) -.5 params.pcn+.5]);

set(gca,'FontSize',18)
xlabel('Time (divisions)','FontSize',20);
ylabel('PCN','FontSize',20);

legend('WT','MUT','Location','NorthEastOutside');

export_fig 'figures/_model_oneRealization.pdf';

