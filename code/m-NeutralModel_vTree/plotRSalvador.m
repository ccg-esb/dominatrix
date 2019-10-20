

clc
clear all
close all

run('lib/addpath_recurse');
addpath_recurse('lib/');
addpath_recurse('src/');


%% PARAMETERS

mut_rate=1e-8;
maxLevels=25;
PCNs=[1,2,4,10,20,40,100];  %plasmid copy numbers
nH=20;

figurePath='../../figures/';
dataPath='../../data/runs/';
expePath=['neutralModel_mutRate',num2str(mut_rate),'_Levels',num2str(maxLevels),'/'];

toFile=1;
cmap=[cbrewer('seq', 'OrRd', length(PCNs))];
numCells=2^maxLevels;
dirName=[dataPath,'',expePath];

%% LOAD DATA

% Experimental data
data_strains={'pBAD::{\it gyrA}', 'pBAD::{\it folA}'};
data_Hs=[0.0086, .874];
data_mutRate=[1.38369504317873E-10, 1.65977010518979E-09];
data_95Low=[6.23933410135799E-11, 8.70183184689011E-10];
data_95Upper=[2.55092181419044E-10, 1.21921602655543e-08];
data_normMutRate=[2.00190193286432e-09,5.47345788943695E-10]; 

% MaximumLikelihood DATA PRODUCED BY rSalvador
datamut=csvread([dirName,'MLmut.csv']); 

%% PLOT H vs mut (lines)

figure(1)
clf('reset');
set(gcf,'DefaultLineLineWidth',2); set(gcf, 'color', 'white');

leg={};
filter=unique([1:2:nH nH]);
indx_norm=datamut(:,1)==1;
normMut=mean(datamut(indx_norm,3)./numCells);
for iPCN=2:length(PCNs)
    
    PCN=PCNs(iPCN);
    indx=find(datamut(:,1)==PCN);
    leg{iPCN-1}=[num2str(PCN),' copies'];
    
    Hs=datamut(indx,2);
    Ms=datamut(indx,3);
    freqMuts=Ms/numCells;
    
    semilogy(Hs(filter), freqMuts(filter)/normMut,'-','Color',cmap(iPCN,:),'LineWidth',2 ); hold all;
end

set(gca, 'YScale', 'log')
semilogy(Hs, ones(1,length(Hs)),'-','Color','k'); hold all;
set(gca,'FontSize',18)
xlabel('Dominance (h)','FontSize',22);
ylabel('-fold change in mutation rate (log)','FontSize',22);

ylim([.002, 4e2]);
xlim([0, 1.05]);

legappend('')
legend(leg,'FontSize',18,'Location','NorthWest')
if toFile
    eval(['export_fig ',figurePath,'H_vs_FreqMut_mutRate',num2str(mut_rate),'_PCNlines.pdf']);
end

%% PLOT H vs mut (inset)

figure('Position', [500 500 240 200])
clf('reset');
set(gcf,'DefaultLineLineWidth',2); set(gcf, 'color', 'white');

leg={};
filter=unique([1:2:nH nH]);
indx_norm=find(datamut(:,1)==1);
normMut=mean(datamut(indx_norm,3)./numCells);

PCN=20;
indx=find(datamut(:,1)==PCN);

Hs=datamut(indx,2);
Ms=datamut(indx,3);
MsLower=datamut(indx,4);
MsUpper=datamut(indx,5);

freqMuts=Ms/numCells;
freqMutsLower=MsLower/numCells;
freqMutsUpper=MsUpper/numCells;

semilogy(Hs(filter), freqMutsLower(filter)/normMut,'-','Color',cmap(iPCN-2,:),'LineWidth',1 ); hold all;
semilogy(Hs(filter), freqMutsUpper(filter)/normMut,'-','Color',cmap(iPCN-2,:) ,'LineWidth',1); hold all;
p_model=semilogy(Hs(filter), freqMuts(filter)/normMut,'-','Color',cmap(iPCN-1,:),'LineWidth',2 ); hold all;

set(gca, 'YScale', 'log')
semilogy([0 2], [1 1],':','Color','k','LineWidth',1); hold all;

set(gca,'FontSize',18)
xlabel('h','FontSize',22);
ylabel('','FontSize',22);
ylim([.02, .5e2]);
xlim([-.05, 1.05]);

%Plot data
for i=1:length(data_Hs)
    errorbar(data_Hs(i), data_mutRate(i)./data_normMutRate(i), data_95Low(i)./data_normMutRate(i), data_95Upper(i)./data_normMutRate(i), 'k.','LineWidth',2);
    plot(data_Hs(i), data_mutRate(i)./data_normMutRate(i),'ko','MarkerFaceColor','k','LineWidth',1,'MarkerSize',10); hold on;
    
    if i<length(data_Hs)
        d=0.05;
        dy=0;
    else
        d=-0.6;
        dy=5;
    end
    text(data_Hs(i)+d, data_mutRate(i)./data_normMutRate(i)+dy, data_strains{i},'FontSize',18)
end


if toFile
    eval(['export_fig ',figurePath,'H_vs_FreqMut_mutRate',num2str(mut_rate),'_PCNdata.pdf']);
end

