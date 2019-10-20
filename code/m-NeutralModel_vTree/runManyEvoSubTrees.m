
%clc
clear all
close all

run('lib/addpath_recurse');
addpath_recurse('lib/');
addpath_recurse('src/');
 

%% PARAMETERS

dataPath='/Users/ESB/ESB_DATA/DOMINATRIX_data/runs/';

maxLevels=25;  %Depth of tree to simulate
PCNs=[1,2,4,10,20,30,40,100,200];  %plasmid copy numbers
mut_rate=1e-8;  %mutation rate
N=100;  %number of simulations
nH=200;

%%
Hs=linspace(1/nH,1,nH);  %Dominances to test
toFile=1;  %Save data


%% RUN SIMULATION

%params=defineParams();

if toFile
    %expePath=['_',datestr(now,'yyyymmdd_hhMMSS'),'_neutralModel_mutRate',num2str(mut_rate),'_Levels',num2str(maxLevels),'_N',num2str(N),'/'];
    expePath=['neutralModel_mutRate',num2str(mut_rate),'_Levels',num2str(maxLevels),'/'];
    
    dirName=[dataPath,'',expePath];
    if ~exist(dirName,'dir')
        mkdir(dirName);
        mkdir([dirName,'data/']);
        mkdir([dirName,'figures/']) 
    end
end

mean_survivors=zeros(length(PCNs), nH);
median_survivors=zeros(length(PCNs), nH);
std_survivors=zeros(length(PCNs), nH);
total_mutations=zeros(1,length(PCNs));

maxCells=2.^maxLevels;
parfor ipcn=1:length(PCNs)
    disp([newline,' PCN=',num2str(PCNs(ipcn))]);
    
    %params=defineParams();
    params=[];
    params.rho=.05; %Plasmid replication rate
    params.pcn=PCNs(ipcn); %Mean plasmid copy number
    params.iniPlasmids = [params.pcn; 0];  %Initial plasmid copy number
    params.mut_rate=mut_rate; %Mutation rate
    
    this_survivors=zeros(N,nH);
    frac_survivors=zeros(N,nH);
    for n=1:N  %replicate simulations
        fprintf('%s', '.'); 
        if mod(n,100)==0
            fprintf('%s',num2str(n));
        end

        %Find mutations
        num_mutations_level=zeros(1,maxLevels);
        for level=1:maxLevels
            num_cells_level=2.^level;
            rs=rand(1,num_cells_level);
            [ipcnt, ~]=find(rs>(1-params.mut_rate.*params.pcn));
            num_mutations_level(level)=length(ipcnt);
            %disp(['Number of mutations (level ',num2str(level),')= ',num2str(num_mutations_level(level)),'/',num2str(num_cells_level)]);
        end

        %Simulate mutated lineages
        numCellsMutatedLineages=0;
        num_mutations_tree=[];
        for M=1:maxLevels
            Ms=num_mutations_level(M);
            for mi=1:Ms
                
                %fprintf('%s', ['(',num2str(M),')']);
                numLevels=maxLevels-M;
                numCellsLevel=2.^numLevels;
                [~, this_frac_survivors, num_mutations_tree]=simTree(params, numLevels, Hs);
                
                if ~isempty(this_frac_survivors)  %If not time-out
                
                    lineage_survivors=numCellsLevel.*this_frac_survivors;
                    this_survivors(n,:)=this_survivors(n,:)+lineage_survivors;

                    numCellsMutatedLineages=numCellsMutatedLineages+numCellsLevel;
                    %disp(['    mut@level:',num2str(M),'   lineage survivors: [',num2str(lineage_survivors),'] CellsInThisLineage:',num2str(numCellsLevel),'']);
               
                else  %Here we should break tree and simulate it
                    this_survivors(n,:)=0;  %TMP
                
                end
            end
        end
        numCellsOther=2.^(maxLevels)-numCellsMutatedLineages;
        
        frac_survivors(n,:)=this_survivors(n,:)./(numCellsMutatedLineages+numCellsOther);
        %disp([num2str(n),'> Total survivors: [',num2str(this_survivors(n,:)),'], NumberOfMutations: ',num2str(sum(num_mutations_level)),'+',num2str(sum(num_mutations_tree)),', TotalCellsInMutantLineages:',num2str(numCellsMutatedLineages),', NonMutantCells:',num2str(numCellsOther)]);

        if toFile
             for kk=1:nH
                 
                 fileName=['sim_PCN',num2str(PCNs(ipcn)),'_H',num2str(Hs(kk)*100),'e-2.txt'];
                 fileID = fopen([dirName,'data/',fileName], 'a');
                 
                 fprintf(fileID, '%.0f\n', this_survivors(n,kk));
                 fclose(fileID); 
             end
           
             %Export mutations
             fileNameMut=['mutations_PCN',num2str(PCNs(ipcn)),'.txt'];
             fileIDmut = fopen([dirName,'data/',fileNameMut], 'a'); 
             fprintf(fileIDmut, '%.0f\t', sum(num_mutations_level));
             for expLevel=1:maxLevels
                fprintf(fileIDmut, '%.0f\t', num_mutations_level(expLevel));
             end
             fprintf(fileIDmut,'\n');
             fclose(fileIDmut); 
             
        end
        
    end
    mean_survivors(ipcn,:)=nanmean(frac_survivors);
    median_survivors(ipcn,:)=nanmedian(frac_survivors);
    std_survivors(ipcn,:)=std(frac_survivors);
end

%%


