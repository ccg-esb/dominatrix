function [freqsT, frac_survivors, mutations_level]=plotTree(params, numLevels, Hs)


if nargin==0
    disp('please run me from runManyTrees');
    return;
end
toPlot=1;

cmap=[cbrewer('seq', 'OrRd', 201)];
nH=length(Hs);

%%
tG1_0 = tree(1);
tR12_0 = tree(1);

tFreq_T = tree(1);

%params.iniPlasmids = [params.pcn-1; 1];  %Initial plasmid copy number: invasion experiment
params.iniPlasmids = [params.pcn; 0];  %Initial plasmid copy number: fluctuation experiment

prev_ns=[1];  %parental node
freqsT=[];

mutations_level=zeros(1,numLevels+1);
for level=1:numLevels+1
    
    freqsT=[];
    ns=[];
    for i=1:prev_ns
        i1=length(ns)+1;
        i2=length(ns)+2;
        
        if level>1
            prevG1=tG1_0.get(prev_ns(i));
            prevR12=tR12_0.get(prev_ns(i));
            
            %Simulate Mutation in lineage
            if rand()>(1-prevG1*params.mut_rate)  %Mutation G1->R12
                prevG1=prevG1-1;
                prevR12=prevR12+1;
                fprintf('%s', ['[',num2str(level),']']);
                mutations_level(level)=mutations_level(level)+1;
            end
        else
            prevG1=params.iniPlasmids(1)-1; %;
            prevR12=1; %params.iniPlasmids(2);
        end
        
        %Replication Dynamics (randomly replicate plasmids)
        [~, rep_y]=simulateReplication(params, [prevG1; prevR12], 1);
        p1=rep_y(end,1);
        p2=rep_y(end,2);
        
        freq=p2/(p1+p2);  %fraction of R12
        
        %Update terminal condition of previous node
        tFreq_T = tFreq_T.set(prev_ns(i), freq);
        
        %Division event (randomly segregate plasmids)
        if level<=numLevels
            [ tFreq_T, ns(i1) ] = tFreq_T.addnode(prev_ns(i), -1);
            [ tFreq_T, ns(i2) ] = tFreq_T.addnode(prev_ns(i), -1);
            
            if params.pcn>1
                if p1>0
                    pG1_0=binornd(floor(p1),0.5,1);
                else
                    pG1_0=0;
                end
                
                if p2>0
                    pR12_0=binornd(floor(p2),0.5,1);
                else
                    pR12_0=0;
                end
                [ tG1_0, ns(i1) ] = tG1_0.addnode(prev_ns(i), pG1_0);
                [ tG1_0, ns(i2) ] = tG1_0.addnode(prev_ns(i), p1-pG1_0);
                [ tR12_0, ns(i1) ] = tR12_0.addnode(prev_ns(i), pR12_0);
                [ tR12_0, ns(i2) ] = tR12_0.addnode(prev_ns(i), floor(p2)-pR12_0);
                
            else
                
                [ tG1_0, ns(i1) ] = tG1_0.addnode(prev_ns(i), p1);
                [ tG1_0, ns(i2) ] = tG1_0.addnode(prev_ns(i), p1);
                [ tR12_0, ns(i1) ] = tR12_0.addnode(prev_ns(i), p2);
                [ tR12_0, ns(i2) ] = tR12_0.addnode(prev_ns(i), p2);
                
            end
        end
        freqsT(i)=freq;
    end
    prev_ns=ns;
end

%COMPUTE FRACTION OF SURVIVORS FOR EACH H
frac_survivors=zeros(1,nH);
for iy=1:nH
    survivors=0;
    for ix=1:length(freqsT)
        if freqsT(ix)>(1-Hs(iy))
            survivors=survivors+1;
        end
    end
    %survivors
    frac_survivors(iy)=survivors/length(freqsT);
end

disp(tFreq_T.tostring)
%% Plot TREE

figure('Position', [500 500 1000 1000])
clf('reset');set(gcf,'DefaultLineLineWidth',1); set(gcf, 'color', 'white');
[vlh, hlh, tlh] = tFreq_T.plot([], 'TextRotation', 90);

colormap(cmap);

aboveTreshold1 = tFreq_T >= .1; % H>0.1
aboveTreshold2 = tFreq_T >= .2; % H>0.2
aboveTreshold3 = tFreq_T >= .3; % H>0.3
aboveTreshold4 = tFreq_T >= .4; % H>0.4
aboveTreshold5 = tFreq_T >= .5; % H>0.5
aboveTreshold6 = tFreq_T >= .6; % H>0.6
aboveTreshold7 = tFreq_T >= .7; % H>0.7
aboveTreshold8 = tFreq_T >= .8; % H>0.8
aboveTreshold9 = tFreq_T >= .9; % H>0.9
aboveTreshold10 = tFreq_T >= 1; % H=1


iterator = aboveTreshold10.depthfirstiterator;
for i = iterator
    %
    set( tlh.get(i), 'Color' , [1 1 1] ); %hide labels
    
    if  aboveTreshold10.get(i)
        set( vlh.get(i), 'Color' , cmap(10,:) )
        set( hlh.get(i), 'Color' , cmap(10,:) )
    elseif  aboveTreshold9.get(i)
        set( vlh.get(i), 'Color' , cmap(9,:) )
        set( hlh.get(i), 'Color' , cmap(9,:) )
    elseif  aboveTreshold8.get(i)
        set( vlh.get(i), 'Color' , cmap(8,:) )
        set( hlh.get(i), 'Color' , cmap(8,:) )
    elseif  aboveTreshold7.get(i)
        set( vlh.get(i), 'Color' , cmap(7,:) )
        set( hlh.get(i), 'Color' , cmap(7,:) )
    elseif  aboveTreshold6.get(i)
        set( vlh.get(i), 'Color' , cmap(6,:) )
        set( hlh.get(i), 'Color' , cmap(6,:) )
    elseif  aboveTreshold5.get(i)
        set( vlh.get(i), 'Color' , cmap(5,:) )
        set( hlh.get(i), 'Color' , cmap(5,:) )
    elseif  aboveTreshold4.get(i)
        set( vlh.get(i), 'Color' , cmap(4,:) )
        set( hlh.get(i), 'Color' , cmap(4,:) )
    elseif  aboveTreshold3.get(i)
        set( vlh.get(i), 'Color' , cmap(3,:) )
        set( hlh.get(i), 'Color' , cmap(3,:) )
    elseif  aboveTreshold2.get(i)
        set( vlh.get(i), 'Color' , cmap(2,:) )
        set( hlh.get(i), 'Color' , cmap(2,:) )
    elseif  aboveTreshold1.get(i)
        set( vlh.get(i), 'Color' , cmap(1,:) )
        set( hlh.get(i), 'Color' , cmap(1,:) )
    end
end
title(['PCN=',num2str(params.pcn)],'FontSize',20);

cb = colorbar('Location','southoutside', 'XTickLabel',{'0','.5', '1'},'XTick', [0,50,100]);
caxis([0 100]);
colormap(cmap);
set(get(cb,'XLabel'),'String','Plasmid fraction','FontSize',14);
set(cb,'Position',[0.65 0.90 0.250 0.025]);% To change size

%%

frac_survivors=zeros(1,nH);

if toPlot==1
    figure('Position', [500 500 700 200])
    clf('reset');set(gcf,'DefaultLineLineWidth',1); set(gcf, 'color', 'white');
    subaxis(1,7,1,1,5,1,'PaddingBottom',0.1)
end

for iy=1:nH
    survivors=0;
    for ix=1:length(freqsT)
        if freqsT(ix)>1-Hs(iy)
            col=cmap(end,:);
            survivors=survivors+1;
        else
            col=[.75 .75 .75];
        end
        if toPlot==1
            plot(ix,Hs(iy),'o','MarkerSize',10,'MarkerFaceColor',col,'MarkerEdgeColor',col); hold on;
        end
    end
    frac_survivors(iy)=survivors/length(freqsT);
end

if toPlot==1
    ylim([-0.15 1.15]);
    ylabel('Dominance','FontSize',16);
    box off
    xticks([]);
    set(gca,'FontSize',14)
    
    % Plot frequency of mutants (horizontal bars)
    subaxis(1,7,6,1,2,1,'PaddingBottom',0.1)
    barh(Hs, frac_survivors,0.8,'FaceColor',[.75 .75 .75]); hold on
    xlabel('Mutant frequency','FontSize',16);
    ylim([-0.15 1.15]);
    xlim([0 1]);
    set(gca,'FontSize',14)
    yticks([]);
    box off
end

end

