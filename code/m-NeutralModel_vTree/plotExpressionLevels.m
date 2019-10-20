
clc
clear all
close all

run('lib/addpath_recurse');
addpath_recurse('lib/');
addpath_recurse('src/');

%% EXPERIMENTAL DATA

xy_gyrA=[0.0857142857142857, 4];
xy_folA=[0.87432498772705, 4];


%% LEVELS SETS
 
figure(1); clf('reset'); set(gcf,'DefaultLineLineWidth',2); set(gcf, 'color', 'white'); 
subaxis(1,1,1,'PaddingRight',0.05,'PaddingBottom',0.075,'PaddingLeft',.01);
set(gcf,'Units','Pixels','Position',[100         918         520 400])
 
maxPCN=20;
hs=linspace(0,1.2,1000);
levels=(0:.1:1);

CT=flipud([1 1 1; cbrewer('seq', 'OrRd', length(levels)-1)]);

area([0.01 max(hs)],[maxPCN+1 maxPCN+1],'FaceColor',CT(1,:),'EdgeColor',CT(1,:)); hold on;
 for ilevel=1:length(levels)-1
    k=1-levels(ilevel);
    PCNs=k./hs;
    area(hs, PCNs,'FaceColor',CT(ilevel+1,:), 'LineWidth',1,'EdgeColor',CT(ilevel+1,:)); hold on;
 end
 
set(gca,'XScale','log');
xticks([.001 .01 .1 1]);
ylim([1,maxPCN+1]);
xlim([.005, max(hs)]);
yticks([1 5 10 15 20])
set(gca,'fontsize',18); 
xlabel('log(h)','FontSize',24);
ylabel('PCN','FontSize',24);
 
colormap(flipud(CT));
h= colorbar;
ylabel(h,'% of phenotypic resistance','FontSize',24)
set(h,'ytick',linspace(0,1,12)+.05);
ytl=num2cell(linspace(0,100,11));
ytl{end}=sprintf(['\x2265 100']);
set(h,'yticklabel',ytl)

% PLOT DATA POINTS
plot(xy_gyrA(1),xy_gyrA(2), 'ok','MarkerFaceColor','k');
text(xy_gyrA(1)+.01,xy_gyrA(2)+.5, '{\it gyrA}','HorizontalAlign','center','VerticalAlign','bottom','FontSize',20,'Color','k');
plot(xy_folA(1),xy_folA(2), 'ow','MarkerFaceColor','w');
text(xy_folA(1)+.01,xy_folA(2)+.5, '{\it folA}','HorizontalAlign','right','VerticalAlign','bottom','FontSize',20,'Color','w');


export_fig 'figures/resistance_logHvPCN.pdf'

