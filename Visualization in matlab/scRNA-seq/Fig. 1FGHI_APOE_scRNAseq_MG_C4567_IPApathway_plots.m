%% APOE_scRNAseq_MG_C4567_IPApathway_plots.m 

clear;
close all;
clc;

%% vs all
datapath = 'F:\APOE-TR-Omics-Data\Single cell RNA-seq\celltype-results\celltype-markers\IPA-vs all\only positive genes';
cd(datapath);

fn = 'MG-subcluster-CP-sort-posgene.xlsx';
[a1,b1,c1] = xlsread(fn,'C4plot');
figfn = 'C4-Feature pathway.emf';
[a1,b1,c1] = xlsread(fn,'C5plot');
figfn = 'C5-Feature pathway.emf';
[a1,b1,c1] = xlsread(fn,'C6plot');
figfn = 'C6-Feature pathway.emf';
[a1,b1,c1] = xlsread(fn,'C7plot');
figfn = 'C7-Feature pathway.emf';

%% 
cps = b1(2:end,1);
zq = flip(a1(:,[3 5]));

%%
figure ('Position',[0.0010    0.0410    1.5360    0.7488]*1000);
for i = 1:length(cps)
    if zq(i,2)>4 
plot(zq(i,1),i,'.','MarkerSize',40);
    elseif zq(i,2)>3 && zq(i,2)<4
plot(zq(i,1),i,'.','MarkerSize',30);
    elseif zq(i,2)>2 && zq(i,2)<3
plot(zq(i,1),i,'.','MarkerSize',20);
    elseif zq(i,2)>1.3 && zq(i,2)<2
plot(zq(i,1),i,'.','MarkerSize',10);
    end
hold on;
end

x_lim = [0 ceil(max(zq(:,1)))];
y_lim = [0 length(cps)+1];
set(gca,'box','off',"XLim",x_lim,"YLim",y_lim,'XGrid','off','YGrid','off','YTick',1:length(cps),'YTickLabel',flip(cps));
xlabel('CP-z score');

cd(datapath);
saveas(gcf,figfn);


