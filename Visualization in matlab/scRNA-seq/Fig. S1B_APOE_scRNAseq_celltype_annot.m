clear;
close all;
clc;

%% 
cd('F:\APOE-TR-Omics-Data\Single cell RNA-seq\anamat');
load('apoescOmic_combinedexp.mat');

%% datapath  & gene list file
datapath  = 'F:\APOE-TR-Omics-Data\Single cell RNA-seq\celltype-results\0-Celltypeannotation';
save_fig_path = datapath;
fn = 'Celltypemarkerlist.xlsx';
cd(datapath);
%% 
[~,b1,~] = xlsread(fn,'Final');
celltype = b1(2:end,1);
[idt,ida1,idb1] = intersect(flip(celltype),cmbsctexp.genelist,'stable');
[x,y] = setdiff(celltype,idt);
idtexp = cmbsctexp.geneexp_ave(idb1,:);
celltypename = cmbsctexp.uniquecelltype;
% cgo  = celltypegram(idtexp,'Standardize','row','colormap',redbluecmap,'celltype','row');  %% ,'celltype','column'
% set(cgo,'RowLabels',idt,'ColumnLabels',celltypename);
%
idtexpz = zscore(idtexp');
figure('Position',[488.0000   41.8000  560.0000  740.8000]);
genename = flip(idt);
geneexpvalues = flip(idtexpz');
h = heatmap(celltypename,flip(idt),flip(idtexpz'));

figfn =  'Celltype-anno-ave.emf';
cd(save_fig_path);
saveas(gcf,figfn);
