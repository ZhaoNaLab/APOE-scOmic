%% APOE_scRNAseq_genes_recruit_Tcells.m
clear;
close all;
clc;

%% gene expression -average across samples
geneexppath = 'F:\APOE-TR-Omics-Data\Single cell RNA-seq\celltype-results\celltype-SCT-average-exp';
cd(geneexppath);
[a0,b0,c0] = xlsread('SCT ave exp.xlsx','AS-SCT-table');
 
clusterlist = b0(2,2:end);
clist = zeros(length(clusterlist),1);
for i=1:length(clusterlist)
    clist(i) = str2num(char(clusterlist(i)));
end

gene_exp = b0(6:end,2:end);
gene_exp0 = zeros(size(gene_exp));
for i=1:size(gene_exp,1)
    for j = 1:size(gene_exp,2)
    gene_exp0(i,j) = str2num(char(gene_exp(i,j)));
    end
end

genelistas = b0(6:end,1);
ngenes = length(genelistas);

ncluster = length(unique(clist));
uniquec = unique(clist);
geneexpas = zeros(ngenes,ncluster);
for i = 1:ncluster
    for j = 1:ngenes
ids = find(clist==uniquec(i));
geneexpas(j,i) = mean(gene_exp0(j,ids));
    end
end

%% genes recruit T cells
gene_recrt = {'Cxcl10','Cxcl14'};

[idt,ida1,idb1] = intersect(gene_recrt,genelistas,'stable');
[x,y] = setdiff(gene_recrt,idt);
idtexp = geneexpas(idb1,:);
clustername = {'C0','C1','C2','C3','C4','C5','C6'};
cgo  = clustergram(idtexp,'Standardize','row','colormap',redbluecmap);  %% ,'cluster','column'
set(cgo,'RowLabels',idt ,'ColumnLabels',clustername);

%% gene expression -average across samples in microglia
geneexppath = 'F:\APOE-TR-Omics-Data\Single cell RNA-seq\celltype-results\celltype-SCT-average-exp';
cd(geneexppath);
[a0,b0,c0] = xlsread('SCT ave exp.xlsx','MG-SCT-table');
 
clusterlist = b0(2,2:end);
clist = zeros(length(clusterlist),1);
for i=1:length(clusterlist)
    clist(i) = str2num(char(clusterlist(i)));
end

gene_exp = b0(6:end,2:end);
gene_exp0 = zeros(size(gene_exp));
for i=1:size(gene_exp,1)
    for j = 1:size(gene_exp,2)
    gene_exp0(i,j) = str2num(char(gene_exp(i,j)));
    end
end

genelist = b0(6:end,1);
ngenes = length(genelist);

clistmg = [0 1 2 3 4 5 6 7];
ncluster = length(clistmg);
geneexpmg = zeros(ngenes,ncluster);

for i = 1:ncluster
    for j = 1:ngenes
ids = find(clist==clistmg(i));
geneexpmg(j,i) = mean(gene_exp0(j,ids));
    end
end

%% genes recruit T cells
gene_recrt = {'Cxcl10','Cxcl16'};

[idt,ida1,idb1] = intersect(gene_recrt,genelist,'stable');
[x,y] = setdiff(gene_recrt,idt);
idtexp = geneexpmg(idb1,:);
clustername = {'C0','C1','C2','C3','C4','C5','C6','C7'};
cgo  = clustergram(idtexp,'Standardize','row','colormap',redbluecmap);  %% ,'cluster','column'
set(cgo,'RowLabels',idt ,'ColumnLabels',clustername);



