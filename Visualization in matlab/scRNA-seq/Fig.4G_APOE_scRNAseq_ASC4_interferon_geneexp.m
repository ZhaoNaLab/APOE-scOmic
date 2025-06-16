%% APOE_scRNAseq_ASC4_interferon_geneexp.m%%
clear;
close all;
clc;

%% gene expression -average across samples
cd('F:\APOE-TR-Omics-Data\Single cell RNA-seq\anamat');
load('apoescOmic_assct.mat');
genelist = assctexp.genelist;
geneexp = assctexp.geneexp;

%% degs in c4
datapath = 'F:\APOE-TR-Omics-Data\Single cell RNA-seq\celltype-aggregate-deg-results\SEX+APOExAGE-selected\IPA';
cd(datapath);
fn = 'AS4-E4-interferon-gene_fcfdr.xls';
[a,b,c] = xlsread(fn);
ifngenes = b(3:end,5);
expfc = a(:,[1,3,5]);
expfdr = a(:,[1,3,5]+1);
nifngenes = length(ifngenes);

%% fc plot 
z = expfc;
q = -log10(expfdr);

zthr = 2;
pthr = 1.3;
xlabels = {'E2' 'E3' 'E4'};
ylabels =ifngenes;
fig_fn = sprintf('%s_interferon_fcqdotplot.emf',fn(1:3));
colorbar_bin=0.5;

z (q<1.3) = nan;q(q<1.3)=nan;
zpfcp_dotplot_neo(z,q,zthr,pthr,xlabels,ylabels,fig_fn,colorbar_bin);
cd(datapath);
saveas(gcf,fig_fn);

%% expression 
[idt,ida,idb] = intersect(ifngenes,genelist,'stable');
c4id = find(strcmp(assctexp.clusterlist,'4')); %% mouseinfo, _2_3_F,_2_3_M,_2_24_F,
% _2_24_M,_3_3_F,_3_3_M,_3_24_F,_3_24_M,_4_3_F,_4_3_M,_4_24_F,_4_24_M
ifnexp1 = geneexp(idb,c4id);
ifnexp2 = zeros(nifngenes,6);
for i = 1:nifngenes
    for j = 1:6
        ifnexp2(i,j) = mean(ifnexp1(i,(j-1)*2+1:(j-1)*2+2));
    end
end

%% heatmap
xlabels = {'E23' 'E224' 'E33' 'E324' 'E43' 'E424'};
idtexpz = zscore(ifnexp2');
figure('Position',[488.0000   41.8000  560.0000  740.8000]);
h = heatmap(xlabels,ifngenes,idtexpz');
% colormap(gca,bluewhitered(256)); 
% colormap redbluecmap;
colormap(linspecer);
% colormap(brewermap([],"BrBG"));
% colormap(brewermap([],"PRGn"));
% colormap(brewermap([],"PiYG"));
colormap(flip(brewermap([],"RdGY")));
% colormap(brewermap([],"PuOr"));
% colormap(flip(brewermap([],"RdBu"))); %% the same as redbluecmap

fn =  'AS-C4-interferon-geneexpheatmap.emf';
cd(datapath);
saveas(gcf,fn);

%% all degs in Cluster 4
datapath = 'F:\APOE-TR-Omics-Data\Single cell RNA-seq\celltype-aggregate-deg-results\SEX+APOExAGE-selected\DEG-AS-subcluster-xlsx';
cd(datapath);
fn1 = 'AS-4 24Mvs3ME2.xlsx';
fn2 = 'AS-4 24Mvs3ME3.xlsx';
fn3 = 'AS-4 24Mvs3ME4.xlsx';

[a1,b1,c1] = xlsread(fn1);
[a2,b2,c2] = xlsread(fn2);
[a3,b3,c3] = xlsread(fn3);

%% deg summary
e2deg = b1(find(a1(:,5)<0.05)+1,1);
e3deg = b2(find(a2(:,5)<0.05)+1,1);
e4deg = b3(find(a3(:,5)<0.05)+1,1);

degs = union_several(e2deg,e3deg,e4deg); %% 118
ndegs = numel(degs);
%% expression 
[idt,ida,idb] = intersect(degs,genelist,'stable');
c4id = find(strcmp(assctexp.clusterlist,'4')); %% mouseinfo, _2_3_F,_2_3_M,_2_24_F,
% _2_24_M,_3_3_F,_3_3_M,_3_24_F,_3_24_M,_4_3_F,_4_3_M,_4_24_F,_4_24_M
degexp1 = geneexp(idb,c4id);
degexp2 = zeros(ndegs,6);
for i = 1:ndegs
    for j = 1:6
        degexp2(i,j) = mean(degexp1(i,(j-1)*2+1:(j-1)*2+2));
    end
end

%% find exp = 0;
zeroids = find(any(degexp2==0,2));
rikgeneids = find(contains(degs,'Rik'));
Gmgeneids = find(contains(degs,'Gm'));
neoids = setdiff(1:ndegs,[zeroids; rikgeneids; Gmgeneids]);
degneo = degs(neoids);
degexp2neo = degexp2(neoids,:);

%% expression norm
degexp2norm = zeros(size(degexp2neo));
for i = 1:numel(degneo)
degexp2norm(i,:) = [degexp2neo(i,1:2)./degexp2neo(i,1)...
    degexp2neo(i,3:4)./degexp2neo(i,3)...
    degexp2neo(i,5:6)./degexp2neo(i,5)];
end
%%
groups = {'e23','e224','e33','e324','e43','e424'};
cgo  = clustergram(degexp2norm','Standardize','column','colormap',redbluecmap,'cluster','row');  %% ,'cluster','column'
set(cgo,'RowLabels',groups,'ColumnLabels',degneo);



