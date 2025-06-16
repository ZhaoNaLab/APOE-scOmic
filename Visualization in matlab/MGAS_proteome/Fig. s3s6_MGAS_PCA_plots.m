%% MGAS_proteome_TP_PCA.m
clear;
close all hidden;
clc;

%% dataapath
datapath = 'F:\APOE-TR-Omics-Data\Single cell type proteomics\Clark analysis-4th round';
save_fig_path = 'F:\APOE-TR-Omics-Data\Single cell type proteomics\Clark analysis-4th round\ana fig';

% fn = 'AS Tube,Sex,Total.protein residuals.csv';
fn = 'MG Tube,Sex,Total.protein residuals.csv';

%% extract data
cd(datapath);
[a1,b1,c1] = xlsread(fn);

prots = b1; 
exp = a1(2:end,:);
expid = a1(1,:);

%% mouseinfo
cd('F:\APOE-TR-Omics-Data\Single cell type proteomics\Ana_mat');
load asexp_ori.mat;
mouseinfo = asexp.mouseinfo;

%%
[~,scores,pvars] = pca(exp');
x = zscore(scores(:,1));
y = zscore(scores(:,2));
z = zscore(scores(:,3));

PC1 = num2str(round(pvars(1)/sum(pvars)*100));
PC2 = num2str(round(pvars(2)/sum(pvars)*100));
PC3 = num2str(round(pvars(3)/sum(pvars)*100));

%% plot PC1 PC2;
figure ('position',[73.8000  241.0000  817.6000  524.0000]);
fig_fn = sprintf('%s_PCA2d.emf',fn(1:2));
h = subplot (231);
gscatter(h,x,y,mouseinfo.Genotype,'rgb','...',15);
title('apoe');
xlabel(['PC1(' PC1 '%)']);ylabel(['PC2(' PC2 '%)']);
h = subplot (232);
gscatter(h,x,y,mouseinfo.Group,'rg','..',15);
title('age');
xlabel(['PC1(' PC1 '%)']);ylabel(['PC2(' PC2 '%)']);

h = subplot (233);
gscatter(h,x,y,mouseinfo.Sex,'rg','..',15);
title('sex');
xlabel(['PC1(' PC1 '%)']);ylabel(['PC2(' PC2 '%)']);

h=subplot (234);
gscatter(h,x,y,mouseinfo.Tube,'rg','..',15);
title('tube');
xlabel(['PC1(' PC1 '%)']);ylabel(['PC2(' PC2 '%)']);

h=subplot (235);
gscatter(h,x,y,mouseinfo.Batch,'rgbmck','......',15);
title('isolaiton date');
xlabel(['PC1(' PC1 '%)']);ylabel(['PC2(' PC2 '%)']);

cd(save_fig_path);
saveas(gcf,fig_fn);

%% % 3d
figure ('position',[.0738    0.2410    1.3392    0.5240]*1000);
fig_fn = sprintf('%s_PCA.emf',fn(1:2));

subplot (231);
gscatter3(x,y,z,mouseinfo.Genotype',{'r','g','b'},{'.','.','.'},15);

title('apoe');
xlabel(['PC1(' PC1 '%)']);ylabel(['PC2(' PC2 '%)']);zlabel(['PC3(' PC3 '%)']);

subplot (232);
gscatter3(x,y,z,mouseinfo.Group',{'r','b'},{'.','.'},15);
title('age');
xlabel(['PC1(' PC1 '%)']);ylabel(['PC2(' PC2 '%)']);zlabel(['PC3(' PC3 '%)']);

subplot (233);
gscatter3(x,y,z,mouseinfo.Sex',{'r','b'},{'.','.'},15);
title('sex');
xlabel(['PC1(' PC1 '%)']);ylabel(['PC2(' PC2 '%)']);zlabel(['PC3(' PC3 '%)']);

subplot (234);
gscatter3(x,y,z,mouseinfo.Tube',{'r','b'},{'.','.'},15);
title('tube');
xlabel(['PC1(' PC1 '%)']);ylabel(['PC2(' PC2 '%)']);zlabel(['PC3(' PC3 '%)']);

subplot (235);
gscatter3(x,y,z,mouseinfo.Batch',{'r','b','g','m','c','k'},{'.','.','.','.','.','.'},15);
title('isolaiton date');
xlabel(['PC1(' PC1 '%)']);ylabel(['PC2(' PC2 '%)']);zlabel(['PC3(' PC3 '%)']);

cd(save_fig_path);
saveas(gcf,fig_fn);

%% 
e2id = find(strcmp(mouseinfo.Genotype,'E2'));
e3id = find(strcmp(mouseinfo.Genotype,'E3'));
e4id = find(strcmp(mouseinfo.Genotype,'E4'));
e2exp = exp(:,e2id);
e3exp = exp(:,e3id);
e4exp = exp(:,e4id);

%% pca E2
[~,scores,pvars] = pca(e2exp');
x = zscore(scores(:,1));
y = zscore(scores(:,2));
z = zscore(scores(:,3));

PC1 = num2str(round(pvars(1)/sum(pvars)*100));
PC2 = num2str(round(pvars(2)/sum(pvars)*100));
PC3 = num2str(round(pvars(3)/sum(pvars)*100));

%% plot PC1 PC2;
figure ('position',[73.8000  241.0000  528.0000  524.0000]);
fig_fn = sprintf('%s_PCA2d_E2.emf',fn(1:2));

h = subplot (221);
gscatter(h,x,y,mouseinfo.Group(e2id),'rg','..',15);
title('age');
xlabel(['PC1(' PC1 '%)']);ylabel(['PC2(' PC2 '%)']);

h = subplot (222);
gscatter(h,x,y,mouseinfo.Sex(e2id),'rg','..',15);
title('sex');
xlabel(['PC1(' PC1 '%)']);ylabel(['PC2(' PC2 '%)']);

h=subplot (223);
gscatter(h,x,y,mouseinfo.Tube(e2id),'rg','..',15);
title('tube');
xlabel(['PC1(' PC1 '%)']);ylabel(['PC2(' PC2 '%)']);

h=subplot (224);
gscatter(h,x,y,mouseinfo.Batch(e2id),'rgbmck','......',15);
title('isolaiton date');
xlabel(['PC1(' PC1 '%)']);ylabel(['PC2(' PC2 '%)']);

cd(save_fig_path);
saveas(gcf,fig_fn);

%% %
figure ('position',[73.8000  241.0000  905.6000  524.0000]);
fig_fn = sprintf('%s_PCA_E2.emf',fn(1:2));

subplot (221);
gscatter3(x,y,z,mouseinfo.Group(e2id)',{'r','b'},{'.','.'},15);
title('age');
xlabel(['PC1(' PC1 '%)']);ylabel(['PC2(' PC2 '%)']);zlabel(['PC3(' PC3 '%)']);

subplot (222);
gscatter3(x,y,z,mouseinfo.Sex(e2id)',{'r','b'},{'.','.'},15);
title('sex');
xlabel(['PC1(' PC1 '%)']);ylabel(['PC2(' PC2 '%)']);zlabel(['PC3(' PC3 '%)']);

subplot (223);
gscatter3(x,y,z,mouseinfo.Tube(e2id)',{'r','b'},{'.','.'},15);
title('tube');
xlabel(['PC1(' PC1 '%)']);ylabel(['PC2(' PC2 '%)']);zlabel(['PC3(' PC3 '%)']);

subplot (224);
gscatter3(x,y,z,mouseinfo.Batch(e2id)',{'r','b','g','m','c','k'},{'.','.','.','.','.','.'},15);
title('isolaiton date');
xlabel(['PC1(' PC1 '%)']);ylabel(['PC2(' PC2 '%)']);zlabel(['PC3(' PC3 '%)']);

cd(save_fig_path);
saveas(gcf,fig_fn);

%% pca E3
[coef,scores,pvars] = pca(e3exp');
x = zscore(scores(:,1));
y = zscore(scores(:,2));
z = zscore(scores(:,3));

PC1 = num2str(round(pvars(1)/sum(pvars)*100));
PC2 = num2str(round(pvars(2)/sum(pvars)*100));
PC3 = num2str(round(pvars(3)/sum(pvars)*100));

%% plot PC1 PC2;
figure ('position',[73.8000  241.0000  528.0000  524.0000]);
fig_fn = sprintf('%s_PCA2d_E3.emf',fn(1:2));

h = subplot (221);
gscatter(h,x,y,mouseinfo.Group(e3id),'rg','..',15);
title('age');
xlabel(['PC1(' PC1 '%)']);ylabel(['PC2(' PC2 '%)']);

h = subplot (222);
gscatter(h,x,y,mouseinfo.Sex(e3id),'rg','..',15);
title('sex');
xlabel(['PC1(' PC1 '%)']);ylabel(['PC2(' PC2 '%)']);

h=subplot (223);
gscatter(h,x,y,mouseinfo.Tube(e3id),'rg','..',15);
title('tube');
xlabel(['PC1(' PC1 '%)']);ylabel(['PC2(' PC2 '%)']);

h=subplot (224);
gscatter(h,x,y,mouseinfo.Batch(e3id),'rgbmck','......',15);
title('isolaiton date');
xlabel(['PC1(' PC1 '%)']);ylabel(['PC2(' PC2 '%)']);

cd(save_fig_path);
saveas(gcf,fig_fn);

%% %
figure ('position',[73.8000  241.0000  905.6000  524.0000]);
fig_fn = sprintf('%s_PCA_E3.emf',fn(1:2));

subplot (221);
gscatter3(x,y,z,mouseinfo.Group(e3id)',{'r','b'},{'.','.'},15);
title('age');
xlabel(['PC1(' PC1 '%)']);ylabel(['PC2(' PC2 '%)']);zlabel(['PC3(' PC3 '%)']);

subplot (222);
gscatter3(x,y,z,mouseinfo.Sex(e3id)',{'r','b'},{'.','.'},15);
title('sex');
xlabel(['PC1(' PC1 '%)']);ylabel(['PC2(' PC2 '%)']);zlabel(['PC3(' PC3 '%)']);

subplot (223);
gscatter3(x,y,z,mouseinfo.Tube(e3id)',{'r','b'},{'.','.'},15);
title('tube');
xlabel(['PC1(' PC1 '%)']);ylabel(['PC2(' PC2 '%)']);zlabel(['PC3(' PC3 '%)']);

subplot (224);
gscatter3(x,y,z,mouseinfo.Batch(e3id)',{'r','b','g','m','c','k'},{'.','.','.','.','.','.'},15);
title('isolaiton date');
xlabel(['PC1(' PC1 '%)']);ylabel(['PC2(' PC2 '%)']);zlabel(['PC3(' PC3 '%)']);

cd(save_fig_path);
saveas(gcf,fig_fn);


%% pca E4
[~,scores,pvars] = pca(e4exp');
x = zscore(scores(:,1));
y = zscore(scores(:,2));
z = zscore(scores(:,3));

PC1 = num2str(round(pvars(1)/sum(pvars)*100));
PC2 = num2str(round(pvars(2)/sum(pvars)*100));
PC3 = num2str(round(pvars(3)/sum(pvars)*100));

%% plot PC1 PC2;
figure ('position',[73.8000  241.0000  528.0000  524.0000]);
fig_fn = sprintf('%s_PCA2d_E4.emf',fn(1:2));

h = subplot (221);
gscatter(h,x,y,mouseinfo.Group(e4id),'rg','..',15);
title('age');
xlabel(['PC1(' PC1 '%)']);ylabel(['PC2(' PC2 '%)']);

h = subplot (222);
gscatter(h,x,y,mouseinfo.Sex(e4id),'rg','..',15);
title('sex');
xlabel(['PC1(' PC1 '%)']);ylabel(['PC2(' PC2 '%)']);

h=subplot (223);
gscatter(h,x,y,mouseinfo.Tube(e4id),'rg','..',15);
title('tube');
xlabel(['PC1(' PC1 '%)']);ylabel(['PC2(' PC2 '%)']);

h=subplot (224);
gscatter(h,x,y,mouseinfo.Batch(e4id),'rgbmck','......',15);
title('isolaiton date');
xlabel(['PC1(' PC1 '%)']);ylabel(['PC2(' PC2 '%)']);

cd(save_fig_path);
saveas(gcf,fig_fn);



%% %
figure ('position',[73.8000  241.0000  905.6000  524.0000]);
fig_fn = sprintf('%s_PCA_E4.emf',fn(1:2));

subplot (221);
gscatter3(x,y,z,mouseinfo.Group(e4id)',{'r','b'},{'.','.'},15);
title('age');
xlabel(['PC1(' PC1 '%)']);ylabel(['PC2(' PC2 '%)']);zlabel(['PC3(' PC3 '%)']);

subplot (222);
gscatter3(x,y,z,mouseinfo.Sex(e4id)',{'r','b'},{'.','.'},15);
title('sex');
xlabel(['PC1(' PC1 '%)']);ylabel(['PC2(' PC2 '%)']);zlabel(['PC3(' PC3 '%)']);

subplot (223);
gscatter3(x,y,z,mouseinfo.Tube(e4id)',{'r','b'},{'.','.'},15);
title('tube');
xlabel(['PC1(' PC1 '%)']);ylabel(['PC2(' PC2 '%)']);zlabel(['PC3(' PC3 '%)']);

subplot (224);
gscatter3(x,y,z,mouseinfo.Batch(e4id)',{'r','b','g','m','c','k'},{'.','.','.','.','.','.'},15);
title('isolaiton date');
xlabel(['PC1(' PC1 '%)']);ylabel(['PC2(' PC2 '%)']);zlabel(['PC3(' PC3 '%)']);

cd(save_fig_path);
saveas(gcf,fig_fn);














