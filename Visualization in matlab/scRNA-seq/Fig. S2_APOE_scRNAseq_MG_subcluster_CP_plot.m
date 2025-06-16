%% APOE_scRNAseq_MG_subcluster_CP_plot.m  05312024
%% plots from IPA results
clear;
close all;
clc;

%% datapath
datapath = 'F:\APOE-TR-Omics-Data\Single cell RNA-seq\celltype-aggregate-deg-results\SEX+APOExAGE-selected\IPA';
save_fig_path = 'F:\APOE-TR-Omics-Data\Single cell RNA-seq\celltype-aggregate-deg-results\SEX+APOExAGE-selected\ana_fig';

%% extract data
cd(datapath);
fn = 'MG0-subcluster-CPURDF-sort.xlsx';
% fn = 'MG1-subcluster-CPURDF-sort.xlsx';
% fn = 'MG2-subcluster-CPURDF-sort.xlsx';
% fn = 'MG3-subcluster-CPURDF-sort.xlsx';
% fn = 'MG4-subcluster-CPURDF-sort.xlsx';
% fn = 'MG5-subcluster-CPURDF-sort.xlsx';
% fn = 'MG6-subcluster-CPURDF-sort.xlsx';
% fn = 'MG7-subcluster-CPURDF-sort.xlsx';

[a1,b1,c1] = xlsread(fn,'CPsort');


% [a1,b1,c1] = xlsread(fn,'CPcustom');
% [a2,b2,c2] = xlsread(fn,'URcustom');

%% 
% [a3,b3,c3] = xlsread(fn,'DFcustom');

%% plots
%% %%% CP %%%%%%%%
% z = a1(:,1:3);
% q = a1(:,7:9);
z = a1(:,2:4);
q = a1(:,8:10);

zthr = 2;
pthr = 1.3;
xlabels = {'E2' 'E3' 'E4'};
ylabels = b1(2:end,1);
fig_fn = sprintf('%s_CP_zqdotplot.emf',fn(1:3));
temp = z(:);
if strcmp(fn (1:3),'MG2')
    colorbar_bin=0.5;
elseif strcmp(fn (1:3),'MG3') || strcmp(fn (1:3),'MG7')
        colorbar_bin=1;
else
colorbar_bin = round((max(temp(~isnan(temp)))-min(temp(~isnan(temp))))/10);
end
colorbar_bin =1;
z (q<1.3) = nan;q(q<1.3)=nan;
zpfcp_dotplot_neo(z,q,zthr,pthr,xlabels,ylabels,fig_fn,colorbar_bin);
cd(save_fig_path);
saveas(gcf,fig_fn);

%% %%% UR %%%%%%%%
if strcmp(fn (1:3),'MG7')
z = a2(:,1:2);
p = a2(:,4:5);
else
z = a2(:,1:3);
p = a2(:,4:6);  
end
zthr = 2;
pthr = 1.3;
xlabels = {'E2' 'E3' 'E4'};
ylabels = b2(2:end,1);
fig_fn = sprintf('%s_UR_zqdotplot.emf',fn(1:3));
colorbar_bin = 1;
z (p<1.3) = nan;p(p<1.3)=nan;
zpfcp_dotplot_neo(z,p,zthr,pthr,xlabels,ylabels,fig_fn,colorbar_bin);
 
cd(save_fig_path);
saveas(gcf,fig_fn);

%% %%% DF %%%%%%%%
z = a3(:,1:3);
p = a3(:,4:6);
zthr = 2;
pthr = 1.3;
xlabels = {'E2' 'E3' 'E4'};
ylabels = b3(2:end,1);
fig_fn = sprintf('%s_DF_zqdotplot.emf',fn(1:3));
colorbar_bin = 0.5;
z (p<1.3) = nan;p(p<1.3)=nan;
zpfcp_dotplot_neo(z,p,zthr,pthr,xlabels,ylabels,fig_fn,colorbar_bin);
cd(save_fig_path);
saveas(gcf,fig_fn);


