%% MGAS_proteome_TP_IPA_plots.m  06032024
clear;
close all hidden;
clc;

%% datapath
datapath = 'F:\APOE-TR-Omics-Data\Single cell type proteomics\Clark analysis-4th round\DEP-count_cutoff=1, prop_cutoff=0.001\IPA';
save_fig_path = 'F:\APOE-TR-Omics-Data\Single cell type proteomics\Clark analysis-4th round\ana fig';
fn = 'AS234-CPUR-sort.xlsx';
% fn = 'MG234-CPUR-sort.xlsx';
cd(datapath);
[a1,b1,c1] = xlsread(fn,'CPsort');
[a2,b2,c2] = xlsread(fn,'URsort');

%% plot 
%% %%% CP %%%%%%%%
z = nan(size(a1,1),3);
if strcmp(fn (1:2),'AS')
z(:,2:3) = a1(:,1:2);
q = a1(:,6:8);
elseif strcmp(fn (1:2),'MG') 
z = a1(:,1:3);
q = a1(:,7:9);
end

zthr = 2;
pthr = 1.3;
xlabels = {'E2' 'E3' 'E4'};
ylabels = b1(2:end,1);
fig_fn = sprintf('%s_CP_zqdotplot.emf',fn(1:5));
temp = z(:);
if strcmp(fn (1:2),'AS')
    colorbar_bin=1;
elseif strcmp(fn (1:2),'MG') 
    colorbar_bin=1;
else
colorbar_bin = round((max(temp(~isnan(temp)))-min(temp(~isnan(temp))))/10);
end
z (q<1.3) = nan;q(q<1.3)=nan;
zpfcp_dotplot_neo(z,q,zthr,pthr,xlabels,ylabels,fig_fn,colorbar_bin);
cd(save_fig_path);
saveas(gcf,fig_fn);





