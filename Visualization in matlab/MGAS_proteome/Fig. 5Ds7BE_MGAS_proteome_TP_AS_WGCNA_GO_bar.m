%% MGAS_proteome_MG_WGCNA_GO_bar.m        WQ 03/19/2024
%% module of interest: red, blue, yellow
%% module gene lists plus log2FC plus FDR also generated
close all;
clear;
clc;

%%
data_path = 'F:\APOE-TR-Omics-Data\Single cell type proteomics\Clark analysis-4th round\WGCNA\AS';
save_fig_path = 'F:\APOE-TR-Omics-Data\Single cell type proteomics\Clark analysis-4th round\ana fig\AS';
save_mat_path = 'F:\APOE-TR-Omics-Data\Single cell type proteomics\Clark analysis-4th round\ana mat\AS';
sig_modules = {'blue' 'yellow' 'red'}; %% 
n_sig_modules = length(sig_modules);

%%
fn = 'AS module-GO-over-representation.xlsx';
for i = 1:n_sig_modules
    cd(data_path);
    [a0,b0,~] = xlsread(fn,char(sig_modules(i)));
    GO_type = b0(2:end,1);
    GO_type_uniq = unique(GO_type);
    n_Go_type = length(GO_type_uniq);
    gene_list_forGO = b0(2:end,9);
    q = a0(:,3);
    GOterms = b0(2:end,3);
    log10_q = -log10(q);

    %%%% plot bar graph for top10
    figure('Position',[427.4000   41.8000  620.8000  740.8000]);
    for j = 1:n_Go_type
        ids = find(contains(GO_type,GO_type_uniq(j)));
        tempq = log10_q(ids);
        [q_neo,ids_neo] = sort(tempq,'descend');
        if length(ids)>=20
            topGO_genelist.(char(GO_type_uniq(j))) = gene_list_forGO(ids(ids_neo(1:20)));
        elseif length(ids)<20
            topGO_genelist.(char(GO_type_uniq(j))) = gene_list_forGO(ids(ids_neo(1:end)));
        end
        tempy = GOterms(ids(ids_neo));
        tempx = tempq(ids_neo);
        subplot(n_Go_type,1,j);
        if length(ids)>=20
            h = bar(flip(tempx(1:15)),0.4);
            xticklabels(flip(tempy(1:15)));

        elseif length(ids)<20
            h = bar(flip(tempx(1:end)),0.4);
            xticklabels(flip(tempy(1:end)));
        end
        hold on;
        ylabel('log10(q)');
        view([90 -90]);
        temp = strcat(sig_modules(i),'-',GO_type_uniq(j));
        title(temp);
    end
    fig_fn = strcat('AS-',sig_modules(i),'-topGOterms.emf');
    cd(save_fig_path);
    saveas(gcf,char(fig_fn));
    cd(save_mat_path);
    temp = strcat('AS-',sig_modules(i),'-topGO_genelist.mat');
    save(char(temp),"topGO_genelist");
end

%% module genes extraction and expression
cd(data_path);
fn = 'geneInfo - AS.csv';
[a1,b1,~] = xlsread(fn);
modulecolors = b1(2:end,2);
genes = b1(2:end,1);
for i = 1:n_sig_modules
temp = sig_modules(i);
ids = find(strcmp(modulecolors,temp));
module_genes.(char(temp)) = genes(ids);
cd(save_mat_path);
    temp = strcat('AS-',sig_modules(i),'.mat');
save(char(temp),"module_genes");
end

%%
cd('F:\APOE-TR-Omics-Data\Single cell type proteomics\Clark analysis-4th round\DEP-count_cutoff=1, prop_cutoff=0.001\log Total.protein\xlsx');
e2fn = 'AS 24Mvs3ME2.xlsx';
e3fn = 'AS 24Mvs3ME3.xlsx';
e4fn = 'AS 24Mvs3ME4.xlsx';

[a2,b2,~] = xlsread(e2fn);
[a3,b3,~] = xlsread(e3fn);
[a4,b4,~] = xlsread(e4fn);

for i = 1:n_sig_modules
temp = module_genes.(char(sig_modules(i)));
[idt2,ida2,idb2] = intersect(temp,b2(2:end,1),'stable');
[idt3,ida3,idb3] = intersect(temp,b3(2:end,1),'stable');
[idt4,ida4,idb4] = intersect(temp,b4(2:end,1),'stable');

temp1 = char(strcat(sig_modules(i),'log2_fcq'));
module_genes.(temp1) = [a2(idb2,[1,5]),a3(idb3,[1,5]),a4(idb4,[1,5])];
cd(save_mat_path);
    temp = strcat('AS-',sig_modules(i),'.mat');
save(char(temp),"module_genes");
end












