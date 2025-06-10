library(dplyr)
library(magrittr)
library(stringr)
library(rlang)

library(car)
library(multcomp)
library(nlme)

library(ggplot2)
library(patchwork)
library(broman)
library(RColorBrewer)

library(limma)
library(SummarizedExperiment)

load_excel <- function() {
  pr = list(AS = readxl::read_xlsx(path='Datasets/Proteomics/Astrocyte/Astrocyte proteomics data_12032021.xlsx', skip=3),
            MG = readxl::read_xlsx(path='Datasets/Proteomics/Microglia/Microglia proteomics data_11182021.xlsx', skip=3))
  colnames(pr[['MG']])[1] = colnames(pr[['AS']])[1] = 'ProteinGroup'
  pr.meta = readxl::read_xlsx('Datasets/Proteomics/ApoE-MG-AS-Proteomics-Mice-Info.xlsx', .name_repair = 'universal') # remove spaces from column names
  
  # sanity check
  stopifnot(colnames(pr[['MG']][9:44]) %in% pr.meta$Mouse.. & 
            pr.meta$Mouse.. %in% colnames(pr[['MG']][9:44]) &
            colnames(pr[['AS']][9:44]) %in% pr.meta$Mouse.. & 
            pr.meta$Mouse.. %in% colnames(pr[['AS']][9:44]))
  return(pr)
}
get_accessions <- function(x) {
  tidyr::separate(data = x, col = 'ProteinAccession', sep='\\|', into = c(NA, 'Accession', NA), remove = FALSE)
}
merge_gene_name <- function(x) {
  # although we normally wouldn't filter the proteomics UNIPROT to just id's we could map to a unique gene,
  # i'm choosing to do so because this would enable us to integrate with rna-seq from the start
  gene_names = read.table('Datasets/Proteomics/gene-names.tsv', sep='\t', header = TRUE)
  x %>% get_accessions %>% left_join(gene_names) %>% relocate(GeneName, .after=GN) %>%
    mutate(GeneNameUpdated = ifelse(is.na(GeneName), GN, GeneName)) %>%
    relocate(GeneNameUpdated, .after=GeneName)
}
load_pr <- function() {
  pr = load_excel()
  pr.meta = readxl::read_xlsx('Datasets/Proteomics/ApoE-MG-AS-Proteomics-Mice-Info.xlsx', .name_repair = 'universal')
  for (celltype in c('AS', 'MG')) {
    pr[[celltype]] <- merge_gene_name(pr[[celltype]])
    
    pr[[celltype]] <- pr[[celltype]] %>% 
      filter(GeneNameUpdated != 'NA' & !is.na(GeneNameUpdated) & GeneNameUpdated != '') %>%
      group_by(GeneNameUpdated) %>%
      arrange(desc(PSMs)) %>%
      slice_head(n=1) %>% # After running DEP, I find that the logFC, AveExpr, are nearly the same for proteins that map to same gene. Therefore I'm keeping the protein with the most PSMs
      as.data.frame()
    rownames(pr[[celltype]]) = pr[[celltype]]$GeneNameUpdated
    dat = pr[[celltype]] %>% dplyr::select(-(ProteinGroup:`Coverage %`))
    dat <- dat[, match(pr.meta$Mouse.., colnames(dat))]
    dat.gene = pr[[celltype]] %>% dplyr::select(ProteinGroup:`Coverage %`)
    pr.meta$Total.protein <- colSums(dat)
    pr[[celltype]] = SummarizedExperiment(list(counts=dat),
                                          rowData=dat.gene, 
                                          colData=pr.meta)
  }
  return(pr)
}

get_normalized <- function(celltype, pr = load_pr(), normalization.method = 'log Total.protein') {
  dat = assay(pr[[celltype]])
  if (normalization.method == 'log') {
    return(log(dat))
  } else if (normalization.method == 'median') {
    return(log(dat / matrixStats::colMedians(as.matrix(dat))))
  } else if (normalization.method == 'log Total.protein') {
    return(log(dat / colData(pr[[celltype]])$Total.protein))
  }
  else {
    return(limma::normalizeBetweenArrays(log(dat), method = normalization.method))
  }
}
get_outdir <- function(normalization.method) {
  if (normalization.method == 'scale') {
    return(paste0('limma-trend ', normalization.method))
  } else {
    return(normalization.method)
  } 
}
load_dge <- function(celltype, pr=NULL) {
  if (is.null(pr)) {
    pr = load_pr()
  }
  y = DGEList(counts=assay(pr[[celltype]]), 
              norm.factors = rep(1, ncol(pr[[celltype]])), # by setting offsets to 0, we don't compute new library sizes ..?
              samples = colData(pr[[celltype]]),
              group = paste0(pr[[celltype]]$Genotype, '.', pr[[celltype]]$Group))
  # Something to be aware of -- colData(pr) has spaces in colnames
  # y$samples does not
  y$samples = load_protein_sample_info(y$samples)
  return(y)
}

load_protein_sample_info <- function(si = readxl::read_xlsx('Datasets/Proteomics/ApoE-MG-AS-Proteomics-Mice-Info.xlsx', .name_repair = 'universal')) {
  si <- as.data.frame(si, row.names=rownames(si))
  si %<>% dplyr::rename(
    Age = Group,
    APOE = Genotype
  )
  si$Age = factor(si$Age, c('Young', 'aged'))
  levels(si$Age) = c('3M', '24M')
  si %<>% mutate(
    Isolation.Date = factor(Isolation.Date),
    APOE = factor(APOE, c('E3', 'E2', 'E4')),
    Sex = factor(Sex, c('M', 'F')),
    group = paste0(APOE, '_', Age),
  ) %>% mutate(
    Isolation.Batch = case_when(
      Isolation.Date %in% paste0('2021-07-2',0:1) ~ 'A',
      Isolation.Date %in% paste0('2021-07-2',2:3) ~ 'B',
      Isolation.Date %in% paste0('2021-07-2',4:5) ~ 'C'     
      #Isolation.Date %in% c('2021-07-20') ~ 'A',
      #Isolation.Date %in% paste0('2021-07-2', 1:4) ~ 'B',
      #Isolation.Date %in% c('2021-07-25') ~ 'C',
    )
  ) %>% mutate(Isolation.Batch = factor(Isolation.Batch)) 
  return(si)
}

plot_gene <- function(x, gene, si=load_protein_sample_info()) {
  x <- x[gene,]
  x %>%
    as.data.frame %>%
    tibble::rownames_to_column('gene') %>%
    tidyr::pivot_longer(-gene, names_to='sample', values_to='value') %>%
    left_join(si %>% dplyr::select(Tube, APOE, Sex, group, Mouse..) %>%
              mutate(Mouse.. = as.character(Mouse..)) %>%
              dplyr::rename(sample = Mouse..)) %>%
    ggplot(aes(x = APOE,
               y=value,
               fill = group,
               #group = group,
               shape = Sex)) + 
    geom_jitter(position=position_jitterdodge(jitter.width=0.1, dodge.width=1),
                size = 2) + 
    facet_wrap(~gene, scales='free_y') +
    scale_shape_manual(values=c(21,25)) + 
    scale_fill_manual(values = group.colors) + 
    #theme(axis.text.x = element_blank()) + 
    guides(fill = guide_legend(override.aes = list(shape = 21)))
}
plot_pca <- function(pca, group, si = load_protein_sample_info(), pcs=c(1,2)) {
  pc1.var = round((pca$eig0 / pca$totvar0)[pcs[1]], 3)
  pc2.var = round((pca$eig0 / pca$totvar0)[pcs[2]], 3)
  PC1 = paste0('PC', pcs[1])
  PC2 = paste0('PC', pcs[2])
  cbind(pca$scores, si) %>%
#    as.data.frame %>%
    ggplot(aes_string(x=eval_tidy(PC1), y=eval_tidy(PC2), fill=eval_tidy(group))) +
    #  scale_color_viridis_c() +
    geom_point(size=3, shape=21) +
    xlab(paste0('PC', pcs[1], ' (', pc1.var*100, '%)')) +
    ylab(paste0('PC', pcs[2], ' (', pc2.var*100, '%)')) 
}
# ------ colors -------
group.colors = brewer.pal(n=6, name='Paired')
names(group.colors) = c('E2_3M', 'E2_24M', 'E3_3M', 'E3_24M', 'E4_3M', 'E4_24M')

plot_gene_heatmap <- function(x, genes, si=load_protein_sample_info(), ...) {
  x <- x[genes,]
  
  #dat = list(E2 = as.matrix(x[, si$APOE == 'E2']),
  #           E3 = as.matrix(x[, si$APOE == 'E3']),
  #           E4 = as.matrix(x[, si$APOE == 'E4']))
  #for (g in names(dat)) {
  #  dat[[g]] <- t(scale(t(dat[[g]])))
  #}
  #x <- bind_cols(dat$E2, dat$E3, dat$E4)
  #x <- x[, match(rownames(si), colnames(x))]
  rownames(x) <- genes
  ha = HeatmapAnnotation(
    df = si %>% dplyr::select(
      #APOE, 
      Sex, 
      #DOB, 
      #Age, 
      Tube, 
      Isolation.Date, 
      #group
    ),
    col = list(
      #Group = group.colors,
      #Age = c(`3M` = 'white', `24M` = 'black')
    )
  )
  
  Heatmap(
    as.matrix(x),
    #t(scale(t(as.matrix(x)))),
    bottom_annotation = ha,
    column_split = factor(paste0(si$APOE, '.', si$Age), levels=c('E2.3M', 'E2.24M', 'E3.3M', 'E3.24M', 'E4.3M', 'E4.24M')),
    cluster_column_slices = FALSE,
    ...
    #column_split = paste0(si$APOE, '_', si$Age)
  )
}

