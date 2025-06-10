# This code implements differentially expressed protein analysis
# Code to load the dataset can be found in 'Analysis/proteomics/util-proteomics.r'

# Limma-trend analysis
source("Analysis/proteomics/util-proteomics.r") # load_pr, load_protein_sample_info
source("Analysis/scRNAseq/util-scRNAseq.r")
library(limma)
si = load_protein_sample_info()
sc <- load_roundtwo_integrated()
sc_mg <- load_roundtwo_integrated_celltype('MG')
#age_comparison <- data.table::fread("Results/DEP-Results/count_cutoff=0, prop_cutoff=0/limma-trend quantile/SEX+APOExAGE/MG/AGE/MG 24Mvs3M.txt")
#colnames(age_comparison)[1] = 'gene'

to_reclassify_as_mf <- colnames(sc_mg)[sc_mg$integrated_snn_res.0.3 %in% c(8,9,11)]
sc_mg$to_reclassify <- colnames(sc_mg) %in% to_reclassify_as_mf
#p <- DimPlot(sc_mg, group.by='to_reclassify') + 
#  labs(title = 'Cells to be reclassified',
#       subtitle = 'subclusters 8,9,11 we consider to not be MG - \nprobably some kind of macrophage')
#cowplot::save_plot(plot=p, filename='Results/DEP-Results/single-cell MG to reclassify.pdf', base_height=5, base_asp=1.2)
#DimPlot(sc, group.by='cell_type')

for (removing_8_9_11 in c(TRUE)) {
  # some clusters of microglia we later determined to be a kind of macrophage and not microglia.
  cell_type <- sc$cell_type
  if (removing_8_9_11) {
    cell_type[match(to_reclassify_as_mf, colnames(sc))] = 'Mf/Mg'
  }
  for (normalization.method in c('log Total.protein')) {
    # After discussion with Junmin Peng's lab, we've decided to apply "log" normalization.
    # We will not normalize by library size here, as that's already been done by Huan.
    for (celltype in c('AS', 'MG')) {
      counts <-  sc@assays$RNA@counts[, cell_type == celltype]
      for (prop_cutoff in c(0.001)) { # cutoff of 0.001 used for manuscript. We also considered 0.01 but felt that this was too high a cut-off.
        pr = load_pr()
        count_cutoff = 1
        if (removing_8_9_11) {
          outdir = file.path(
            "Results",
            "DEP-Results",
            'Removing 8,9,11 from microglia',
            paste0('count_cutoff=', count_cutoff, ', ',
                   'prop_cutoff=', prop_cutoff),
            get_outdir(normalization.method)
          )
        } else {
          outdir = file.path(
            "Results",
            "DEP-Results",
            'Without removing 8,9,11 from microglia',
            paste0('count_cutoff=', count_cutoff, ', ',
                   'prop_cutoff=', prop_cutoff),
            get_outdir(normalization.method)
          )
        }
        keep <- sparseMatrixStats::rowMeans2(counts > count_cutoff) > prop_cutoff
        proteins_keep <- intersect(rownames(pr[[celltype]]), rownames(counts)[keep])
        pr[[celltype]] <- pr[[celltype]][proteins_keep,]
        
        dat = get_normalized(celltype=celltype, pr=pr, normalization.method = normalization.method)
        write.csv(dat, paste0(outdir, '/', 'SEX+APOExAGE/', celltype, '/', celltype, ' total_protein_normalized_prop_cutoff=', prop_cutoff, '.csv'))
        si$Total.protein = colSums(assay(pr[[celltype]]))
        design = model.matrix(~0 + group + Tube + Sex, data=si)
        colnames(design) = stringr::str_remove(colnames(design), 'group|Tube')
        fit <- lmFit(dat, design=design)
         
        # ---- SEX + APOExAGE contrasts ----
        contrasts <- makeContrasts(
          `E2vsE3|24M` = (E2_24M - E3_24M),
          `E4vsE3|24M` = (E4_24M - E3_24M),
          `E2vsE4|24M` = (E2_24M - E4_24M),
          `E2vsE3|3M` = (E2_3M - E3_3M),
          `E4vsE3|3M` = (E4_3M - E3_3M),
          `E2vsE4|3M` = (E2_3M - E4_3M),
          `24Mvs3M|E2` = (E2_24M - E2_3M),
          `24Mvs3M|E3` = (E3_24M - E3_3M),
          `24Mvs3M|E4` = (E4_24M - E4_3M),
          `24Mvs3M.E2vsE3` = (E2_24M - E3_24M) - (E2_3M - E3_3M),
          `24Mvs3M.E4vsE3` = (E4_24M - E3_24M) - (E4_3M - E3_3M),
          #`24Mvs3M.E2vsE4` = (E2_24M - E4_24M) - (E2_3M - E4_3M),
          E2vsE3 = 0.5*(E2_24M + E2_3M) - 0.5*(E3_24M + E3_3M),
          E4vsE3 = 0.5*(E4_24M + E4_3M) - 0.5*(E3_24M + E3_3M),
          E2vsE4 = 0.5*(E2_24M + E2_3M) - 0.5*(E4_24M + E4_3M),
          `24Mvs3M` = 1/3*(E2_24M + E3_24M + E4_24M) - 1/3*(E2_3M + E3_3M + E4_3M),
          FvsM = SexF, #1/3*(X2021.07.21 + X2021.07.23 + X2021.07.24),
          levels=colnames(design)
        )
        dir.create(paste0(outdir, '/', 'SEX+APOExAGE/', celltype), recursive = TRUE)
        dir.create(paste0(outdir, '/', 'SEX+APOExAGE/', celltype, '/', 'APOE|AGE'))
        for (cont in 1:6) {# APOE | age
          fit2 <- contrasts.fit(fit, contrasts=contrasts[,cont])
          fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)
          write.table(topTable(fit2,number=2e4), file=paste0(outdir, '/', 'SEX+APOExAGE/', celltype, '/', 'APOE|AGE/', celltype, ' ', colnames(contrasts)[cont], '.txt'), sep='\t', col.names=T, row.names=T, quote=F)
        }
        dir.create(paste0(outdir, '/', 'SEX+APOExAGE/', celltype, '/', 'AGE|APOE'))
        for (cont in 7:9) {
          fit2 <- contrasts.fit(fit, contrast=contrasts[,cont])
          fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)
          write.table(topTable(fit2,number=2e4), file=paste0(outdir, '/', 'SEX+APOExAGE/', celltype, '/',  'AGE|APOE', '/', celltype, ' ', colnames(contrasts)[cont], '.txt'), sep='\t', col.names=T, row.names=T, quote=F)
        }
        dir.create(paste0(outdir, '/', 'SEX+APOExAGE/', celltype, '/', 'APOExAGE'))
        fit2 <- contrasts.fit(fit, contrast=contrasts[,10:11]) # APOE x age
        fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)
        write.table(topTable(fit2,number=2e4), file=paste0(outdir, '/', 'SEX+APOExAGE/', celltype, '/', 'APOExAGE/', celltype, ' APOExAGE.txt'), sep='\t', col.names=T, row.names=T, quote=F)
        dir.create(paste0(outdir, '/', 'SEX+APOExAGE/', celltype, '/', 'APOE'))
        for (cont in 12:14) {
          fit2 <- contrasts.fit(fit, contrast=contrasts[,cont])
          fit2 <- eBayes(fit2, trend=TRUE, robust=TRUE)
          write.table(topTable(fit2,number=2e4), file=paste0(outdir, '/', 'SEX+APOExAGE/', celltype, '/', 'APOE/', celltype, ' ', colnames(contrasts)[cont], '.txt'), sep='\t', col.names=T, row.names=T, quote=F)   
        }
        dir.create(paste0(outdir, '/', 'SEX+APOExAGE/', celltype, '/', 'AGE'))
        fit2 <- contrasts.fit(fit, contrast=contrasts[,15])
        fit2 <- eBayes(fit2, trend=TRUE, robust=TRUE)
        write.table(topTable(fit2,number=2e4), file=paste0(outdir, '/', 'SEX+APOExAGE/', celltype, '/', 'AGE/', celltype, ' ', colnames(contrasts)[15], '.txt'), sep='\t', col.names=T, row.names=T, quote=F)   
        dir.create(paste0(outdir, '/', 'SEX+APOExAGE/', celltype, '/', 'SEX'))
        fit2 <- contrasts.fit(fit, contrast=contrasts[,16])
        fit2 <- eBayes(fit2, trend=TRUE, robust=TRUE)
        write.table(topTable(fit2,number=2e4), file=paste0(outdir, '/', 'SEX+APOExAGE/', celltype, '/', 'SEX/', celltype, ' ', colnames(contrasts)[16], '.txt'), sep='\t', col.names=T, row.names=T, quote=F)
      }
    }
  }
}

#  ################################################################################
#  # AGE + APOE x SEX
#  ################################################################################
#  y$samples %<>% mutate(
#    group = paste0(outdir, '/', APOE, '_', Sex)
#  )
#  design = model.matrix(~0 + group + Age + Isolation.Batch, data=y$samples)
#  colnames(design) = stringr::str_remove(colnames(design), 'group|Isolation.Batch')
#  logCPM = cpm(y, log=TRUE)
#  fit = lmFit(logCPM, design=design)
#   contrasts <- makeContrasts(
#     `E2vsE3|M` = (E2_M - E3_M),
#     `E4vsE3|M` = (E4_M - E3_M),
#     `E2vsE4|M` = (E2_M - E4_M),
#     `E2vsE3|F` = (E2_F - E3_F),
#     `E4vsE3|F` = (E4_F - E3_F),
#     `E2vsE4|F` = (E2_F - E4_F),
#     `FvsM|E2` = (E2_F - E2_M),
#     `FvsM|E3` = (E3_F - E3_M),
#     `FvsM|E4` = (E4_F - E4_M),
#     `FvsM.E2vsE3` = (E2_F - E3_F) - (E2_M - E3_M),
#     `FvsM.E4vsE3` = (E4_F - E3_F) - (E4_M - E3_M),
#     E2vsE3 = 0.5*(E2_F + E2_M) - 0.5*(E3_F + E3_M),
#     E4vsE3 = 0.5*(E4_F + E4_M) - 0.5*(E3_F + E3_M),
#     E2vsE4 = 0.5*(E2_F + E2_M) - 0.5*(E4_F + E4_M),
#     `FvsM` = 1/3*(E2_F + E3_F + E4_F) - 1/3*(E2_M + E3_M + E4_M),
#     `24Mvs3M` = Age24M, 
#     levels=colnames(design)
#  )
#  dir.create(paste0(outdir, '/', 'AGE+APOExSEX/', celltype), recursive = TRUE)
#  dir.create('APOE|SEX')
#  for (cont in 1:6) {# APOE | age
#    fit2 <- contrasts.fit(fit, contrast=contrasts[,cont])
#    fit2 <- eBayes(fit2, trend=TRUE, robust=TRUE)
#    write.table(topTable(fit2,number=2e4), file=paste0(outdir, '/', 'APOE|SEX/', celltype, ' ', colnames(contrasts)[cont], '.txt'), sep='\t', col.names=T, row.names=T, quote=F)
#  }
#  dir.create('SEX|APOE')
#  for (cont in 7:9) {
#    fit2 <- contrasts.fit(fit, contrast=contrasts[,cont])
#    fit2 <- eBayes(fit2, trend=TRUE, robust=TRUE)
#    write.table(topTable(fit2,number=2e4), file=paste0(outdir, '/', 'SEX|APOE/', celltype, ' ', colnames(contrasts)[cont], '.txt'), sep='\t', col.names=T, row.names=T, quote=F)
#  }
#  dir.create('APOExSEX')
#  fit2 <- contrasts.fit(fit, contrast=contrasts[,10:11]) # APOE x age
#  fit2 <- eBayes(fit2, trend=TRUE, robust=TRUE)  
#  write.table(topTable(fit2,number=2e4), file=paste0(outdir, '/', 'APOExSEX/', celltype, ' APOExSEX.txt'), sep='\t', col.names=T, row.names=T, quote=F)
#  dir.create('APOE')
#  for (cont in 12:14) {
#    fit2 <- contrasts.fit(fit, contrast=contrasts[,cont])
#    fit2 <- eBayes(fit2, trend=TRUE, robust=TRUE)  
#    write.table(topTable(fit2,number=2e4), file=paste0(outdir, '/', 'APOE/', celltype, ' ', colnames(contrasts)[cont], '.txt'), sep='\t', col.names=T, row.names=T, quote=F)   
#  }
#  dir.create('SEX')
#  fit2 <- contrasts.fit(fit, contrast=contrasts[,15])
#  fit2 <- eBayes(fit2, trend=TRUE, robust=TRUE)  
#  write.table(topTable(fit2,number=2e4), file=paste0(outdir, '/', 'SEX/', celltype, ' ', colnames(contrasts)[15], '.txt'), sep='\t', col.names=T, row.names=T, quote=F)   
#  dir.create('AGE')
#  fit2 <- contrasts.fit(fit, contrast=contrasts[,16])
#  fit2 <- eBayes(fit2, trend=TRUE, robust=TRUE)  
#  write.table(topTable(fit2,number=2e4), file=paste0(outdir, '/', 'AGE/', celltype, ' ', colnames(contrasts)[16], '.txt'), sep='\t', col.names=T, row.names=T, quote=F)
#  
#  ################################################################################
#  # APOE + AGE x SEX
#  ################################################################################
#  y$samples %<>% mutate(
#    group = paste0(outdir, '/', Sex, '_', Age)
#  )
#  design = model.matrix(~0 + group + APOE + Isolation.Batch, data=y$samples)
#  colnames(design) = stringr::str_remove(colnames(design), 'group|Isolation.Batch')
#  logCPM = cpm(y, log=TRUE)
#  fit = lmFit(logCPM, design=design)
#  contrasts <- makeContrasts(
#    `24Mvs3M|M` = (`M_24M` - `M_3M`),
#    `24Mvs3M|F` = (`F_24M` - `F_3M`),
#    `FvsM|24M` = (`F_24M` - `M_24M`),
#    `FvsM|3M` = (`F_3M` - `M_3M`),
#    `24Mvs3M.FvsM` = (`M_24M` - `M_3M`) - (`F_24M` - `F_3M`),
#    E2vsE3 = APOEE2,
#    E4vsE3 = APOEE4,
#    E2vsE4 = APOEE2 - APOEE4,
#    `FvsM` = 1/2*(`F_24M` + `F_3M`) - 1/2*(`M_24M` + `M_3M`),
#    `24Mvs3M` = 1/2*(`F_24M` + `M_24M`) - 1/2*(`F_3M` + `M_3M`),
#    levels=colnames(design)
#  )
#  dir.create(paste0(outdir, '/', 'APOE+AGExSEX/', celltype), recursive = TRUE)
#  dir.create('AGE|SEX')
#  for (cont in 1:2) {# APOE | age
#    fit2 <- contrasts.fit(fit, contrast=contrasts[,cont])
#    fit2 <- eBayes(fit2, robust=TRUE, trend=TRUE)
#    write.table(topTable(fit2,number=2e4), file=paste0(outdir, '/', 'AGE|SEX/', celltype, ' ', colnames(contrasts)[cont], '.txt'), sep='\t', col.names=T, row.names=T, quote=F)
#  }
#  dir.create('SEX|AGE')
#  for (cont in 3:4) {
#    fit2 <- contrasts.fit(fit, contrast=contrasts[,cont])
#    fit2 <- eBayes(fit2, robust=TRUE, trend=TRUE)
#    write.table(topTable(fit2,number=2e4), file=paste0(outdir, '/', 'SEX|AGE/', celltype, ' ', colnames(contrasts)[cont], '.txt'), sep='\t', col.names=T, row.names=T, quote=F)
#  }
#  dir.create('AGExSEX')
#  fit2 <- contrasts.fit(fit, contrast=contrasts[,5])
#  fit2 <- eBayes(fit2, robust=TRUE, trend=TRUE)
#  write.table(topTable(fit2,number=2e4), file=paste0(outdir, '/', 'AGExSEX/', celltype, ' AGExSEX.txt'), sep='\t', col.names=T, row.names=T, quote=F)
#  dir.create('APOE')
#  for (cont in 6:8) {
#    fit2 <- contrasts.fit(fit, contrast=contrasts[,cont])
#    fit2 <- eBayes(fit2, robust=TRUE, trend=TRUE)
#    write.table(topTable(fit2,number=2e4), file=paste0(outdir, '/', 'APOE/', celltype, ' ', colnames(contrasts)[cont], '.txt'), sep='\t', col.names=T, row.names=T, quote=F)   
#  }
#  dir.create('SEX')
#  fit2 <- contrasts.fit(fit, contrast=contrasts[,9])
#  fit2 <- eBayes(fit2, robust=TRUE, trend=TRUE)
#  write.table(topTable(fit2,number=2e4), file=paste0(outdir, '/', 'SEX/', celltype, ' ', colnames(contrasts)[9], '.txt'), sep='\t', col.names=T, row.names=T, quote=F)   
#  dir.create('AGE')
#  fit2 <- contrasts.fit(fit, contrast=contrasts[,10])
#  fit2 <- eBayes(fit2, robust=TRUE, trend=TRUE)
#  write.table(topTable(fit2,number=2e4), file=paste0(outdir, '/', 'AGE/', celltype, ' ', colnames(contrasts)[10], '.txt'), sep='\t', col.names=T, row.names=T, quote=F)
#}
##  contrasts <- makeContrasts(
##    `E2vsE3|24M` = (E2_24M - E3_24M),
##    `E4vsE3|24M` = (E4_24M - E3_24M),
##    `E2vsE4|24M` = (E2_24M - E4_24M),
##    `E2vsE3|3M` = (E2_3M - E3_3M),
##    `E4vsE3|3M` = (E4_3M - E3_3M),
##    `E2vsE4|3M` = (E2_3M - E4_3M),
##    `24Mvs3M|E2` = (E2_24M - E2_3M),
##    `24Mvs3M|E3` = (E3_24M - E3_3M),
##    `24Mvs3M|E4` = (E4_24M - E4_3M),
##    `24Mvs3M.E2vsE3` = (E2_24M - E3_24M) - (E2_3M - E3_3M),
##    `24Mvs3M.E4vsE3` = (E4_24M - E3_24M) - (E4_3M - E3_3M),
##    E2vsE3 = 0.5*(E2_24M + E2_3M) - 0.5*(E3_24M + E3_3M),
##    E4vsE3 = 0.5*(E4_24M + E4_3M) - 0.5*(E3_24M + E3_3M),
##    E2vsE4 = 0.5*(E2_24M + E2_3M) - 0.5*(E4_24M + E4_3M),
##    `24Mvs3M` = 1/3*(E2_24M + E3_24M + E4_24M) - 1/3*(E2_3M + E3_3M + E4_3M),
##    FvsM = SexF, #1/3*(X2021.07.21 + X2021.07.23 + X2021.07.24),
##    levels=colnames(design)
##  ) 
##  y.list[[celltype]] = y
##  fit.list[[celltype]] = fit
##  result.list[[celltype]][['24M']] = glmQLFTest(fit.list[[celltype]], contrast=contrasts[,c(1,2)])
##  result.list[[celltype]][['3M']]  = glmQLFTest(fit.list[[celltype]], contrast=contrasts[,c(4,5)])
##  result.list[[celltype]][['APOE']]  = glmQLFTest(fit.list[[celltype]], contrast=contrasts[,c(1,2,4,5)])
##saveRDS(y.list, 'Datasets/protein-DGEList.rds')
##saveRDS(fit.list, 'Datasets/protein-QLF-fit-list.rds')
#
