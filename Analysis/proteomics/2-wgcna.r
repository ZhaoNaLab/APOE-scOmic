# Complete to perform WGCNA analysis.
# Code to load proteomics dataset is in 'Analysis/proteomics/util-proteomics.r'

source('Analysis/util/wgcna-wrappers.r') # blockwise.wrapper, moduleTraitCor, plotDendroAndColors.wrapper, plotModuleTraitCorrelation.wrapper
source('Analysis/plot/heatmap.r') # draw_heatmap
source("Analysis/proteomics/util-proteomics.r")  # load_protein_sample_info, load_pr
source("Analysis/scRNAseq/util-scRNAseq.r") # load_roundtwo_integrated, load_roundtwo_integrated_celltype
options(stringsAsFactors = FALSE)
library(batchtools)
library(WGCNA)
library(ComplexHeatmap)
library(clusterProfiler)
library(org.Mm.eg.db)

si <- load_protein_sample_info()
sc <- load_roundtwo_integrated()
x_mg <- load_roundtwo_integrated_celltype(celltype='MG')
to_reclassify_as_mf <- colnames(x_mg)[x_mg$integrated_snn_res.0.3 %in% c(8,9,11)]
x_mg$to_reclassify <- colnames(x_mg) %in% to_reclassify_as_mf
pr <- load_pr()
#write.csv(cbind(as.data.frame(rowData(pr$AS)), assay(pr$AS)), 'Datasets/AS full table.csv')
#write.csv(cbind(as.data.frame(rowData(pr$MG)), assay(pr$MG)), 'Datasets/MG full table.csv')

# This function helps prepare proteomics data by filtering proteome to just the ones that are abundant in the scRNAseq for each celltype
filtered_problem <- function(data, to_reclassify=FALSE, latent.vars = NULL, normalization.method=FALSE, prop_cutoff, count_cutoff) {
  library(Seurat)
  source("Analysis/proteomics/util-proteomics.r")
  #sc <- load_roundtwo_integrated()
  cell_type = sc$cell_type
  if (to_reclassify) {
    cell_type[match(to_reclassify_as_mf, colnames(sc))] = 'Mf/Mg'
  }
  
  counts <- sc@assays$RNA@counts[, cell_type == celltype]
  keep <- sparseMatrixStats::rowMeans2(counts > count_cutoff) > prop_cutoff
  proteins_keep <- intersect(rownames(data), rownames(counts)[keep])
  
  data <- data[proteins_keep,]
  dat = assay(data)
  if (normalization.method == 'log') {
    dat = log(dat)
  } else if (normalization.method == 'median') {
    dat <-  log(dat / matrixStats::colMedians(as.matrix(dat)))
  } else {
    dat = limma::normalizeBetweenArrays(log(dat), method = normalization.method)
  }
  
  if (!is.null(latent.vars)) {
    fit = lm(as.formula(paste('t(dat) ~', paste(latent.vars, collapse='+'))), data=colData(data))
    dat = t(fit$residuals)
  }
  x = SummarizedExperiment(assays=list(data=dat),
                           colData=colData(data),
                           rowData=rowData(data))
  return(x)
}

# ----------- abundance filtering of the proteins based on scRNAseq mRNA abundances -------
for(celltype in c('MG', 'AS')) {
  normalization.method = 'log'
  prop_cutoff = 0.001
  count_cutoff = 1
  latent.vars <- c('Tube', 'Sex', 'Total.protein')
  outdir = file.path(
    "Results",
    "WGCNA-Results",
    'Removing 8,9,11 from microglia',
    paste0('count_cutoff=', count_cutoff, ', ',
           'prop_cutoff=', prop_cutoff),
    paste0(get_outdir(normalization.method), ' latent.vars=', paste0(latent.vars, collapse=',')),
    celltype
  )
  dir.create(outdir, recursive = TRUE)
  # run through the problem definition to obtain below x, 6562 by 36
  x <- filtered_problem(data=pr[[celltype]], to_reclassify = TRUE, latent.vars=latent.vars, normalization.method = normalization.method, count_cutoff=count_cutoff, prop_cutoff=prop_cutoff) #6562 36 (MG); 7212 36 (AS)
  write.csv(assay(x), file.path(outdir, paste0(celltype, ' ', paste0(latent.vars, collapse=','), ' residuals.csv')))
  
  si <- load_protein_sample_info()
  si$Total.protein = colData(pr[[celltype]])$Total.protein
  # test with MG problem
  disableWGCNAThreads()
  sft <- WGCNA::pickSoftThreshold(t(assay(x)), dataIsExpr = TRUE, corFnc = WGCNA::bicor, networkType = 'signed', verbose = 5, powerVector = 1:30, RsquaredCut = 0.80, blockSize = nrow(x))
  
  plot_sft(sft)
  
  minModuleSize = 20
  mergeCutHeight = 0.25
  net <- blockwise.wrapper(
    data = NA, job=NA, instance=x, 
    beta=sft$powerEstimate, # 23 - MG, 18 - AS (prop.cutoff=0.001)
    networkType = 'signed', # a _signed_ network
    loadTOM = TRUE,
    minModuleSize = minModuleSize,
    mergeCutHeight = mergeCutHeight,
    saveTOMFileBase = file.path(outdir, 'allSamplesTOM')
  )
  net.name <- paste0(celltype, ' minModuleSize=', minModuleSize, 'mergeCutHeight=', mergeCutHeight)
  plotDendroAndColors.wrapper(net=net, path=outdir, prefix=paste('dendrogram', net.name))
  design <- model.matrix(as.formula(paste0('~0 + group +', paste0(latent.vars, collapse='+'))), data=si)
  myContrasts <- limma::makeContrasts(
    `Old.vs.Young|E2` = groupE2_24M - groupE2_3M,
    `Old.vs.Young|E3` = groupE3_24M - groupE3_3M,
    `Old.vs.Young|E4` = groupE4_24M - groupE4_3M,
    levels=colnames(design)
  )
  
  datTraits = as.data.frame(colData(pr[[celltype]]))
  datTraits <- datTraits[, c('Tube', 'Sex')]
  datTraits.binary <- binarizeCategoricalColumns.pairwise(
    datTraits,
    levelOrder = c(#Genotype = c('E2', 'E3', 'E4'),
                   #Group = c('Young', 'aged'),
                   Tube = c('Regular', 'LB'),
                   Sex = c('M', 'F')))
  datTraits.binary[[paste0(celltype, '.total.protein')]] = log(colSums(assay(pr[[celltype]])))
  datTraits.binary <- cbind(datTraits.binary, as.data.frame(design %*% myContrasts))
  # 05/29/2024 -- thanks to Qiao we found a mistake with using this function. We need to subset the above contrasts for E2, E3, E4 to E2, E3, E4, respectively.
  datExpr <- t(assay(x))
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  moduleColors = labels2colors(net$colors)
  name <- paste(celltype, net.name)
  path <- outdir
  MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs = orderMEs(MEs0)
  write.table(MEs, paste0(path, '/', 'MEs - ', net.name, '.txt'), sep='\t', row.names=T, col.names=T, quote=F)
  moduleTraitCor = cbind(
    cor(MEs, as.matrix(datTraits.binary[, c('Tube.LB.vs.Regular', 'Sex.M.vs.F', paste0(celltype, '.total.protein'))]), use="p", method=cor_method),
    data.frame(
      Old.vs.Young.E2= cor(MEs[si$APOE == 'E2',], (datTraits.binary[si$APOE == 'E2', 'Old.vs.Young|E2']), use="p", method=cor_method),
      Old.vs.Young.E3= cor(MEs[si$APOE == 'E3',], (datTraits.binary[si$APOE == 'E3', 'Old.vs.Young|E3']), use="p", method=cor_method),
      Old.vs.Young.E4= cor(MEs[si$APOE == 'E4',], (datTraits.binary[si$APOE == 'E4', 'Old.vs.Young|E4']), use="p", method=cor_method)
  )) %>%
    as(Class=c("matrix", "array"))
  
  moduleTraitPvalue <- cbind(
    corPvalueStudent(as.matrix(moduleTraitCor[, c('Tube.LB.vs.Regular', 'Sex.M.vs.F', paste0(celltype, '.total.protein'))]), nSamples=nSamples),
    data.frame(
      Old.vs.Young.E2 = corPvalueStudent(moduleTraitCor[, 'Old.vs.Young.E2'], nSamples=nrow(datTraits.binary[si$APOE == 'E2',])),
      Old.vs.Young.E3 = corPvalueStudent(moduleTraitCor[, 'Old.vs.Young.E3'], nSamples=nrow(datTraits.binary[si$APOE == 'E3',])),
      Old.vs.Young.E4 = corPvalueStudent(moduleTraitCor[, 'Old.vs.Young.E4'], nSamples=nrow(datTraits.binary[si$APOE == 'E4',]))
    )
  )  %>%
    as(Class=c("matrix", "array"))
  pdf(paste0(path, '/', 'MTC - ', net.name, '.pdf'), wi=9, he=11+0.1*(ncol(net$MEs)-20)*(ncol(net$MEs) > 20))
  textMatrix=paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep="")
  dim(textMatrix) = dim(moduleTraitCor)
  par(mar=c(6,8.5+0.1*(ncol(net$MEs)-20)*(ncol(net$MEs) > 20),3,3))
  #  par(mar=c(1,1,1,1))
  labeledHeatmap(Matrix=moduleTraitCor,
                 xLabels=names(datTraits.binary),
                 yLabels=names(MEs),
                 ySymbols=names(MEs),
                 colorLabels=FALSE,
                 colors=blueWhiteRed(50),
                 textMatrix=textMatrix,
                 setStdMargins=FALSE,
                 cex.text=0.5,
                 zlim=c(-1,1),
                 main=paste("Module-trait relationships -", net.name))
  dev.off()
  
  #plotModuleTraitCorrelation.wrapper(
  #  t(assay(x)),
  #  datTraits = datTraits.binary,
  #  net=net,
  #  path=outdir,
  #  cor_method = 'pearson',
  #  name=paste(celltype, net.name)
  #)
  
  MEs <- read.table(file.path(outdir, list.files(outdir, pattern='MEs.*.txt')))
  p <- cbind(MEs, si) %>% 
    mutate(APOE = relevel(APOE, 'E2')) %>% 
    tidyr::pivot_longer(starts_with('ME'), names_to='module', values_to='ME') %>% 
   # filter(module == 'MEred') %>%
    ggplot(aes(x = interaction(APOE, Age, lex.order = TRUE), y=ME)) +
    facet_wrap(~module)
  module_colors <- unique(labels2colors(net$colors))
  names(module_colors) <- paste0('ME', module_colors)
  cowplot::save_plot(
    plot = p + 
    geom_boxplot(mapping = aes(fill=module),
                 position = position_dodge(),
                 outlier.size = 0) +
    geom_jitter(width = 0.1,
                fill='white', 
                shape=21) + 
    scale_fill_manual(values=module_colors) + 
    theme(axis.text.x = element_text(angle=45, hjust=1)) + 
    labs(
      y = 'Module Eigengene',
      x = 'APOE Genotype x Age' 
    ),
    filename = file.path(outdir, 'ME boxplots.pdf'),
    base_width = 10,
    base_height=8
  )
  
  # ----- Analysis of variance -----
  MEs <- read.table(file.path(outdir, list.files(outdir, pattern='MEs.*.txt')))
  fit = lm(as.matrix(MEs) ~ APOE + Age + APOE:Age, data=si)
  aov_fit <- aov(fit)
  anova(aov_fit)
  writeLines(
    capture.output(summary(aov_fit)),
    file.path(outdir, 'aov_summary.txt')
  )
  
  for(color in module_colors) {
    genes=rownames(x)[labels2colors(net$colors) == color]
    xh <- assay(x)[genes,]
    rownames(xh) <- genes
    ha = HeatmapAnnotation(
      df = si %>% dplyr::select(
        Sex, 
        Tube, 
        Isolation.Date, 
      ),
      col = list(
        #Group = group.colors,
        #Age = c(`3M` = 'white', `24M` = 'black')
      )
    )
    
    h <- Heatmap(
          as.matrix(xh),
          #t(scale(t(as.matrix(x)))),
          bottom_annotation = ha,
          column_split = factor(paste0(si$APOE, '.', si$Age), levels=c('E2.3M', 'E2.24M', 'E3.3M', 'E3.24M', 'E4.3M', 'E4.24M')),
          cluster_column_slices = FALSE,
          row_names_gp = gpar(fontsize=5)
          #column_split = paste0(si$APOE, '_', si$Age)
        )
    cowplot::save_plot(
      plot = draw_heatmap(h),
      filename = file.path(outdir, paste0(color, '.pdf')),
      base_height= min(nrow(xh)*0.05 + 4, 10),
      base_width = min(ncol(xh)*0.05*5 + 4, 8)
      #base_asp=1.25
    )
  }
  
  module_enrichment_results <- list()
  for(color in module_colors) {
    genes=rownames(x)[labels2colors(net$colors) == color]
    module_enrichment_results[[color]] <- enrichGO(
      gene = genes,
      OrgDb = org.Mm.eg.db, 
      keyType = 'SYMBOL', 
      ont = 'ALL', 
      pvalueCutoff = 1,
      pAdjustMethod = 'bonferroni', 
      universe = rownames(x), 
      minGSSize = 10
    ) 
  }
  purrr::map(module_enrichment_results, ~.x@result) %>%
    writexl::write_xlsx(path=file.path(outdir, paste0(celltype, ' module-GO-over-representation.xlsx')))
  
  obtain_module_membership(t(assay(x)), datTraits.binary, MEs, labels2colors(net$colors), outdir, celltype)
  
  GO_DATA <- clusterProfiler:::get_GO_data(org.Mm.eg.db, 'ALL', 'SYMBOL')
  geneSets = DOSE:::getGeneSet(GO_DATA)
  pValueCutoff = 1
  pAdjustMethod = 'bonferroni'
  
  module_gsea_results <- list()
  geneInfo <- read.csv(file.path(outdir, paste0('geneInfo - ', celltype, '.csv')))
  for(color in module_colors[4:length(module_colors)]) {
    geneList <- geneInfo[[paste0('MM.', color)]][geneInfo$moduleColor == color]
    names(geneList) <- geneInfo$geneSymbol[geneInfo$moduleColor == color]
    geneList <- sort(geneList, decreasing = TRUE)
    
    tmp_res <- fgsea::fgsea(
      stats = geneList,
      pathways = geneSets, 
      minSize = 10,
      maxSize = min(500, length(geneList) - 1),
      eps = 1e-10, 
      gseaParam = 1,
      nproc = 1
    )
    p.adj <- p.adjust(tmp_res$pval, method=pAdjustMethod)
    qvalues <- DOSE:::calculate_qvalue(tmp_res$pval)
    Description <- DOSE:::TERM2NAME(tmp_res$pathway, GO_DATA)
    res <- data.frame(
      ID = as.character(tmp_res$pathway),
      Description = unname(Description),
      setSize = tmp_res$size,
      enrichmentScore = tmp_res$ES,
      NES = tmp_res$NES,
      pvalue = tmp_res$pval,
      p.adjust = p.adj,
      qvalue = qvalues,
      stringsAsFactors = FALSE
    )
    
    res <- res[!is.na(res$pvalue),]
    res <- res[ res$pvalue <= pValueCutoff, ]
    res <- res[ res$p.adjust <= pValueCutoff, ]
    idx <- order(res$p.adjust, -abs(res$NES), decreasing = FALSE)
    res <- res[idx, ]
    module_gsea_results[[color]] <- res
  }
  
  writexl::write_xlsx(module_gsea_results, path=file.path(outdir, paste0(celltype, ' module-GO-GSEA.xlsx')))
}
