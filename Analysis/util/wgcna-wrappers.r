compute_sft <- function(x, ...) {
  powers <- c(c(1:10), seq(from=12, to=20, by=2))
  WGCNA::pickSoftThreshold(SummarizedExperiment::assay(x), powerVector=powers, verbose=5, ...)
}

plot_sft <- function(sft, filename = 'Network-Topology-Analysis.png') {
  library(magrittr)
  library(patchwork)
  p1 <- sft$fitIndices %>% 
    ggplot(aes(x = Power, y = -sign(slope)*SFT.R.sq, label = Power)) +
    geom_text() +
    geom_abline(aes(color = 'red', slope = 0, intercept = 0.85)) +
    ylim(c(-0.75, 1)) + 
    xlab('Soft Threshold (power)') + 
    ylab('Scale Free Topology Model Fit, signed R^2') +
    labs(title = paste('Scale independence. Beta =', sft$powerEstimate)) +
    theme_minimal() +
    theme(legend.position = 'none')

  p2 <- sft$fitIndices %>% 
    ggplot(aes(x = Power, y =mean.k., label = Power)) +
    geom_text() +
    xlab('Soft Threshold (power)') + 
    ylab('Mean Connectivity') +
    labs(title = 'Mean Connectivity') +
    theme_minimal()

  if (!is.null(filename)) {
    cowplot::save_plot(plot = p1 + p2, filename = filename, base_asp=2)
  }
  return(p1 + p2)
}

pvalue_code <- function(p.values, corr = FALSE, na = FALSE, cutpoints = c(0,  0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "), ...) {
  symnum(p.values, corr = corr, na = na, cutpoints = cutpoints, symbols = symbols, ...)
}

## wrapper for plot_sft function
sft.wrapper <- function(data, job, instance, path, ...) {
  powers <- c(c(1:10), seq(from=12, to=20, by=2))
  sft <- WGCNA::pickSoftThreshold(SummarizedExperiment::assay(data), powerVector=powers, verbose=5, ...)
  plot_sft(sft, filename = paste0(path, '/', 'Network-Topology-Analysis.png'))
}

## wrapper for blockwiseModules function
# data is samples x genes
blockwise.wrapper <- function(data, job, instance, future.globals.maxSize = 20*1024^3, beta, minModuleSize, mergeCutHeight, deepSplit=2, TOMType = 'signed',
                              networkType='signed',
                              reassignThreshold=0, numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = TRUE, saveTOMFileBase = 'allSamplesTOM',
                              maxBlockSize = NULL, corType = 'bicor', nThreads = 16, verbose = 3, ...) {
  #library(future)
  options(stringsAsFactors = FALSE)
  #plan(multicore, workers=nThreads)
  #options(future.globals.maxSize=future.globals.maxSize)
  enableWGCNAThreads()

  if (class(instance)[1] == 'SummarizedExperiment') { datExpr = t(SummarizedExperiment::assay(instance)) }
  if (class(instance)[1] %in% c('matrix', 'array', 'data.frame')) { datExpr = instance }

  if (is.null(maxBlockSize)) { maxBlockSize = ncol(datExpr) }

  net = blockwiseModules(datExpr, power = beta,
                         TOMType = TOMType, minModuleSize = minModuleSize,
                         reassignThreshold = reassignThreshold, mergeCutHeight = mergeCutHeight,
                         numericLabels = numericLabels, pamRespectsDendro = pamRespectsDendro,
                         saveTOMs = saveTOMs,
                         saveTOMFileBase = saveTOMFileBase,
                         maxBlockSize = maxBlockSize,
                         corType = corType,
                         nThreads = nThreads,
                         deepSplit=deepSplit,
                         verbose = verbose,
                         ...)
  return(net)
}

## wrapper for blockwiseConsensusModules function
# data is a multi-set 
blockwiseConsensus.wrapper <- function(data, job, instance, future.globals.maxSize = 20*1024^3, beta, minModuleSize, mergeCutHeight, deepSplit=2, TOMType = 'signed', 
                                       networkType='signed',
                                       reassignThreshold=0, numericLabels = TRUE, pamRespectsDendro = FALSE, saveIndividualTOMs=FALSE, maxBlockSize = max(unlist(lapply(instance, function(x) ncol(x$data)))), corType = 'bicor', nThreads = 8, verbose = 3, ...) {
  library(future)
  options(stringsAsFactors = FALSE)
  plan(multicore, workers=nThreads)
  options(future.globals.maxSize=future.globals.maxSize)
  enableWGCNAThreads()

  net = blockwiseConsensusModules(instance, power = beta,
                                  TOMType = TOMType, minModuleSize = minModuleSize,
                                  reassignThreshold = reassignThreshold, mergeCutHeight = mergeCutHeight,
                                  numericLabels = numericLabels, pamRespectsDendro = pamRespectsDendro,
                                  saveIndividualTOMs=saveIndividualTOMs,
                                  maxBlockSize = maxBlockSize,
                                  corType = corType,
                                  nThreads = nThreads,
                                  deepSplit=deepSplit,
                                  verbose = verbose,
                                  ...)
  return(net)
}


## wrapper for plotDendroAndColors
plotDendroAndColors.wrapper <- function(net, path=NULL, prefix='Dendrogram', dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05, ...) {
  if (!is.null(path)) { pdf(paste0(path, '/', prefix, '.pdf'), wi=8, he=6) }
  mergedColors <- labels2colors(net$colors)
  WGCNA::plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],"Module colors", dendroLabels = dendroLabels, hang = hang, addGuide = addGuide, guideHang=guideHang, ...)
  if (!is.null(path)) { dev.off() }
}

## Wrapper for module-trait correlations.
# net  : network
# meta : data frame of traits you want to correlate with. Ensure that rownames can re-arrange net$MEs
#      into the correct order.
plotModuleTraitCorrelation.wrapper <- function(datExpr, datTraits, net, name, path='./', cor_method = c("pearson", "kendall", "spearman"), ...) {
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  moduleColors = labels2colors(net$colors)
  MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs = orderMEs(MEs0)
  write.table(MEs, paste0(path, '/', 'MEs - ', name, '.txt'), sep='\t', row.names=T, col.names=T, quote=F)
  moduleTraitCor = cor(MEs, as.matrix(datTraits), use="p", method=cor_method)
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
  pdf(paste0(path, '/', 'MTC - ', name, '.pdf'), wi=9, he=11+0.1*(ncol(net$MEs)-20)*(ncol(net$MEs) > 20))
  textMatrix=paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep="")
  dim(textMatrix) = dim(moduleTraitCor)
  par(mar=c(6,8.5+0.1*(ncol(net$MEs)-20)*(ncol(net$MEs) > 20),3,3))
#  par(mar=c(1,1,1,1))
  labeledHeatmap(Matrix=moduleTraitCor,
                 xLabels=names(datTraits),
                 yLabels=names(MEs),
                 ySymbols=names(MEs),
                 colorLabels=FALSE,
                 colors=blueWhiteRed(50),
                 textMatrix=textMatrix,
                 setStdMargins=FALSE,
                 cex.text=0.5,
                 zlim=c(-1,1),
                 main=paste("Module-trait relationships -", name))
  dev.off()
}

## Wrapper for consensus module-trait correlations.
# net  : network
# meta : data frame of traits you want to correlate with. Ensure that rownames can re-arrange net$MEs
#      into the correct order.
plotConsensusModuleTraitCorrelation.wrapper <- function(net, meta, add=0.5, ...) {
  library(ggplot2)
  library(ggdendro)
  library(patchwork)
  for (i in 1:length(net$multiMEs)) {
    net$MEs <- net$multiMEs[[i]]$data
    meta <- meta[[1]]$data
    net$MEs <- net$MEs[rownames(meta),]
    moduleTraitCor <- WGCNA::cor(net$MEs, meta)
    moduleTraitPValues <- WGCNA::corPvalueStudent(moduleTraitCor, nrow(meta))
    moduleTraits <- left_join(moduleTraitCor %>% reshape2::melt() %>% rename(Module = Var1, Trait = Var2, Cor = value),
                              moduleTraitPValues %>% reshape2::melt() %>% rename(Module = Var1, Trait = Var2, p.value = value), by = c('Module', 'Trait'))
    moduleTraits %<>% mutate(Signif = pvalue_code(p.value, legend=FALSE))
    model <- hclust(dist(t(net$MEs)), 'ave')
    dhc <- as.dendrogram(model)
    ddata <- dendro_data(dhc, type = 'rectangle')
    p <- ggplot() + 
    geom_segment(data = ddata$segments, aes(x = x, y = y, xend = xend, yend=yend)) +
    theme_dendro() + 
    scale_x_continuous(expand = expansion(add = add)) #, 1 / ( 4 * nrow(moduleTraitCor) + 0.01) )))
    p2 <- moduleTraits %>%
    ggplot(aes(x = factor(Module, ddata$labels$label), y = Trait, fill = Cor, label = Signif)) +
    geom_tile(show.legend=FALSE) +
    geom_text(size=4) + 
    scale_fill_gradient2(low='blue', mid='white', high='red', limits=c(-1,1)) + 
    theme_minimal() +
    xlab(net$name) + 
    ylab(NULL) + 
    theme(axis.text.x = element_text(angle = 90))
    p <- p / p2
#  pdf(paste0(path, '/', prefix, '.pdf'), width=15, height=10)
    print(p)
#  dev.off()
  }
}

obtain_module_membership <- function(datExpr, datTraits, MEs, moduleColors, path, name) {
  nSamples = nrow(datExpr)
  modNames = substring(names(MEs), 3)
  geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
  names(geneModuleMembership) = paste("MM", modNames, sep="")
  names(MMPvalue) = paste("p.MM", modNames, sep="")
  geneTraitSignificance = as.data.frame(sapply(datTraits, FUN = function(x) cor(datExpr, x, use="p")))
  rownames(geneTraitSignificance) <- colnames(datExpr)
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
  names(geneTraitSignificance) = paste("GS.", names(datTraits), sep="")
  names(GSPvalue) = paste("p.GS.", names(datTraits), sep="")

  geneInfo0 = data.frame(geneSymbol = colnames(datExpr),
                         moduleColor = moduleColors,
                         geneTraitSignificance,
                         GSPvalue)
  # Order modules by their significance for the first trait of datTraits
  modOrder = order(-abs(cor(MEs, datTraits[,1], use='p')))
  # Add module membership information in the chosen order
  for (mod in 1:ncol(geneModuleMembership))
  {
    oldNames = names(geneInfo0)
    geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                           MMPvalue[, modOrder[mod]]);
    names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                         paste("p.MM.", modNames[modOrder[mod]], sep=""))
  }
  geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0[[paste0('GS.', names(datTraits)[1])]]));
  geneInfo = geneInfo0[geneOrder, ]
  write.csv(geneInfo, file = paste0(path, '/', "geneInfo - ", name, ".csv"), row.names=F, quote=F)
}

# Another module trait correlation (this type using WGCNA moduletrait correlation)
moduleTraitCorrelation.wrapper <- function(datExpr, datTraits, net, name, path='./', biomart="ensembl", host="http://www.ensembl.org", dataset="mmusculus_gene_ensembl", attributes=c("chromosome_name", "start_position", "end_position", "external_gene_name", "entrezgene_id", "ensembl_gene_id"), attribute_given = 'ensembl_gene_id', organism="mouse", GOcollection=NULL) {
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  moduleColors = labels2colors(net$colors)
  MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs = orderMEs(MEs0)
  write.table(MEs, paste0(path, '/', 'MEs - ', name, '.txt'), sep='\t', row.names=T, col.names=T, quote=F)
  moduleTraitCor = cor(MEs, as.matrix(datTraits), use="p")
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
  pdf(paste0(path, '/', 'Module-Trait-Correlation - ', name, '.pdf'), wi=8, he=11)
  textMatrix=paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep="")
  dim(textMatrix) = dim(moduleTraitCor)
  par(mar=c(6,8.5,3,3))
  labeledHeatmap(Matrix=moduleTraitCor,
                 xLabels=names(datTraits),
                 yLabels=names(MEs),
                 ySymbols=names(MEs),
                 colorLabels=FALSE,
                 colors=blueWhiteRed(50),
                 textMatrix=textMatrix,
                 setStdMargins=FALSE,
                 cex.text=0.5,
                 zlim=c(-1,1),
                 main=paste("Module-trait relationships -", name))
  dev.off()
  modNames = substring(names(MEs), 3)
  geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
  names(geneModuleMembership) = paste("MM", modNames, sep="")
  names(MMPvalue) = paste("p.MM", modNames, sep="")
  geneTraitSignificance = as.data.frame(sapply(datTraits, FUN = function(x) cor(datExpr, x, use="p")))
  rownames(geneTraitSignificance) <- colnames(datExpr)
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
  names(geneTraitSignificance) = paste("GS.", names(datTraits), sep="")
  names(GSPvalue) = paste("p.GS.", names(datTraits), sep="")

  # Use biomaRt to annotate Ensembl IDs with HUGO gene symbol and entrez ID.
  tryCatch({
    source('/research/labs/moleneurosci/bug/projects/single-cell-analysis/util/load_biomart.r')
    ensembl <- useEnsembl(biomart, dataset=dataset, host = host)
    bm <- getBM(attributes=attributes, mart=ensembl)
    probes = colnames(datExpr)
    probes2bm = match(probes, bm[[attribute_given]])
    sum(is.na(probes2bm))
    geneInfo0 = data.frame(Ensemble_ID = probes,
                           geneSymbol = bm$external_gene_name[probes2bm],
                           Chr = bm$chromosome_name[probes2bm],
                           Start= bm$start_position[probes2bm],
                           Stop= bm$end_position[probes2bm],
                           LocusLinkID = bm$entrezgene[probes2bm],
                           moduleColor = moduleColors,
                           geneTraitSignificance,
                           GSPvalue)
# Order modules by their significance for the first trait of datTraits
    modOrder = order(-abs(cor(MEs, datTraits[,1], use='p')))
# Add module membership information in the chosen order
    for (mod in 1:ncol(geneModuleMembership))
    {
        oldNames = names(geneInfo0)
      geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                                                      MMPvalue[, modOrder[mod]]);
        names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                                                    paste("p.MM.", modNames[modOrder[mod]], sep=""))
    }
# Order the genes in the geneInfo variable first by module color, then by #geneTraitSignificance
    geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0[[paste0('GS.', names(datTraits)[1])]]));
    geneInfo = geneInfo0[geneOrder, ]
    write.csv(geneInfo, file = paste0(path, '/', "geneInfo - ", name, ".csv"), row.names=F, quote=F)
  }, error = function(cond) {
    obtain_module_membership(datExpr, datTraits, MEs, moduleColors, path, name)
  })

  tryCatch({
    # Annotate modules using anRichment package
    LocusLinkID = bm$entrezgene[probes2bm]

    ## Run enrichment analysis
    library("anRichment")
    #unloadNamespace('dplyr') # unload dplyr
    if (is.null(GOcollection)) {
      GOcollection = buildGOcollection(organism = organism)
    }
    GOenrichment=enrichmentAnalysis(classLabels=moduleColors,
                                    identifiers=LocusLinkID,
                                    refCollection=GOcollection,
                                    useBackground="given",
                                    threshold=1e-4,
                                    thresholdType=c("Bonferroni", "FDR", "nominal"),
                                    nBestDataSets=10,
                                    getOverlapEntrez=FALSE,
                                    getOverlapSymbols=TRUE,
                                    maxReportedOverlapGenes=100,
                                    ignoreLabels="grey",
                                    getFDR=TRUE,
                                    getBonferroniCorrection=TRUE)
    collectGarbage()
    write.csv(GOenrichment$enrichmentTable, file=paste0(path, '/', "AnRichment-enrichmentTable - ", name, ".csv"), row.names=FALSE, quote=F)
  }, error = function(cond) {
    print("enrichmentAnalysis errored. Recommend using clusterProfiler")  
  })
}
