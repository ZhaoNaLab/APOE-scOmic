# This script is for replicating analysis of the scRNAseq mouse brain APOE dataset.
# We begin with a merged Seurat object created from the filtered 10X secondary outputs, saved in 'Datasets/apoetr.rds'
# The below code utilizes the batchtools R package to help parallelize the analysis.
# See https://batchtools.mlr-org.com/ for details on how to use batchtools.
# A separate script will be made available to parallelize without the usage of batchtools.
# These scripts were run using R version 4.0.3, Seurat v4.

# At a high level, analysis was done with two rounds of integration, which are defined by separate `addProblem` statements below. 3a however is an additional step that is defined by a separate function in `util_for_manuscript.r`
# 1. integrate.filtered: Initial round of integration
# 2. apoetr.subcluster:  Integration within each celltype found from clusters of integrate.filtered result.
# 3. remove_and_reintegrate: Re-do integration within each celltype, after removing doublet or low-quality clusters.
# 3a. load_roundtwo_integrated_celltype: After visualizing the results of `remove_and_reintegrate`, we removed further clusters and repeated dimensionality reduction and clustering.
# 4. reintegrate.celltypes:  Merge results of `load_roundtwo_integrated_celltype`, re-perform integration to compute the final, merged dataset.

library(batchtools)
source('scRNAseq/util/batchtools-wrappers.r') # for loading the adjcor_seu algorithm definition
makeExperimentRegistry('ExperimentRegistries/reg-integration', packages = c("Seurat", "batchtools", "magrittr", "dplyr", "ggplot2")) # Run only once (then comment out)
reg = loadRegistry('ExperimentRegistries/reg-integration', writeable=TRUE)
addAlgorithm('seurat', adjcor_seu, reg = reg)

# 1. Initial round of integration
addProblem('integrate.filtered', data = NA, fun = function(data, job, percent.mt.t, threshold) {
  if (threshold > 0.5 & threshold > 0 & threshold < 1) {
    threshold = 1-threshold
  } 
  apoetr <- readRDS('Datasets/apoetr.rds')
  apoetr.split <- SplitObject(apoetr, split.by='flow.cell')
  apoetr.split <- lapply(apoetr.split, function(x) {
    subset(x, subset = percent.mt < percent.mt.t &
           nFeature_RNA < quantile(apoetr$nFeature_RNA, 1-threshold) & 
           nFeature_RNA > quantile(apoetr$nFeature_RNA, threshold) &
           nCount_RNA < quantile(apoetr$nCount_RNA, 1-threshold) &
           nCount_RNA > quantile(apoetr$nCount_RNA, threshold)
    )}
  )
  apoetr <- merge(apoetr.split[[1]], apoetr.split[2:length(apoetr.split)])
},
reg = reg)

addExperiments(
  prob.designs = list(integrate.filtered = CJ(percent.mt.t = 20, threshold = 0.05)),
  algo.designs = list(seurat = CJ(
    adjust.method = 'SCT',
    batch.method = 'seurat',
    dims = 30,
    nfeatures = 3000,
    workers = 2,
    maxSize = 50 * 1000^3,
    reference = 1,
    vars.to.regress = 'percent.mt',
    batch = 'orig.ident'
    #sorted=FALSE)
    )
  ))

# Uncomment to run the job
#submitJobs(reg = reg) # adjust resources as needed. Job can last up to 14 hours.


# 2. Integration within each celltype, to identify spurious clusters of cells to remove manually.
addProblem('apoetr.subcluster', data=NA, fun=function(data, job, celltype){
  reg = batchtools::loadRegistry('ExperimentRegistries/reg-integration')
  job.id = batchtools::findExperiments(prob.name = 'integrate.filtered', 
                                       prob.pars = threshold == 0.05 & percent.mt.t == 20,
                                       reg = reg)
  re = batchtools::loadResult(job.id) # load integration result after mitochondria % filter
  DefaultAssay(re) <- 'RNA'
  re = DietSeurat(re, assays='RNA')
  re = NormalizeData(re) 
  re <- subset(re, subset = seurat_clusters %in% c(0:20, 22:24)) # Get rid of doublet clusters
  re@meta.data %<>% mutate(
    cell_type = case_when(
      seurat_clusters %in% c(0, 5, 7, 8, 4, 9, 19, 2, 23, 11, 22) ~ 'MG',
      seurat_clusters %in% c(1, 6, 24) ~ 'AS',
      seurat_clusters %in% c(10) ~ 'OL',
      seurat_clusters %in% c(3, 12, 20) ~ 'EC',
      seurat_clusters %in% c(17) ~ 'PC',
      seurat_clusters %in% c(18) ~ 'SC',
      seurat_clusters %in% c(13, 15) ~ "Neu",
      seurat_clusters %in% c(14) ~ 'CP',
      seurat_clusters %in% c(16) ~ 'Lym'
    )
  ) %>% as.data.frame(row.names=colnames(re))
  re = subset(re, subset=cell_type == celltype)
  re$percent.rp = PercentageFeatureSet(re, pattern='^Rp(s|l)')
  return(re)
}, reg = reg)

addExperiments(prob.designs = list(apoetr.subcluster = CJ(celltype=c('MG', 'AS', 'Lym', 'Neu', 'OL', 'EC', 'PC', 'SC', 'CP'))),
               algo.designs=list(seurat = CJ(
                 adjust.method=list('SCT'),
                 batch.method='seurat',
                 dims=30,
                 nfeatures=3000,
                 workers=2,
                 maxSize=50 * 1000^3,
                 reference=1,
                 vars.to.regress=list(c('percent.mt', 'sex', 'percent.rp')),
                 batch='flow.cell',
                 sorted=FALSE
               )),
               reg = reg)

# Uncomment to run the job
#submitJobs(reg = reg) # adjust resources as needed. MG job can last up to 6 hours.


# 3. Re-integrate the cells after removing spurious clusters.
addProblem('remove_and_reintegrate', data=NULL, fun=function(data, job, cell_type, resolution, to_remove, ...){
  reg = batchtools::loadRegistry('ExperimentRegistries/reg-integration/')
  job_id = findExperiments(prob.name='apoetr.subcluster', 
                           prob.pars = celltype == cell_type, 
                           reg=reg) %>% unlist
  re = batchtools::loadResult(id=job_id, reg=reg)
  re = FindClusters(re, resolution=resolution, graph.name='SCT.seurat_snn')
  to_keep = unique(re$seurat_clusters)
  to_keep = to_keep[!to_keep %in% to_remove]
  re = subset(re, seurat_clusters %in% to_keep)
  print(re)
  print(paste0('Clusters to keep: ', unique(re$seurat_clusters)))
  DefaultAssay(re) <- 'RNA'
  re = DietSeurat(re, assays='RNA')
  re$APOE.cell = paste0('E', re$APOE, '_', 'F', re$flow.cell)
  return(re)
}, reg = reg)

addExperiments(prob.designs=list(remove_and_reintegrate=data.table(
  cell_type = c('MG', 'AS', 'EC', 
                'CP', 'Lym', 
                'Neu', 'OL', 
                'PC', 'SC'),
  resolution = c(0.3, 0.3, 0.2, 
                 0.3, 0.2, 
                 0.3, 0.3, 
                 0.3, 0.3),
  to_remove = list(c(10, 11, 15), c(4, 8), c(6,9,3,10), 
                   c(4,7,6), c(1,9), 
                   c(4,5,6,7,8,10), c(4,14,15,11), 
                   c(1,3), c(3,4))
  )),
  algo.designs = list(seurat = CJ(
    adjust.method = 'SCT',
    batch.method = 'seurat',
    dims = 15,
    nfeatures = 3000,
    workers = 6,
    maxSize = 50 * 1000^3,
    reference = 1,
    vars.to.regress = list(c('percent.mt', 'percent.rp', 'sex')),
    batch = 'flow.cell',
    sorted = FALSE)
  ),
  reg = reg
)

# Uncomment to run the job
#submitJobs(reg = reg) # can take up to 6h for MG

# 3a. Define `load_roundtwo_integrated_celltype` in 'scRNAseq/util/util_scRNAseq.r', a function that takes results from above, further identifies clusters that appear spurious, and then re-performs simple dimensionality reduction and clustering one last time.


# 4. Final merging and integration of the cells following above quality control steps.
addProblem('reintegrate.celltypes', data=NA, fun=function(data, job, ...) {
  source("scRNAseq/util_scRNAseq.r") # load_roundtwo_integrated_celltype()
  x = load_roundtwo_integrated_celltype()
  x = lapply(names(x), function(n) {
    DefaultAssay(x[[n]]) <- 'RNA'
    DietSeurat(x[[n]],counts=TRUE,data=TRUE,assays=c('RNA'))
    x[[n]]$cell_type = n
    return(x[[n]])
    })
  x = Reduce(merge, x)
  return(x)
}, reg=reg)


addExperiments(prob.designs = list(reintegrate.celltypes = CJ()),
               algo.designs = list(seurat = CJ(
                 adjust.method = 'SCT',
                 batch.method = 'seurat',
                 dims = 30,
                 nfeatures = 3000,
                 workers = 2,
                 maxSize = 50 * 1000^3,
                 reference = 1,
                 vars.to.regress = list(c('percent.mt', 'percent.rp','sex')), 
                 batch = 'flow.cell',
                 sorted = FALSE)),
               reg=reg)

# Uncomment to run the job
#submitJobs(reg = reg) # ~14 hours
