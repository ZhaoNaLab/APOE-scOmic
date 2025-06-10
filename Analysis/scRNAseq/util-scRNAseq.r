library(tidyverse)
library(magrittr)
library(ggplot2)
library(patchwork)
library(Seurat)

celltypes = c('MG', 'AS', 'EC', 'CP', 'Lym', 'Neu', 'OL', 'PC', 'SC', 'NP')
genotype_colors <- c(
  E2 = RColorBrewer::brewer.pal(n=8,name='Set1')[[2]],
  E3 = RColorBrewer::brewer.pal(n=8,name='Set1')[[3]],
  E4 = RColorBrewer::brewer.pal(n=8,name='Set1')[[1]]
)
age_colors <- c(
  `3M` = RColorBrewer::brewer.pal(n=8,name='Set3')[[1]],
  `24M` = RColorBrewer::brewer.pal(n=8,name='Set3')[[2]]
)
sex_colors <- c(
  `M` = RColorBrewer::brewer.pal(n=8,name='Set2')[[1]],
  `F` = RColorBrewer::brewer.pal(n=8,name='Set2')[[2]]
)
flowcell_colors <- c(
  `1` = RColorBrewer::brewer.pal(n=8, name='Dark2')[[1]],
  `2` = RColorBrewer::brewer.pal(n=8, name='Dark2')[[2]]
)


load_sample_info <- function(si = readxl::read_xlsx(path='scRNAseq/scRNAseq-metadata.xlsx', sheet='INPUT-FOR-METADATA')) {
  si <- si[1:24,]
  colnames(si)[1] <- 'orig.ident'
  si$Collection.Date <- str_split_fixed(si$`Sample Collection Date and Time...4`, ' ', 2)[,1]
  si <- si %>%
    dplyr::select(orig.ident, Collection.Date, Gender, APOE, age.in.months, `Flow Cell`, Lane) %>%
    dplyr::rename(sex = Gender, flow.cell = `Flow Cell`) %>%
    dplyr::mutate(flow.cell = as.factor(flow.cell),
                  Lane = as.factor(Lane),
                  APOE = as.factor(APOE))
  return(si)
}


# Rename or relevel factors
sample_info_factors <- function(x) {
  info_helper <- function(.x) { 
    .x %>% 
      mutate(
        orig.ident = factor(orig.ident, paste0('BU-APOE-', 1:24)),
        APOE = factor(paste0('E', APOE), c('E3', 'E2', 'E4')),
        age.in.months = factor(paste0(age.in.months, 'M'), levels=c('3M', '24M')),
        sex = factor(sex, c('M', 'F')),
        flow.cell = factor(flow.cell),
        APOE.age = factor(paste0(APOE, '_', age.in.months), levels=c('E2_3M', 'E2_24M', 'E3_3M', 'E3_24M', 'E4_3M', 'E4_24M')),
        APOE.age.sex = factor(paste0(APOE, '_', age.in.months, '_', sex), levels=c('E2_3M_M', 'E2_3M_F', 'E2_24M_M', 'E2_24M_F',
                                                                                   'E3_3M_M', 'E3_3M_F', 'E3_24M_M', 'E3_24M_F',
                                                                                   'E4_3M_M', 'E4_3M_F', 'E4_24M_M', 'E4_24M_F'
                                                                                   )),
        SampleID = factor(paste0(stringr::str_extract(orig.ident, '[0-9]+'), '.', APOE.age.sex))
      )
  }
  
  sample_info <- load_sample_info() %>% info_helper
  x <- info_helper(x)
  x %>% mutate(SampleID = factor(SampleID, levels = sample_info$SampleID))
}


load_roundtwo_integrated_celltype <- function(celltype=NULL) {
  library(batchtools)
  reg <- loadRegistry("ExperimentRegistries/reg-integration")

  # Helper function to get the job.id
  jobid_celltype <- function(celltype) {
    findExperiments(prob.pars = cell_type == celltype, algo.pars = length(vars.to.regress) == 3 & all(c('percent.mt', 'percent.rp', 'sex') %in% vars.to.regress) & batch == 'flow.cell', reg=reg_roundtwo)
  }

  # Helper function for re-doing dimensionality reduction
  dimreduc_celltype <- function(x, celltype) {
    id <- jobid_celltype(celltype)
    resolution <- unwrap(getJobPars(id))[,resolution][[1]]
    x %>%
      RunPCA() %>%
      FindNeighbors(dims=1:15) %>%
      FindClusters(resolution = resolution) %>%
      RunUMAP(reduction='pca', dims=1:15)
  }

  # Helper function for removing unwanted subclusters of celltype result
  subset_celltype <- function(x, celltype) {
    if (celltype == 'AS') {
      x <- subset(x, seurat_clusters != 7)
    } else if (celltype == 'CP') {
      x <- subset(x, seurat_clusters != 7)
    } else if (celltype == 'NP') {
      x <- subset(x, seurat_clusters == 10)
    } else if (celltype == 'MG') {
      x <- subset(x, seurat_clusters != 10)
    }
    x$seurat_clusters <- droplevels(x$seurat_clusters)
    return(x)
  }

  # If celltype is NULL, load the entire list
  if (is.null(celltype)) {
    ids <- findExperiments(algo.pars = length(vars.to.regress) == 3 & all(c('percent.mt', 'percent.rp', 'sex') %in% vars.to.regress) & batch == 'flow.cell')
    x <- reduceResultsList(ids)
    names(x) <- unlist(unwrap(getJobPars(ids))[,cell_type])
    for (celltype in names(x)) {
      x[[celltype]] <- dimreduc_celltype(x[[celltype]], celltype)
    }
    x[['NP']] <- subset_celltype(x[['MG']], 'NP')
    for (celltype in c('AS', 'CP', 'MG')) {
      x[[celltype]] = subset_celltype(x[[celltype]], celltype)
    }
  } else {
    if (celltype == 'NP') {
      id <- jobid_celltype('MG')
    } else {
      id <- jobid_celltype(celltype)
    }
    x <- loadResult(id)
    x <- dimreduc_celltype(x, celltype)
    x <- subset_celltype(x, celltype)
  }
  return(x)
}


load_roundtwo_integrated <- function() {
  library(batchtools)
  reg = loadRegistry('ExperimentRegistries/reg-integration')
  job.id = findExperiments(prob.name = 'reintegrate.celltypes')
  loadResult(job.id)
}


load_markergenes <- function(
  celltype = NULL,
  latent_vars = c('flow.cell', 'sex', 'cngeneson')
  ) {
  library(batchtools)
  reg_markers = loadRegistry('ExperimentRegistries/findmarkers', writeable=TRUE)
  
  if (is.null(celltype)) {
    job_id = findExperiments(findDone(), prob.name='apoetr', algo.pars = all(latent.vars == latent_vars), reg=reg_markers) 
  } else {
    job_id = findExperiments(findDone(), prob.name='apoetr.subcluster', prob.pars = Cell_type == celltype, algo.pars = all(latent.vars == latent_vars), reg=reg_markers) 
  }
  
  loadResult(job_id)
}


create_group_names <- function(meta.data) {
  meta.data %>%
    mutate(
      orig.ident = factor(orig.ident, paste0('BU-APOE-', 1:24)),
      APOE = paste0('E', APOE),
      age.in.months = paste0(age.in.months, 'M'),
      flow.cell = paste0('FC', flow.cell),
    )
}

#'
#' x : SeuratObject, is the given celltype
#' celltype : string, is the celltype of given SeuratObject x.
average_expression <- function(x, celltype, assays='integrated', slot='scale.data') {
  dat = AverageExpression(x, 
                          assays=assays, 
                          slot=slot,
                          group.by=c('APOE', 'age.in.months','sex', 'orig.ident')
  )$integrated
  dat = rbind(
    celltype   = celltype, 
    APOE       = str_split_fixed(colnames(dat), '_', 4)[,1], 
    age        = str_split_fixed(colnames(dat), '_', 4)[,2],
    sex        = str_split_fixed(colnames(dat), '_', 4)[,3],
    orig.ident = str_split_fixed(colnames(dat), '_', 4)[,4],
    dat)
  dat = as.data.frame(dat) %>%
    tibble::rownames_to_column('gene')
  return(dat)
}
