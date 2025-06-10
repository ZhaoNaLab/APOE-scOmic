WRAPPERS='scRNAseq/util'
source(file.path(WRAPPERS, 'transfer_wrappers.r'))
source(file.path(WRAPPERS, 'multiplet_wrappers.r'))
source("scRNAseq/util/future.r") # setup_future

# Generalized adjustment and batch correction.
adjcor_seu <- function(data, job, instance, assay=instance@active.assay, vars.to.regress=NULL, adjust.method=NULL, batch=NULL, batch.method=NULL, reduction.method='umap', cluster.method=NULL, nfeatures=3000, workers=NULL, maxSize=NULL, dims=1:10, reference=NULL, resolution=0.8, adjust.in.parallel=TRUE, plan='multisession', ...) {
  source('scRNAseq/seu.r')
  DefaultAssay(object=instance) <- assay
  print(instance)
  library(future.apply)

  if (is.null(batch)) { # | is.null(batch.method)) {
    x <- adjust_seu(x=instance, assay=assay, vars.to.regress=vars.to.regress, adjust.method=adjust.method)
  } else {
    rlog::log_info("seu.r adjcor_seu: Applying normalization / scaling (adjustment) by batch...")
    instance@meta.data[[batch]] <- as.factor(x=instance@meta.data[[batch]])
    x <- SplitObject(object=instance, split.by=batch)
    workers = ifelse(is.null(workers), length(x), workers)
    
    # 2024-02-23: "linear" method (ie ScaleData) doesn't need adjustment berfore integration. That step comes after. However SCTransform scales data before integration and not after.
    if (!is.null(adjust.method)) {
      if (adjust.method == 'SCT') {
        if (adjust.in.parallel) {
          setup_future(x, plan=plan, workers=workers)
          rlog::log_info(paste('seu.r adjcor_seu: Using multisession future with', workers, 'workers.'))
          x <- future_lapply(x, function(y) sct_seu(x=y, assay=assay, vars.to.regress=vars.to.regress, workers=workers), future.seed=TRUE)
          # now that x is bigger, we need to update the size of the global
          options(future.globals.maxSize=object.size(x)*2)
        } else {
          x <- lapply(x, function(y) sct_seu(x=y, assay=assay, vars.to.regress=vars.to.regress, workers=workers))
        }
        for (i in 1:length(x)) {
          levels(x[[i]][['SCT']]) <- names(x)[[i]]
        }
      } else if (adjust.method == 'linear') {
        x <- lapply(x, NormalizeData)
        x <- lapply(x, FindVariableFeatures)
        x <- future_lapply(x, linear_seu, assay=assay, vars.to.regress=vars.to.regress)
      }
    }

    # input x is a list split by batch
    # output x is re-merged SeuratObject
    x <- batch_seu(x, assay=assay, adjust.method=adjust.method, batch=batch, batch.method=batch.method, vars.to.regress=vars.to.regress, reference=reference, workers=workers, dims=dims, ...)

    if (!is.null(batch.method)) {
      if(batch.method == 'seurat' & adjust.method == 'linear') {
        x <- linear_seu(x, assay=assay, vars.to.regress=vars.to.regress, split.by=batch)
      }
    }
  }
  x <- reduction_seu(x, batch.method=batch.method, reduction.method='pca', dims=dims)
  x <- clustering_seu(x, adjust.method=adjust.method, batch.method=batch.method, cluster.method=cluster.method, dims=dims, resolution=resolution)
  x <- reduction_seu(x, batch.method=batch.method, reduction.method=reduction.method, dims=dims)
  return(x)
}

adjust_seu <- function(x, assay=x@active.assay, vars.to.regress=NULL, adjust.method=NULL, workers=NULL) {
  print("Now running adjust_seu ...") 
  if (is.null(adjust.method)) {
    return(x)
  } 
  if (adjust.method == 'linear') {
    return(linear_seu(x=x, assay=assay, vars.to.regress=vars.to.regress))
  } else if (adjust.method == 'combat') {
    return(x) 
  } else if (adjust.method == 'SCT') {
    return(sct_seu(x=x, assay=assay, vars.to.regress=vars.to.regress, workers=workers))
  }
}

# Linearly adjust covariates.
linear_seu <- function(x, vars.to.regress, assay=x@active.assay, split.by=NULL) {
  #DefaultAssay(object=x) <- assay
  #x <- NormalizeData(object=x, assay=x@active.assay)
  #x <- FindVariableFeatures(x, assay=assay)
  x <- ScaleData(object=x, vars.to.regress=vars.to.regress, assay=x@active.assay, split.by=split.by)
  #x <- FindVariableFeatures(object=x) # 2024-02-22 uhhhhhh yeah why are we doing this AFTER scaledata?
  return(x)
}

# Adjust with SCTransform
sct_seu <- function(x, assay=x@active.assay, vars.to.regress, workers) {
  DefaultAssay(object=x) <- assay
  library(glmGamPoi)
  library(future)
  #if (nbrOfFreeWorkers() < 1) {
  #  plan(sequential)
  #  rlog::log_info('seu.r sct_seu: Using sequential plan ...')
  #} else {
  #  plan(multicore, workers=workers) # 04/14/2023 -- removing this because it slows things down.
  #  options(future.globals.maxSize=object.size(x)*2)
  #  rlog::log_info(paste('seu.r sct_seu: Using multicore future with', workers, 'workers.'))
  #}
  x <- SCTransform(object=x, assay=x@active.assay, vars.to.regress=vars.to.regress, method='glmGamPoi', vst.flavor='v2')
  return(x)
}

batch_seu <- function(x, assay, adjust.method=NULL, batch=NULL, batch.method=NULL, vars.to.regress=NULL, nfeatures=3000, reference=NULL, workers=4, dims=1:30, ...) {
  rlog::log_info('seu.r batch_seu: began batch_seu.')
  if (is.null(batch)) {
  #if (is.null(batch.method) | is.null(batch)) {
    return(x)
  }
  if (length(dims) == 1) { dims = 1:dims }
  if (is.null(batch.method)) {
  warning("batch.method set to NULL. Running FindVariableFeatures on each batch, SelectIntegrationFeatures, then merging.")
    x <- lapply(x, FindVariableFeatures)
    var.features <- SelectIntegrationFeatures(x)
    x <- Reduce(merge, x)
    VariableFeatures(x) <- var.features
  } else if (batch.method == 'harmony') {
    if (adjust.method == 'linear') {
      var.features <- SelectIntegrationFeatures(object.list=x, nfeatures=nfeatures)
      scale.dat <- Reduce(cbind, lapply(x, function(y) GetAssayData(y, slot='scale.data')))
      x <- Reduce(merge,x)
      scale.dat <- scale.dat[rownames(x), colnames(x)]
      x <- SetAssayData(object=x, slot='scale.data', new.data=scale.dat)
      VariableFeatures(object=x, assay=assay) <- var.features
    } else if (adjust.method == 'SCT') {
      var.features <- SelectIntegrationFeatures(object.list=x, nfeatures=nfeatures)
      x <- Reduce(f=merge, x=x)
      VariableFeatures(object=x, assay='SCT') <- var.features
    }
    x <- harmony_seu(x=x, batch=batch)
  } else if (batch.method == 'seurat') {
    if (is.null(adjust.method) | adjust.method == 'linear') {
      features <- SelectIntegrationFeatures(object.list=x, nfeatures=nfeatures)
      normalization.method = 'LogNormalize'
    } else if (adjust.method == 'SCT') {
      features <- SelectIntegrationFeatures(object.list=x, nfeatures=nfeatures)
      x <- PrepSCTIntegration(object.list=x, anchor.features=features)
      normalization.method = 'SCT'
    }
    x <- integrate_seu(x=x, batch=batch, features=features, normalization.method=normalization.method, reference=reference, dims=dims, workers=workers, ...)
  } else if (batch.method == 'combat') {
    x <- combat_seu(x=x, batch=batch, vars.to.regress=vars.to.regress)
    x <- FindVariableFeatures(object=x)
  }
  return(x)
}
# Correct batch effect with ComBat
combat_seu <- function(x, batch, vars.to.regress) {
  dat <- GetAssayData(object=x, slot='counts')
  f <- as.formula(paste0('~1+', paste0(vars.to.regress, collapse='+')))
  mod <- model.matrix(formula=f, data=x@meta.data)
  cdat <- sva::ComBat(dat=dat, batch=batch, mod=mod)
  x <- SetAssayData(object=x, slot='scale.data', new.data = cdat)
  print('Returned ComBat corrected data to scale.data slot of RNA assay.')
  return(x)
}

# Batch correction with harmony.
harmony_seu <- function(x, batch) {
  x <- RunPCA(object=x)
  x <- harmony::RunHarmony(object=x,  group.by.vars=batch, reduction='pca', assay.use=x@active.assay)
  return(x)
}

# Batch correction with seurat.
integrate_seu <- function(x, batch, features, normalization.method, reference=NULL, workers=4, dims=1:30, ...)  {
  rlog::log_info("integrate_seu")
  setup_future(x, plan='multicore', workers=workers)
  if (length(dims) == 1) { dims = 1:dims }
  rlog::log_info(paste('Using dims:', paste0(dims,collapse=',')))
  # TODO: 2024-02-22 FindIntegrationAnchors, scale=FALSE: because we almost certainly scaled the data prior
  anchors <- FindIntegrationAnchors(object.list=x, anchor.features=features, normalization.method=normalization.method, reference=reference, dims=dims, scale=FALSE, ...)
  setup_future(c(anchors, x), plan='multicore', workers=workers)
  x <- IntegrateData(anchorset=anchors, normalization.method=normalization.method, dims=dims, ...)
  return(x)
}

# perform dimensionality reduction
reduction_seu <- function(x, batch.method=NULL, reduction.method='umap', dims=1:10) {
  rlog::log_info("reduction_seu")
  if (length(dims) == 1) { dims = 1:dims }
  rlog::log_info(paste('Using dims:', dims))
  if (reduction.method == 'pca') { # previously i checked to see if pca present. i'm okay with recomputing by default.
    if (length(VariableFeatures(x)) == 0) {
      x <- FindVariableFeatures(x) 
    }
    if (ncol(GetAssayData(x, slot='scale.data')) == 0) {
      x <- ScaleData(x)
    }
    x <- RunPCA(object=x)
    return(x)
  }
  # Compute a UMAP on desired dim reduc
  if (is.null(batch.method)) {
    reduction='pca'
  } else if (batch.method == 'harmony') {
    reduction='harmony'
  } else {
    reduction='pca'
  }
  if (reduction.method == 'umap') {
    if (length(dims) == 1) {  dims = 1:dims }
    x <- RunUMAP(object=x, reduction=reduction, dims=dims)
    x$UMAP_1 <- Embeddings(object=x, reduction=reduction.method)[,1]
    x$UMAP_2 <- Embeddings(object=x, reduction=reduction.method)[,2]
  }
  return(x)
}

# perform clustering
clustering_seu <- function(x, adjust.method=NULL, batch.method=NULL, cluster.method=NULL, dims=1:10, resolution=0.8) {
  rlog::log_info("clustering_seu")
  if (length(dims) == 1) { dims = 1:dims }
  rlog::log_info(paste('Using dims:', dims))
  if (is.null(adjust.method)) {
    adjust.method='RNA'
  }
  if (is.null(batch.method)) {
    batch.method = ''
  }
  x <- FindNeighbors(object=x, graph.name=paste0(paste0(adjust.method, '.', batch.method), c('_nn', '_snn')), dims=dims)
  for (res in resolution) {
    x <- FindClusters(object=x, graph.name=paste0(adjust.method, '.', batch.method, '_snn'), resolution=res)
  }
  return(x)
}


