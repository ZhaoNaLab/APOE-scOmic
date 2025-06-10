load_degs <- function(
  celltype,
  model = c('SEX+APOExAGE', 'AGE+APOExSEX', 'APOE+AGExSEX'),
  mode = c('RNA', 'protein'),
  method = c('edgeR', 'limma', 'limma-trend log', 'limma-trend median', 'limma-trend quantile'),
  interaction = TRUE
) {
  model = model[1]
  mode = mode[1]
  method = method[1]
  if (mode == 'RNA') {
    file.dir = 'Results/scRNAseq/'
  } else {
    file.dir = paste0('Results/proteomics/', method) #, '-DEG-analysis')
  }
  comparisons = list.files(file.path(file.dir, model, celltype), pattern = '[A-Z]+', full.names = TRUE)
  if (!interaction) {
    comparisons = comparisons[!stringr::str_detect(stringr::str_remove(comparisons, paste0(stringr::str_replace(model, '\\+', '\\\\+'), '/')), 'x')]
  }
  file.list = lapply(comparisons, function(f) list.files(path=f,  pattern='.txt',  recursive = T, full.names=T)) %>% unlist
  deg.list = lapply(file.list, 
                    function(f) { 
                      tmp = data.table::fread(f) %>% dplyr::rename(gene=V1)
                      colnames(tmp) = paste(mode, colnames(tmp), sep='_')
                      tmp$mode = mode
                      tmp$method = method
                      return(tmp)
                    })
  names(deg.list) = stringr::str_remove( stringr::str_remove( file.list,  '.txt'),  '^.*/') 
  return(deg.list)
}

inner_degs <- function(
  rna,
  pro,
  comparison = intersect(names(rna), names(pro))[1]
) {
  if(length(comparison) == 2) {
    rna_comparison = comparison[1]
    pro_comparison = comparison[2]
  } else if (length(comparison) == 1 | length(comparison) > 2) {
    rna_comparison = pro_comparison = comparison
  } else {
    stop('No comparisons.')
  }
  dplyr::bind_cols(
    rna[[rna_comparison]] %>% dplyr::select(!c(mode, method)) %>% dplyr::filter(RNA_gene %in% pro[[pro_comparison]]$protein_gene) %>% arrange(RNA_gene),
    pro[[pro_comparison]] %>% dplyr::select(!c(mode, method)) %>% dplyr::filter(protein_gene %in% rna[[rna_comparison]]$RNA_gene) %>% arrange(protein_gene)
  )
}
