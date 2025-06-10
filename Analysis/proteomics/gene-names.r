# This script helps get gene names for peptides in the proteomics data.
# It utilizes the UniProt mouse fasta file which was downloaded from https://www.uniprot.org/proteomes/UP000000589 on March 15, 2023

library(org.Mm.eg.db)
library(AnnotationDbi)
library(biomaRt)
source("Analysis/proteomics/util.r")

pr = load_excel()
accessions = union(as.data.frame(pr$MG) %>% get_accessions %>% dplyr::select(Accession) %>% unlist %>% unname,
                  as.data.frame(pr$AS) %>% get_accessions %>% dplyr::select(Accession) %>% unlist %>% unname)
#ensembl = select(org.Mm.eg.db, keys = keys(org.Mm.eg.db), columns = c('ENSEMBL', 'SYMBOL', 'UNIPROT'))
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = 'https://useast.ensembl.org')

mouse_attributes = getBM(
  attributes = c(
    'external_gene_name',
    'uniprot_isoform',
    'uniprotsptrembl',
    'uniprotswissprot'
  ), 
  mart = mouse)
more_attributes = getBM(
  attributes = c(
    'external_gene_name',
    'uniprot_gn_symbol',
    'uniprot_gn_id'
  ), 
  mart = mouse)

# ---- Uniprot fasta ----
library(seqinr)
fasta_to_table <- function(fasta) {
  library(stringr)
  uniprot_table = data.frame(
    name = lapply(fasta, attr, which='name') %>% unlist,
    annot = lapply(fasta, attr, which='Annot') %>% unlist
  )
  uniprot_table <- uniprot_table %>%
    mutate(
      GN = str_extract(annot, 'GN=.* ') %>%
             str_remove('^GN=') %>% 
             str_remove(' .*') %>%
             str_trim()
    ) %>%
    tidyr::separate(
      col=name,
      sep='\\|',
      into = c(NA, 'Accession1', 'Accession_Mouse')
    ) %>%
    mutate(
      Accession2 = str_remove(Accession_Mouse, '_[A-Z].*$')
    )
}

uniprot = read.fasta(file = 'Datasets/UniProt/mouse-2023.03.15-19.15.47.29.fasta')
uniprot_table <- fasta_to_table(uniprot)

missing = !(
  accessions %in% uniprot_table$Accession1 |
  accessions %in% uniprot_table$Accession2 |
  accessions %in% mouse_attributes$uniprot_isoform | 
  accessions %in% mouse_attributes$uniprotsptrembl |
  accessions %in% mouse_attributes$uniprotswissprot |
  accessions %in% more_attributes$uniprot_gn_id
)
summary(missing)

accessions[missing]

missing.fasta = read.fasta("Datasets/Proteomics/missing-accessions.fasta")
missing_table = fasta_to_table(missing.fasta)
missing_table <- missing_table %>% filter(!str_detect(Accession1, '-[0-9]$'))
found = c(F,F,T,T,T,T,F,F,F,F,F,T,T,F,T,F,T,T,T,T,T,T,T,F,T,T,T) # double-check!

missing_accessions = data.frame(
  missing = accessions[missing][found],
  found = missing_table$Accession1,
  GN = missing_table$GN
)


get_gene_name <- function(accessions, uniprot_table, mouse_attributes, missing_accessions) {
  genes = uniprot_table[match(accessions, uniprot_table$Accession1), 'GN']
  genes[is.na(genes)] = uniprot_table[match(accessions[is.na(genes)], uniprot_table$Accession2), 'GN']
  genes[is.na(genes)] = mouse_attributes[match(accessions[is.na(genes)], mouse_attributes$uniprot_isoform), 'external_gene_name']
  genes[is.na(genes)] = mouse_attributes[match(accessions[is.na(genes)], mouse_attributes$uniprotsptrembl), 'external_gene_name']
  genes[is.na(genes)] = mouse_attributes[match(accessions[is.na(genes)], mouse_attributes$uniprotswissprot), 'external_gene_name']
  genes[is.na(genes)] = mouse_attributes[match(accessions[is.na(genes)], mouse_attributes$uniprot_gn_id), 'external_gene_name']
  genes[is.na(genes)] = missing_accessions[match(accessions[is.na(genes)], missing_accessions$missing), 'GN']
  names(genes) = accessions
  return(genes)
}
genes = get_gene_name(accessions, uniprot_table, mouse_attributes, missing_accessions)
write.table(
  data.frame(Accession=accessions, GeneName=genes, row.names = NULL),
  'Datasets/Proteomics/gene-names.tsv',
  sep='\t',
  col.names = TRUE,
  quote=FALSE
)
