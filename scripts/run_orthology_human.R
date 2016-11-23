#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(biomaRt)
human <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl')

attributes <- c('ensembl_transcript_id',
                'ensembl_gene_id',
                'mmusculus_homolog_ensembl_gene',
                'mmusculus_homolog_perc_id_r1')

attributes <- c(attributes,
                'mmusculus_homolog_orthology_type',
                'mmusculus_homolog_subtype',
                'mmusculus_homolog_perc_id')#, 'hgnc_symbol', 'hgnc_id')

ortho.mouse <- getBM(attributes,
                    filters='with_mmusculus_homolog',
                    values=TRUE,
                    mart = human,
                    bmHeader=FALSE)

ortho.mouse.one2one <- subset(ortho.mouse, mmusculus_homolog_orthology_type == 'ortholog_one2one')
genelist <- readLines(args[1])
genelist <- gsub('exon:', '', genelist)
genelist <- gsub('\\..*', '', genelist)
ortho.mouse.subset <- subset(ortho.mouse.one2one, ensembl_gene_id %in% genelist | ensembl_transcript_id %in% genelist)
ortho.mouse.subset.genes <- ortho.mouse.subset[, c('ensembl_gene_id', 'mmusculus_homolog_ensembl_gene')]
ortho.mouse.subset.genes <- ortho.mouse.subset.genes[!duplicated(ortho.mouse.subset.genes), ]

write.table(ortho.mouse.subset, file=args[2], row.names=FALSE, sep='\t', quote=FALSE)
write.table(ortho.mouse.subset.genes, file=paste(args[2], 'genes', sep='.'), row.names=FALSE, sep='\t', quote=FALSE)
