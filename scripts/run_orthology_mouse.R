#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(biomaRt)
mouse <- useMart('ensembl', dataset = 'mmusculus_gene_ensembl')

attributes <- c('ensembl_transcript_id',
                'ensembl_gene_id',
                'hsapiens_homolog_ensembl_gene',
                'hsapiens_homolog_perc_id_r1')

attributes <- c(attributes,
                'hsapiens_homolog_orthology_type',
                'hsapiens_homolog_subtype',
                'hsapiens_homolog_perc_id')#, 'hgnc_symbol', 'hgnc_id')

ortho.human <- getBM(attributes,
                    filters='with_hsapiens_homolog',
                    values=TRUE,
                    mart = mouse,
                    bmHeader=FALSE)

ortho.human.one2one <- subset(ortho.human, hsapiens_homolog_orthology_type == 'ortholog_one2one')
genelist <- readLines(args[1])
genelist <- gsub('exon:', '', genelist)
genelist <- gsub('\\..*', '', genelist)
ortho.human.subset <- subset(ortho.human.one2one, ensembl_transcript_id %in% genelist | ensembl_gene_id %in% genelist)
ortho.human.subset.genes <- ortho.human.subset[, c('ensembl_gene_id', 'hsapiens_homolog_ensembl_gene')]
ortho.human.subset.genes <- ortho.human.subset.genes[!duplicated(ortho.human.subset.genes), ]

write.table(ortho.human.subset, file=args[2], row.names=FALSE, sep='\t', quote=FALSE)
write.table(ortho.human.subset.genes, file=paste(args[2], 'genes', sep='.'), row.names=FALSE, sep='\t', quote=FALSE)
