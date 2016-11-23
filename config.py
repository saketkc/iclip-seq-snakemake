####################################################Input##################################################

## Specify barcode length
BARCODE_LENGTH = 9

## Specify genome build and species of species to which the fastqs belong
GENOME_BUILD = 'hg38'
SPECIES = 'human'

## Specify genome build of other species for liftover
GENOME_BUILD_OTHER = 'mm10'

## Parent location of raw .fq files
RAWDATA_DIR ='/home/cmb-06/as/skchoudh/data/HuR_Mouse_Human_liver/rawdata'

## List of .fq filess
## IMPORTANT: Assumes '.fq' suffix
SAMPLES = ['Liver-iCLIP-Index3_FCD2GYDACXX_L2_CHKSE13080180-Index3_1',
           'Liver-iCLIP-Index5_FCD2GYDACXX_L2_CHKSE13080180-Index5_1',
           'Liver-iCLIP-Index8_FCD2GYDACXX_L2_CHKSE13080180-Index8_1']

## Absolute path to the scripts directory
SRC_DIR = '/home/cmb-panasas2/skchoudh/github_projects/clip_seq_pipeline/scripts'

## Absolute path to root directory where results will be created
ANALYSIS_DIR = '/home/cmb-panasas2/skchoudh/HuR_results/human/clip_seq'
###########################################################################################################



##############################################GENOME specific##############################################
GENOMES_DIR='/home/cmb-panasas2/skchoudh/genomes'

GENOME_FASTA = GENOMES_DIR + '/' + GENOME_BUILD + '/' + 'fasta' + '/' + GENOME_BUILD + '.fa'

GENOME_SIZES = GENOMES_DIR + '/' + GENOME_BUILD + '/' + 'fasta' + '/' + GENOME_BUILD + '.sizes'

STAR_INDEX = GENOMES_DIR + '/' + GENOME_BUILD + '/' + 'star_annotated'

UTR5_BED = GENOMES_DIR + '/' + GENOME_BUILD + '/' + 'annotation' + '/' + 'gencode.v25.5UTRs.bed'

UTR3_BED = GENOMES_DIR + '/' + GENOME_BUILD + '/' + 'annotation' + '/' + 'gencode.v25.3UTRs.bed'

CDS_BED = GENOMES_DIR + '/' + GENOME_BUILD + '/' + 'annotation' + '/' + 'gencode.v25.CDS.bed'

INTRON_BED = GENOMES_DIR + '/' + GENOME_BUILD + '/' + 'annotation' + '/' + 'gencode.v25.introns.bed'

MIRNA_BED = GENOMES_DIR + '/' + GENOME_BUILD + '/' + 'miRNA' + '/'+ 'GrCh38.miRNA.bed6'

GENE_BED = GENOMES_DIR + '/' + GENOME_BUILD + '/' + 'annotation' + '/' + 'gencode.v25.genes.bed'

LIFTOVER_CHAIN = GENOMES_DIR + '/' + GENOME_BUILD + 'liftover' + '/' + 'hg38ToMm10.over.chain'

LINCRNA_BED = GENOMES_DIR + '/' + GENOME_BUILD + '/' + 'lncRNA' + '/'+ 'gencode.v25.long_noncoding_RNAs.named.bed'
GENE_BED = GENOMES_DIR + '/' + GENOME_BUILD + '/' + 'annotation' + '/' + 'gencode.v25.genes.bed'
GENE_NAMES = GENOMES_DIR + '/' + GENOME_BUILD + '/' + 'annotation' + '/' + GENOME_BUILD + '_gene_names.tsv'
###########################################################################################################


################################################OTHER GENOME################################################

UTR5_BED_OTHER = GENOMES_DIR + '/' + GENOME_BUILD_OTHER + '/' + 'annotation' + '/' + 'gencode.vM11.5UTRs.bed'

UTR3_BED_OTHER = GENOMES_DIR + '/' + GENOME_BUILD_OTHER + '/' + 'annotation' + '/' + 'gencode.vM11.3UTRs.bed'

CDS_BED_OTHER = GENOMES_DIR + '/' + GENOME_BUILD_OTHER + '/' + 'annotation' + '/' + 'gencode.vM11.CDS.bed'

INTRON_BED_OTHER = GENOMES_DIR + '/' + GENOME_BUILD_OTHER + '/' + 'annotation' + '/' + 'gencode.vM11.introns.bed'

MIRNA_BED_OTHER = GENOMES_DIR + '/' + GENOME_BUILD_OTHER + '/' + 'miRNA' + '/' + 'GrCh38.miRNA.bed'

GENE_BED_OTHER = GENOMES_DIR + '/' + GENOME_BUILD_OTHER + '/' + 'annotation' + '/' + 'gencode.vM11.genes.bed'

LIFTOVER_CHAIN_OTHER = GENOMES_DIR + '/' + GENOME_BUILD + 'liftover' + '/' + 'mm10ToHg38.over.chain'
###########################################################################################################


################################################DO NOT EDIT#################################################
LIFT_PREFIX = GENOME_BUILD + 'To' + GENOME_BUILD_OTHER
############################################################################################################
