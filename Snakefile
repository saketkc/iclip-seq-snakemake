include:
    #'config_cluster.py'
    'config_human.py'

workdir: ANALYSIS_DIR

rule all:
    input:
        expand('preprocessed/{sample}_btrimmed.fq', sample=SAMPLES),
        expand('mapped/tracks/{sample}.bigwig', sample=SAMPLES),
        expand('mapped/beds/{sample}.bed', sample=SAMPLES),
        expand('mapped/bams/{sample}.sorted.bam', sample=SAMPLES),
        expand('mapped/peaks/{sample}.piranha_bam_05.bed', sample=SAMPLES),
        expand('mapped/peaks/{sample}.peaks_05.bed', sample=SAMPLES),
        expand('mapped/annotated_beds/{sample}.annotated.bed', sample=SAMPLES),
        'mapped/annotated_beds/'+GENOME_BUILD+'_union.annotated.bed',
        'mapped/peaks/'+GENOME_BUILD+'_union.peaks.bed',
        'mapped/annotated_peaks/'+GENOME_BUILD+'_intersected.peaks.targetgenes.bed',
        'mapped/peaks/'+GENOME_BUILD+'_intersected.peaks.strandspecific.bed',
        'mapped/annotated_peaks/'+GENOME_BUILD+'_intersected.peaks.types.lengths.bed',
        'mapped/plots/'+GENOME_BUILD+'_all_hist.png',
        'mapped/plots/'+GENOME_BUILD+'_all_utr_regions.png',
        ## UNCOMMENT these if you want to do other specie analysis
        #'mapped/annotated_peaks/'+GENOME_BUILD+'_union.peaks.scored.bed',
        #'mapped/annotated_peaks/'+LIFT_PREFIX+'_union_lifted.peaks.bed',
        #'mapped/annotated_peaks/'+LIFT_PREFIX+'_union_unmapped.peaks.bed',
        #'mapped/peaks/'+LIFT_PREFIX+'_unmapped.peaks.bed',
        #'mapped/peaks/'+LIFT_PREFIX+'_lifted.peaks.bed',
        #'mapped/plots/'+LIFT_PREFIX+'_lifted_hist.png',
        #'mapped/plots/'+LIFT_PREFIX+'_unmapped_hist.png',
        #'mapped/plots/'+LIFT_PREFIX+'_lifted_utr_regions.png',
        #'mapped/plots/'+GENOME_BUILD+'_specific_utr_regions.png',
        #'mapped/annotated_beds/'+LIFT_PREFIX+'_union.annotated.mapped.bed',
        #'mapped/annotated_beds_merged/'+LIFT_PREFIX+'_union.annotated.mapped.bed',
        #unmapped = 'mapped/annotated_beds_merged/'+LIFT_PREFIX+'_union.annotated.unmapped.bed',

rule perform_qc:
    input: expand('{rawdata_dir}/{specie}/{sample}.fq', specie=SPECIES, rawdata_dir=RAWDATA_DIR, sample=SAMPLES)
    params:
        out_dir = 'qc'
    output:
        html = 'qc/{sample}_fastqc.html',
        zip = 'qc/{sample}_fastqc.zip',
    shell:
        r'''
            fastqc -o {params.out_dir} -f fastq {input}
        '''

rule trim_barcodes:
    input: expand('{rawdata_dir}/{specie}/{sample}.fq', specie=SPECIES, rawdata_dir=RAWDATA_DIR, sample=SAMPLES)
    output:
        fastq = 'preprocessed/{sample}_btrimmed.fq',
        barcodes = 'preprocessed/{sample}.barcodes.txt'
    params:
        out_dir = 'preprocessed'
    shell:
        r'''
        mkdir -p {params.out_dir}/barcodes && export LC_ALL=en_US.UTF-8 && python {SRC_DIR}/collapse_barcode_fastq.py -l {BARCODE_LENGTH} -o {params.out_dir} {input}
        '''

rule trim_adapters:
    input: 'preprocessed/{sample}_btrimmed.fq'
    output: 'preprocessed/{sample}_btrimmed_trimmed.fq'
    shell:
        r'''
         trim_galore -o preprocessed {input}
        '''

rule map_star:
    input: 'preprocessed/{sample}_btrimmed_trimmed.fq'
    output: 'mapped/bams/{sample}.sorted.bam'
    params:
        prefix = 'mapped/bams/{sample}',
        unmapped = 'unmapped/fastq/{sample}',
        starlogs = 'mapped/starlogs'
    threads: 16
    shell:
        r'''
        STAR --runThreadN {threads} --genomeDir {STAR_INDEX}\
        --alignEndsType EndToEnd\
        --outFileNamePrefix {params.prefix}\
        --outSAMtype BAM SortedByCoordinate\
        --outReadsUnmapped {params.unmapped}\
        --readFilesIn {input} && mv {params.prefix}Aligned.sortedByCoord.out.bam {output} && mkdir -p {params.starlogs} && mv {params.prefix}Log.final.out {params.prefix}Log.out {params.prefix}Log.progress.out {params.starlogs}
        '''

rule collapse_barcode_by_header:
    input:
        bams = 'mapped/bams/{sample}.sorted.bam',
        barcodes = 'preprocessed/{sample}.barcodes.txt',
    output:
        usort = 'mapped/bams/collapsed_unsorted/{sample}.bam',
        sort = 'mapped/bams/collapsed/{sample}.bam'
    shell:
        r'''
        export LC_ALL=en_US.UTF-8 && python {SRC_DIR}/collapse_barcode_by_header.py --barcodes {input.barcodes}\
        --inbam {input.bams} --outbam {output.usort} && samtools sort {output.usort} -o {output.sort}

        '''

rule bam_to_bed:
    input: 'mapped/bams/collapsed/{sample}.bam'
    output: 'mapped/beds/{sample}.bed'
    shell:
        r'''
        bamToBed -i {input} | sort -k1,1 -k3,3n -k2,2n -k4,4 > {output}
        '''

rule call_peaks_on_bam:
    input: 'mapped/bams/collapsed/{sample}.bam'
    output:'mapped/peaks/{sample}.piranha_bam_05.bed'
    shell:
        ##TODO What is sorting rule for bam??
        r'''
        Piranha -u 5 -s -b 1 {input} > {output}
        '''

rule create_wigs:
    input:
        bed = 'mapped/peaks/{sample}.piranha_bam_05.bed',
        bam = 'mapped/bams/collapsed/{sample}.bam'
    params:
        prefix = '{sample}',
        outdir = 'mapped/tracks'
    output:
        bedgraph = 'mapped/tracks/{sample}.bedGraph',
        bigwig = 'mapped/tracks/{sample}.bigwig'
    shell:
        r'''
        export LC_ALL=en_US.UTF-8 && python {SRC_DIR}/create_wigs.py\
                                        --bam {input.bam}\
                                        --outdir {params.outdir}\
                                        --genome {GENOME_SIZES}\
                                        --prefix {params.prefix}\
                                        --normalize\
                                        --strand
        '''

rule extract_peaks:
    input: 'mapped/peaks/{sample}.piranha_bam_05.bed'
    output: 'mapped/peaks/{sample}.peaks_05.bed'
    shell:
        r'''
        cut -f1,2,3,4,5,6 {input} > {output}
        '''

rule lift_to_other:
    input: 'mapped/peaks/{sample}.piranha_bam_05.bed'
    output:
        mapped = 'mapped/peaks/{sample}.piranha_bam_05.'+LIFT_PREFIX+'.bed',
        unmapped = 'mapped/peaks/{sample}.piranha_bam_05.'+LIFT_PREFIX+'.unmapped.bed',
    params:
        chain = expand('{genome_dir}/{genome_build}/liftover/{lift_prefix}.over.chain', genome_dir=GENOMES_DIR,
                       genome_build=GENOME_BUILD,
                       genome_build_other=GENOME_BUILD_OTHER,
                       lift_prefix=LIFT_PREFIX)
    shell:
        r'''
        liftOver {input} {params.chain} {output.mapped} {output.unmapped}
        '''

rule intersect_peaks:
    input: expand('mapped/peaks/{sample}.peaks_05.bed', sample=SAMPLES)
    output: 'mapped/peaks/'+GENOME_BUILD+'_intersected.peaks.bed'
    shell:
        # Not stranded
        r'''
        bedtools intersect -a {input[1]} -b {input[0]} {input[2]} > {output}
        '''

rule union_peaks:
    input: expand('mapped/peaks/{sample}.peaks_05.bed', sample=SAMPLES)
    output:
        #'mapped/peaks/'+GENOME_BUILD+'_union.peaks.uncollapsed.bed',
        'mapped/peaks/'+GENOME_BUILD+'_union.peaks.bed'
    shell:
        # Not stranded
        r'''
        cat {input[0]} {input[1]} {input[2]} | sort -k1,1 -k2,2n > {output}
        '''

rule annotate_union_peaks:
    input: 'mapped/peaks/'+GENOME_BUILD+'_union.peaks.bed'
    output: 'mapped/annotated_peaks/'+GENOME_BUILD+'_union.peaks.types.bed'
    params:
        gene_annotations = GENE_BED
    shell:
        r'''
        export LC_ALL=en_US.UTF-8 && python {SRC_DIR}/determine_target_genes.py --annotation {params.gene_annotations}\
        --bed {input}\
        --outbed {output}
        '''

rule unique_union_peaks:
    input: 'mapped/annotated_peaks/'+GENOME_BUILD+'_union.peaks.types.bed'
    output: 'mapped/annotated_peaks/'+GENOME_BUILD+'_union.peaks.unique.bed'
    shell:
        r'''
        awk '{{print $1"\t"$2"\t"$3"\t"$1"_"$2"_"$3"_"$4"\t"".""\t"$6}}' {input} | sort -k1,1 -k2,2n | uniq -u > {output}
        '''

rule count_union_replicates:
    input:
        'mapped/annotated_peaks/'+GENOME_BUILD+'_union.peaks.unique.bed',
        expand('mapped/peaks/{sample}.peaks_05.bed', sample=SAMPLES)
    output: 'mapped/annotated_peaks/'+GENOME_BUILD+'_union.peaks.scored.bed'
    shell:
        r'''
        python {SRC_DIR}/count_replicates.py --master {input[0]} --outbed {output} {input[1]} {input[2]} {input[3]}
        '''

rule lift_union_peaks:
    input: 'mapped/annotated_peaks/'+GENOME_BUILD+'_union.peaks.scored.bed'
    output:
        mapped = 'mapped/annotated_peaks/'+LIFT_PREFIX+'_union_lifted.peaks.bed',
        unmapped = 'mapped/annotated_peaks/'+LIFT_PREFIX+'_union_unmapped.peaks.bed',
    params:
        chain = GENOMES_DIR+'/'+GENOME_BUILD+'/liftover/'+LIFT_PREFIX+'.over.chain',
        unmapped = 'mapped/annotated_peaks/'+LIFT_PREFIX+'_union_unmapped_all.peaks.bed',
        unmapped_info = 'mapped/annotated_peaks/'+LIFT_PREFIX+'_union_unmapped_info.peaks.bed',
    shell:
        r'''
        liftOver {input} {params.chain} {output.mapped} {params.unmapped} &&\
        awk 'NR % 2 == 0 {{print;}}' {params.unmapped} > {output.unmapped} &&\
        awk 'NR % 2 == 1 {{print;}}' {params.unmapped} > {params.unmapped_info}
        '''




rule intersect_peaks_strandspecific:
    input: expand('mapped/peaks/{sample}.peaks_05.bed', sample=SAMPLES)
    output: 'mapped/peaks/'+GENOME_BUILD+'_intersected.peaks.strandspecific.bed'
    shell:
        # Not stranded
        r'''
        bedtools intersect -s -a {input[1]} -b {input[0]} {input[2]} > {output}
        '''

rule create_wigs_interesected:
    input:
        bed = 'mapped/peaks/'+GENOME_BUILD+'_intersected.peaks.bed'
    params:
        prefix = GENOME_BUILD+'_intersected.strandspecific',
        outdir = 'mapped/tracks'
    output:
        bedgraph = 'mapped/tracks/'+ GENOME_BUILD+'_intersected.bedGraph',
        bigwig = 'mapped/tracks/'+ GENOME_BUILD+'_intersected.bigwig'
    shell:
        r'''
        export LC_ALL=en_US.UTF-8 && python {SRC_DIR}/create_wigs.py\
                                        --bed {input.bed}\
                                        --outdir {params.outdir}\
                                        --genome {GENOME_SIZES}\
                                        --prefix {params.prefix}\
                                        --normalize\
                                        --strand
        '''

rule create_wigs_interesected_ss:
    input:
        bed = 'mapped/peaks/'+GENOME_BUILD+'_intersected.peaks.strandspecific.bed'
    params:
        prefix = GENOME_BUILD+'_intersected.strandspecific',
        outdir = 'mapped/tracks'
    output:
        bedgraph = 'mapped/tracks/'+ GENOME_BUILD+'_intersected.strandspecific.bedGraph',
        bigwig = 'mapped/tracks/'+ GENOME_BUILD+'_intersected.strandspecific.bigwig'
    shell:
        r'''
        export LC_ALL=en_US.UTF-8 && python {SRC_DIR}/create_wigs.py\
                                        --bed {input.bed}\
                                        --outdir {params.outdir}\
                                        --genome {GENOME_SIZES}\
                                        --prefix {params.prefix}\
                                        --normalize\
                                        --strand
        '''

rule lift_intersected_peaks:
    input: 'mapped/peaks/'+GENOME_BUILD+'_intersected.peaks.bed'
    output:
        mapped = 'mapped/peaks/'+LIFT_PREFIX+'_lifted.peaks.bed',
        unmapped = 'mapped/peaks/'+LIFT_PREFIX+'_unmapped.peaks.bed',
    params:
        chain = GENOMES_DIR+'/'+GENOME_BUILD+'/liftover/'+LIFT_PREFIX+'.over.chain',
        unmapped = 'mapped/peaks/'+LIFT_PREFIX+'_unmapped_all.peaks.bed',
        unmapped_info = 'mapped/peaks/'+LIFT_PREFIX+'_unmapped_info.peaks.bed',
    shell:
        r'''
        liftOver {input} {params.chain} {output.mapped} {params.unmapped} &&\
        awk 'NR % 2 == 0 {{print;}}' {params.unmapped} > {output.unmapped} &&\
        awk 'NR % 2 == 1 {{print;}}' {params.unmapped} > {params.unmapped_info}
        '''

rule determine_target_genes_all:
    input: 'mapped/peaks/'+GENOME_BUILD+'_intersected.peaks.bed'
    output:
        'mapped/annotated_peaks/'+GENOME_BUILD+'_intersected.peaks.targetgenes.bed',
        'mapped/annotated_peaks/'+GENOME_BUILD+'_intersected.peaks.targetgenes.genelist.ensembl',
        'mapped/annotated_peaks/'+GENOME_BUILD+'_intersected.peaks.targetgenes.genelist.genename',
    params:
        gene_annotations = GENE_BED
    shell:
        r'''
        export LC_ALL=en_US.UTF-8 && python {SRC_DIR}/determine_target_genes.py --annotation {params.gene_annotations}\
        --bed {input}\
        --outbed {output[0]}
        '''

rule pick_nogenehit_coordinates:
    input: 'mapped/annotated_peaks/'+GENOME_BUILD+'_intersected.peaks.targetgenes.bed'
    output: 'mapped/annotated_peaks/'+GENOME_BUILD+'_intersected.peaks.nogenehit.great.bed'
    shell:
        r'''
        export LC_ALL=en_US.UTF-8 && python {SRC_DIR}/filter_noncds_regions.py --bed {input}\
        --outbed {output} --nogenehits

        '''

rule pick_intergenic_coordinates:
    input: 'mapped/annotated_peaks/'+GENOME_BUILD+'_intersected.peaks.targetgenes.bed'
    output: 'mapped/annotated_peaks/'+GENOME_BUILD+'_intersected.peaks.intergenic.great.bed'
    shell:
        r'''
        export LC_ALL=en_US.UTF-8 && python {SRC_DIR}/filter_noncds_regions.py --bed {input}\
        --outbed {output} --intergenic

        '''

rule run_ensembl_orthology:
    input: 'mapped/annotated_peaks/'+GENOME_BUILD+'_intersected.peaks.targetgenes.genelist.ensembl'
    output: 'mapped/orthology/'+GENOME_BUILD+'_intersected.peaks.orthology'
    shell:
        r'''
        Rscript --vanilla {SRC_DIR}/run_orthology_{SPECIES}.R {input} {output}
        '''

rule determine_peak_type:
    input: 'mapped/peaks/'+GENOME_BUILD+'_intersected.peaks.bed'
    output: 'mapped/annotated_peaks/'+GENOME_BUILD+'_intersected.peaks.types.bed'
    params:
        cds = CDS_BED,
        utr5 = UTR5_BED,
        utr3 = UTR3_BED,
        intron = INTRON_BED

    shell:
        r'''
        export LC_ALL=en_US.UTF-8 && python {SRC_DIR}/determine_peak_type.py --cds {params.cds}\
        --utr5 {params.utr5}\
        --utr3 {params.utr3}\
        --intron {params.intron}\
        --bed {input}\
        --outbed {output}
        '''

rule determine_peak_type_lengths:
    input: 'mapped/peaks/'+GENOME_BUILD+'_intersected.peaks.bed'
    output: 'mapped/annotated_peaks/'+GENOME_BUILD+'_intersected.peaks.types.lengths.bed'
    params:
        cds = CDS_BED,
        utr5 = UTR5_BED,
        utr3 = UTR3_BED,
        intron = INTRON_BED

    shell:
        r'''
        export LC_ALL=en_US.UTF-8 && python {SRC_DIR}/determine_peak_type.py --cds {params.cds}\
        --utr5 {params.utr5}\
        --utr3 {params.utr3}\
        --intron {params.intron}\
        --bed {input}\
        --outbed {output}\
        --length
        '''

rule plot_utr_distribution_all:
    input: 'mapped/annotated_peaks/'+GENOME_BUILD+'_intersected.peaks.types.bed'
    output: 'mapped/plots/'+GENOME_BUILD+'_all_utr_regions.png'
    shell:
        r'''
        export LC_ALL=en_US.UTF-8 && python {SRC_DIR}/plot_utr_distribution.py --infile {input} --outfile {output} --title {GENOME_BUILD}_all
        '''

rule plot_peak_length_hist_all:
    input: 'mapped/annotated_peaks/'+GENOME_BUILD+'_intersected.peaks.types.bed'
    output: 'mapped/plots/'+GENOME_BUILD+'_all_hist.png'
    params:
        prefix = 'mapped/plots/'+GENOME_BUILD+'_all_hist'
    shell:
        r'''
        export LC_ALL=en_US.UTF-8 && python {SRC_DIR}/plot_peak_hist.py --inbed {input} --outprefix {params.prefix} --title {GENOME_BUILD}_all
        '''
rule determine_peak_type_lifted:
    input: 'mapped/peaks/'+LIFT_PREFIX+'_lifted.peaks.bed',
    output: 'mapped/annotated_peaks/'+LIFT_PREFIX+'_lifted.peaks.types.bed'
    params:
        cds = CDS_BED_OTHER,
        utr5 = UTR5_BED_OTHER,
        utr3 = UTR3_BED_OTHER,
        intron = INTRON_BED_OTHER

    shell:
        r'''
        export LC_ALL=en_US.UTF-8 && python {SRC_DIR}/determine_peak_type.py --cds {params.cds}\
        --utr5 {params.utr5}\
        --utr3 {params.utr3}\
        --intron {params.intron}\
        --bed {input}\
        --outbed {output}
        '''


rule plot_utr_distribution_lifted:
    input: 'mapped/annotated_peaks/'+LIFT_PREFIX+'_lifted.peaks.types.bed'
    output: 'mapped/plots/'+LIFT_PREFIX+'_lifted_utr_regions.png'
    shell:
        r'''
        export LC_ALL=en_US.UTF-8 && python {SRC_DIR}/plot_utr_distribution.py --infile {input} --outfile {output} --title {LIFT_PREFIX}
        '''

rule plot_peak_length_hist_lifted:
    input: 'mapped/annotated_peaks/'+LIFT_PREFIX+'_lifted.peaks.types.bed'
    output: 'mapped/plots/'+LIFT_PREFIX+'_lifted_hist.png'
    params:
        prefix = 'mapped/plots/'+LIFT_PREFIX+'_lifted_hist'
    shell:
        r'''
        export LC_ALL=en_US.UTF-8 && python {SRC_DIR}/plot_peak_hist.py --inbed {input} --outprefix {params.prefix} --title {LIFT_PREFIX}
        '''

rule determine_target_genes_lifted:
    input: 'mapped/peaks/'+LIFT_PREFIX+'_lifted.peaks.bed',
    output:
        'mapped/annotated_peaks/'+LIFT_PREFIX+'_lifted.peaks.targetgenes.bed',
        'mapped/annotated_peaks/'+LIFT_PREFIX+'_lifted.peaks.targetgenes.genelist.ensembl',
        'mapped/annotated_peaks/'+LIFT_PREFIX+'_lifted.peaks.targetgenes.genelist.genename',
    params:
        gene_annotations = GENE_BED_OTHER
    shell:
        r'''
        export LC_ALL=en_US.UTF-8 && python {SRC_DIR}/determine_target_genes.py --annotation {params.gene_annotations}\
        --bed {input}\
        --outbed {output[0]}
        '''

rule determine_peak_type_specific:
    input: 'mapped/peaks/'+LIFT_PREFIX+'_unmapped.peaks.bed',
    output: 'mapped/annotated_peaks/'+LIFT_PREFIX+'_unmapped.peaks.types.bed'
    params:
        cds = CDS_BED,
        utr5 = UTR5_BED,
        utr3 = UTR3_BED,
        intron = INTRON_BED

    shell:
        r'''
        export LC_ALL=en_US.UTF-8 && python {SRC_DIR}/determine_peak_type.py --cds {params.cds}\
        --utr5 {params.utr5}\
        --utr3 {params.utr3}\
        --intron {params.intron}\
        --bed {input}\
        --outbed {output}
        '''

rule plot_utr_distribution_specific:
    input: 'mapped/annotated_peaks/'+LIFT_PREFIX+'_unmapped.peaks.types.bed'
    output: 'mapped/plots/'+GENOME_BUILD+'_specific_utr_regions.png'
    shell:
        r'''
        export LC_ALL=en_US.UTF-8 && python {SRC_DIR}/plot_utr_distribution.py --infile {input} --outfile {output} --title {GENOME_BUILD}_specific
        '''

rule plot_peak_length_hist_specific:
    input: 'mapped/annotated_peaks/'+LIFT_PREFIX+'_unmapped.peaks.types.bed'
    output: 'mapped/plots/'+LIFT_PREFIX+'_unmapped_hist.png'
    params:
        prefix = 'mapped/plots/'+LIFT_PREFIX+'_unmapped_hist'
    shell:
        r'''
        export LC_ALL=en_US.UTF-8 && python {SRC_DIR}/plot_peak_hist.py --inbed {input} --outprefix {params.prefix} --title {GENOME_BUILD}_specific
        '''

rule determine_target_genes_specific:
    input: 'mapped/peaks/'+LIFT_PREFIX+'_unmapped.peaks.bed',
    output:
        'mapped/annotated_peaks/'+LIFT_PREFIX+'_unmapped.peaks.targetgenes.bed',
        'mapped/annotated_peaks/'+LIFT_PREFIX+'_unmapped.peaks.targetgenes.genelist.ensembl',
        'mapped/annotated_peaks/'+LIFT_PREFIX+'_unmapped.peaks.targetgenes.genelist.genename',
    params:
        gene_annotations = GENE_BED
    shell:
        r'''
        export LC_ALL=en_US.UTF-8 && python {SRC_DIR}/determine_target_genes.py --annotation {params.gene_annotations}\
        --bed {input}\
        --outbed {output[0]}
        '''

rule annotate_all_beds:
    input: 'mapped/peaks/{sample}.peaks_05.bed'
    output: 'mapped/annotated_beds/{sample}.annotated.bed'
    shell:
        r'''export LC_ALL=en_US.UTF-8 && python {SRC_DIR}/annotate_by_gene_region.py --cds {CDS_BED} \
            --utr5 {UTR5_BED} \
            --utr3 {UTR3_BED} \
            --intron {INTRON_BED} \
            --mirna {MIRNA_BED} \
            --lincrna {LINCRNA_BED} \
            --gene {GENE_BED} \
            --genename {GENE_NAMES} \
            --bed {input} \
            --outbed {output}
            '''

rule annotate_union_bed:
    input: 'mapped/peaks/'+GENOME_BUILD+'_union.peaks.bed'
    output: 'mapped/annotated_beds/'+GENOME_BUILD+'_union.annotated.bed'
    shell:
        r'''export LC_ALL=en_US.UTF-8 && python {SRC_DIR}/annotate_by_gene_region.py --cds {CDS_BED} \
            --utr5 {UTR5_BED} \
            --utr3 {UTR3_BED} \
            --intron {INTRON_BED} \
            --mirna {MIRNA_BED} \
            --lincrna {LINCRNA_BED} \
            --gene {GENE_BED} \
            --genename {GENE_NAMES} \
            --bed {input} \
            --outbed {output}
            '''

rule liftover_annotated_union_bed:
    input: 'mapped/annotated_beds/'+GENOME_BUILD+'_union.annotated.bed'
    output:
        mapped = 'mapped/annotated_beds/'+LIFT_PREFIX+'_union.annotated.mapped.bed',
        unmapped = 'mapped/annotated_beds/'+LIFT_PREFIX+'_union.annotated.unmapped.bed',
    params:
        chain = GENOMES_DIR+'/'+GENOME_BUILD+'/liftover/'+LIFT_PREFIX+'.over.chain',
        unmapped = 'mapped/annotated_beds/'+LIFT_PREFIX+'_union.annotated.unmapped_withinfo.bed',
        unmapped_info = 'mapped/annotated_beds/'+LIFT_PREFIX+'_union.annotated.unmapped_info.bed',
    shell:
        r'''
        liftOver {input} {params.chain} {output.mapped} {params.unmapped} &&\
        awk 'NR % 2 == 0 {{print;}}' {params.unmapped} > {output.unmapped} &&\
        awk 'NR % 2 == 1 {{print;}}' {params.unmapped} > {params.unmapped_info}
        '''

rule merge_union_bed:
    input: 'mapped/annotated_beds/'+GENOME_BUILD+'_union.annotated.bed'
    output:
        'mapped/annotated_beds_merged/'+GENOME_BUILD+'_union.annotated.bed3',
        'mapped/annotated_beds_merged/'+GENOME_BUILD+'_union.annotated.bed'
    shell:
        r'''export LC_ALL=en_US.UTF-8 && bedtools merge -i {input} > {output[0]} && python {SRC_DIR}/annotate_by_gene_region.py --cds {CDS_BED} \
            --utr5 {UTR5_BED} \
            --utr3 {UTR3_BED} \
            --intron {INTRON_BED} \
            --mirna {MIRNA_BED} \
            --lincrna {LINCRNA_BED} \
            --gene {GENE_BED} \
            --genename {GENE_NAMES} \
            --bed {output[0]} \
            --outbed {output[1]} &&\
            bedSort {output[1]} {output[1]}
            '''

rule liftover_merged_union_bed:
    input: 'mapped/annotated_beds_merged/'+GENOME_BUILD+'_union.annotated.bed'
    output:
        mapped = 'mapped/annotated_beds_merged/'+LIFT_PREFIX+'_union.annotated.mapped.bed',
        unmapped = 'mapped/annotated_beds_merged/'+LIFT_PREFIX+'_union.annotated.unmapped.bed',
    params:
        chain = GENOMES_DIR+'/'+GENOME_BUILD+'/liftover/'+LIFT_PREFIX+'.over.chain',
        unmapped = 'mapped/annotated_beds/'+LIFT_PREFIX+'_union.annotated.unmapped_withinfo.bed',
        unmapped_info = 'mapped/annotated_beds/'+LIFT_PREFIX+'_union.annotated.unmapped_info.bed',
    shell:
        r'''
        liftOver {input} {params.chain} {output.mapped} {params.unmapped} &&\
        awk 'NR % 2 == 0 {{print;}}' {params.unmapped} > {output.unmapped} &&\
        awk 'NR % 2 == 1 {{print;}}' {params.unmapped} > {params.unmapped_info} &&\
        bedSort {output.mapped} {output.mapped} &&\
        bedSort {output.unmapped} {output.unmapped}
        '''

