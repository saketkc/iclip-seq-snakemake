#!/usr/bin/env python
"""
collapse_barcode_bam.py
Modified from https://github.com/YeoLab/gscripts/blob/master/gscripts/clipseq/barcode_collapse.py
MIT License

Copyright (c) 2016, Saket Choudhary
Copyright (c) 2014, Gabriel Pratt, Olga Botvinnik, Michael Lovci, Patrick Liu, Gene Yeo

"""
from collections import Counter, OrderedDict
import csv
import os
import click
import pysam
import sys
import ntpath


def path_leaf(path):
    head, tail = ntpath.split(path)
    return head, tail

def trim_fileext(filepath):
    """
    Return file path sans the extension and the extension
    """
    return os.path.splitext(filepath)[0]

def hamming(str1, str2):
    """Hamming distance between two strings"""
    return sum(a!=b and not( a=='N' or b=='N' ) for a,b in zip(str1, str2))

def collapse_reads(reads, outbam, barcode_header_mapping, alternate_barcode_mapping=None, max_ham_distance=0):
    barcode_counter = Counter()
    removed_barcode_counter = Counter()
    unique_reads = 0
    for read in reads:
        query_name = read.query_name

        ##NOTE: This is a hack to get it working
        ## STAR doesn't seem to output '/1' for paired end data
        ##TODO: Check if this works for paired end data
        try:
            current_barcode = barcode_header_mapping[query_name]
        except KeyError:
            try:
                current_barcode = barcode_header_mapping[query_name+'/1']
            except KeyError:
                ##TODO For some arbitrary reason the bam contains the read names
                ## as 'SRRXXXX.1232123' missing the flowchannel information
                ## I am not sure why that happens, but the identifiers of these types
                ## seem to be themselves unique so we can rely on a fuzzy match
                ## which is actually not so fuzzy since we just check if the stored barcode information
                ## starts with this string
                ## Fuzzy search
                """
                lookup_count = 0
                lookup_key = None
                for key in barcode_header_mapping.keys():
                    if key.startswith(query_name+' '):
                        lookup_count+=1
                        lookup_key = key
                if lookup_count == 1:
                    current_barcode = barcode_header_mapping[lookup_key]
                else:
                    sys.stderr.write('Failed to lookup barcode for: {}\n'.format(query_name))
                    sys.stderr.write('Barcode mapping: {}\n'.format('\n'.join(barcode_header_mapping.keys())))
                    sys.exit(1)
                """
                current_barcode = alternate_barcode_mapping[query_name]

        if current_barcode in barcode_counter:
            removed_barcode_counter[current_barcode] += 1
        else:
            write_read = True
            ## Is the current barcode already an existing barcode?
            ## Then don't output this read!
            for existing_barcode in barcode_counter:
                if hamming(current_barcode, existing_barcode) <= max_ham_distance:
                    write_read = False
            if write_read:
                unique_reads+=1
                outbam.write(read)
        barcode_counter[current_barcode] += 1
    return barcode_counter, removed_barcode_counter, removed_barcode_counter

def is_bam_sorted(inbam):
    """Check if input is sorted bam by coordinate"""
    if not inbam.header['HD']['SO'] == 'coordinate':
        return False
    return True

def extract_five_prime_pileup(inbam, chrom, position):
    """
    Given a set of reads, at a certain
    (chrom, startposition) it outputs something like this:
        <chr> <star_pos> <Number of reads with 5' end mapping to + strand at (chr, star_pos)> <Number of reads with 5' end mapping to - strand at <chrom, start_pos>
    """
    positive_mappings = 0
    negative_mappings = 0

    for read in inbam.fetch(chrom, position):
        if read.is_reverse:
            negative_mappings +=1
        else:
            positive_mappings +=1
    return (positive_mappings, negative_mappings)


"""
def write_fpp_file(inbam, parsed_positions, outfpp):
    with open(outfpp, 'w') as out_f:
        for chrom, pos in parsed_positions:
            pos_mappings, neg_mappings =
            out_f.write('{} {} {} {}\n'.format())
"""

def write_fpp_record(fpp_file, record):
    fpp_file.write('{}\n'.format(('\t').join(str(x) for x in record)))

@click.command()
@click.option('--inbam', help='Input BAM', required=True)
@click.option('--outbam', help='Output BAM', required=True)
@click.option('--barcodes', help='Two column barcode map file', required=True)
#@click.option('--barcodelength', help='Barcode length', type=int, required=True)



def collapse_bam_barcodes(inbam, outbam, barcodes):
    """
    barcode_map = None
    with open(barcodes, 'r') as barcode_f:
        reader = csv.reader(barcode_f, delimiter='\t')
        barcode_map = {row[0]: row[1] for row in rows}
    """
    pysam.index(inbam)
    inbam = pysam.Samfile(inbam, 'rb')
    """
    if not is_bam_sorted(inbam):
        raise RuntimeError('{} is not sorted. Use "samtools sort "'.format(inbam.filename))
    """
    outbam = pysam.Samfile(outbam, 'wb', template=inbam)
    outfpp = trim_fileext(outbam.filename).decode() + '.fpp'
    previous_position = None
    previous_chrom = None
    outfpp_fp = open(outfpp, 'w')

    barcode_header_mapping = {}
    with open(barcodes) as barcodes_fh:
        reader = csv.reader(barcodes_fh, delimiter='\t')
        barcode_header_mapping = dict([(row[0], row[1]) for row in reader])

    alternate_barcode_mapping = {k.split(' ')[0]:v for k,v in barcode_header_mapping.items()}
    positive_mappings = OrderedDict()
    negative_mappings = OrderedDict()
    all_barcodes_counter = Counter()
    removed_barcodes_counter = Counter()
    previous_position = None
    fpp_recording = []
    pos_counts = 0
    neg_counts = 0

    for i, read in enumerate(inbam.fetch()):
        chrom = read.rname
        ## Is the read mapping to negative strand?
        ## Where is the 5' end?
        if read.is_reverse:
            strand = '-'
        else:
            strand = '+'

        five_prime_loc = read.positions[0] if strand == '+' else read.positions[-1]
        current_position = read.positions[0]

        five_prime_position = (chrom, five_prime_loc)
        if previous_position != current_position:
            ## Delete all previous positive mappings
            #print((previous_position, current_position, chrom, previous_chrom))
            ## Delete all previous negative mappings
            for (t_chrom, t_pos), reads in list(negative_mappings.items()):
                if t_pos >= current_position:
                    break
                all_barcodes, removed_barcodes, counts = collapse_reads(reads, outbam, barcode_header_mapping, alternate_barcode_mapping)
                #neg_counts+=counts
                #all_barcodes_counter += all_barcodes
                #removed_barcodes_counter += removed_barcodes
                del negative_mappings[(t_chrom, t_pos)]
            if previous_chrom != chrom:
                for (t_chrom, t_pos), reads in list(negative_mappings.items()):
                    all_barcodes, removed_barcodes, counts = collapse_reads(reads, outbam, barcode_header_mapping, alternate_barcode_mapping)
                    #neg_counts+=counts
                    #all_barcodes_counter += all_barcodes
                    #removed_barcodes_counter += removed_barcodes
                    del negative_mappings[(t_chrom, t_pos)]
                assert len(negative_mappings) == 0

            for (t_chrom, t_pos), reads in list(positive_mappings.items()):
                all_barcodes, removed_barcodes, counts = collapse_reads(reads, outbam, barcode_header_mapping, alternate_barcode_mapping)
                #pos_counts+= counts
                #all_barcodes_counter += all_barcodes
                #removed_barcodes_counter += removed_barcodes
                del positive_mappings[(t_chrom, t_pos)]

        if strand == '-':
            if five_prime_position not in list(negative_mappings.keys()):
                negative_mappings[five_prime_position] = []
            negative_mappings[five_prime_position].append(read)
        if strand == '+':
            if five_prime_position not in list(positive_mappings.keys()):
                positive_mappings[five_prime_position] = []
            positive_mappings[five_prime_position].append(read)
        previous_position = current_position
        previous_chrom = chrom
        #record = [previous_chrom, previous_position, pos_counts, neg_counts]
        #write_fpp_record(outfpp_fp, record)


    for (t_chrom, t_pos), reads in list(positive_mappings.items()):
        all_barcodes, removed_barcodes, counts = collapse_reads(reads, outbam, barcode_header_mapping, alternate_barcode_mapping)
        #pos_counts+=counts
        #all_barcodes_counter += all_barcodes
        #removed_barcodes_counter += removed_barcodes
        del positive_mappings[(t_chrom, t_pos)]
    ## Delete all previous negative mappings
    for (t_chrom, t_pos), reads in list(negative_mappings.items()):
        all_barcodes, removed_barcodes, counts = collapse_reads(reads, outbam, barcode_header_mapping, alternate_barcode_mapping)
        #neg_counts+=counts
        #all_barcodes_counter += all_barcodes
        #removed_barcodes_counter += removed_barcodes
        del negative_mappings[(t_chrom, t_pos)]
    neg_mapping_counts = len(negative_mappings)
    pos_mapping_counts = len(positive_mappings)
    record = [chrom, current_position, pos_mapping_counts, neg_mapping_counts]
    write_fpp_record(outfpp_fp, record)
    inbam.close()
    outbam.close()

    outfpp_fp.close()



if __name__ == '__main__':
    collapse_bam_barcodes()
