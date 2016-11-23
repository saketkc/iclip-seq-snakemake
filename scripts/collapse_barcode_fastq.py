#!/usr/bin/env python
"""Collapse barcodes from fastq
Usage:
    python remove_barcodes.py -l 5 <list_of_fq>
"""
import csv
import os
import click
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import ntpath
#from collections import Counter


def path_leaf(path):
    head, tail = ntpath.split(path)
    return head, tail

def split_path(filepath):
    """Return filename and extension of file"""
    return os.path.splitext(filepath)

def move_barcode_to_header(input_f, output_f, barcode_file_prefix, length):
    """Trim barcode and write to new fastq
    We rely on the fact that the idea experiment would have 5' end of each read
    starting with a random barcode. PCR duplicates will have identical barcodes
    and identical sequences. Other way round, you should expect different barcodes
    for sequences that are not barcoded.

    So the way to do this, is to move the barcode sequence to the readname header
    We do not want to play around with the format for read names, so we put it in the
    place where the index sequence is


    Arguments
    ---------
    input: string
        input fastq

    output: string
        output fastq

    barcode_file: string
        location to store collected barcodes

    length: int
        length of barcode
    """
    barcodes = []
    barcode_name_map = {}
    barcode_stats_file = barcode_file_prefix + '_stats.tsv'.format(length)
    barcode_map_file = barcode_file_prefix + '.barcodes.txt'.format(length)
    with open(input_f) as in_handle:
        with open(output_f, 'w') as out_handle, open(barcode_map_file, 'w') as map_handle:
            map_writer = csv.writer(map_handle, delimiter='\t')
            for title, seq, qual in FastqGeneralIterator(in_handle):
                barcode, trimmed_seq = seq[:length], seq[length:]
                #title_string, title_index = title.split('/')
                #title_start, title_index_seq = title_string.split('#')
                #new_title = '{}#{}/{}'.format(title_start, barcode, title_index)
                trimmed_qual = qual[length:]
                barcodes.append(barcode)
                map_writer.writerow([str(title), str(barcode)])
                out_handle.write('@{}\n{}\n+\n{}\n'.format(title, trimmed_seq, trimmed_qual))

@click.command()
@click.option('-l', '--length', help='Barcode length', required=True, type=int)
@click.option('-s', '--suffix', help='Suffix of output file', default='btrimmed', type=str)
@click.option('-o', '--outdir', help='Directory to output all files', type=str)
@click.argument('inputs', nargs=-1, required=True)

def remove_barcodes(inputs, length, suffix, outdir):
    """Remove barcodes for all input files
    Arguments
    ---------
    inputs: list
        List of fastqs
    length: int
        barcode length
    suffix: string
        Suffix for output fastqs
    outdir: string
        Directory where all files are written
    """
    if not outdir:
        outdir = os.getcwd()
    for input_f in inputs:
        if not os.path.isfile(input_f):
            raise RuntimeError('{} is not a valid file'.format(input_f))
        out_name, out_ext = split_path(input_f)
        output_f = os.path.join(outdir, path_leaf(out_name)[1]) + '_' + suffix + out_ext
        barcode_file_prefix = os.path.join(outdir, path_leaf(out_name)[1])
        move_barcode_to_header(input_f, output_f, barcode_file_prefix, length)

if __name__ == '__main__':
    remove_barcodes()
