#!/usr/bin/env python
"""Script to create bedGraph and bigwig normalized to
Reads per million mapped reads"""

__VERSION__ = '0.2'

import click
import os
import subprocess
import fileinput


@click.command()
@click.version_option(version=__VERSION__)
@click.option('--bam', help='Input bam')
@click.option('--bed', help='Input bam')
@click.option('--outdir', help='Output dir', required=True)
@click.option('--prefix', help='Prefix for output files', required=True)
@click.option('--trackname', help='Trackname')
@click.option('-g', '--genome', help='Genome fasta', required=True)
@click.option('--normalize', is_flag=True, help='Normalize to reads per million mapped reads')
@click.option('--strand', is_flag=True, help='Split by strand')
def create_wigs(bam, bed, prefix, outdir,
                trackname, genome, normalize, strand):
    bam = os.path.abspath(bam)
    samscale_cmd = 'samtools view -c -F 0x4 {}'.format(bam)
    if not trackname:
        trackname = prefix

    prefix = os.path.join(os.path.abspath(outdir), prefix)

    if normalize and bam:
        p = subprocess.Popen(samscale_cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()
        if stderr:
            raise ValueError('samtools error: %s' % stderr)
        readcount = float(stdout)
        scale = 1/ (readcount / 1e6)
        print('Read Count: {}'.format(readcount))
    else:
        scale = 1

    print('Scale: {}'.format(scale))
    if bam:
        bedgraph_cmd = 'bedtools genomecov -ibam {} -bg -scale {} > {}.bedGraph'.format(bam, scale, prefix)
        os.system(bedgraph_cmd)
    elif bed:
        bedgraph_cmd = 'bedtools genomecov -i {} -g {} -bg -scale {} > {}.bedGraph'.format(bam, genome, scale, prefix)
        os.system(bedgraph_cmd)

    bedsort_cmd = 'bedSort {0}.bedGraph {0}.bedGraph'.format(prefix)
    os.system(bedsort_cmd)
    bigwig_cmd = 'bedGraphToBigWig {0}.bedGraph {1} {0}.bigwig'.format(prefix, genome)
    os.system(bigwig_cmd)
    if strand:
        if bam:
            bedgraph_cmd = 'bedtools genomecov -ibam {} -bg -scale {} -strand + > {}.pos.bedGraph'.format(bam, scale, prefix)
            os.system(bedgraph_cmd)
            bedgraph_cmd = 'bedtools genomecov -ibam {} -bg -scale {} -strand - > {}.neg.bedGraph'.format(bam, scale, prefix)
            os.system(bedgraph_cmd)
        elif bed:
            bedgraph_cmd = 'bedtools genomecov -i {} -g {} -bg -scale {} -strand + > {}.pos.bedGraph'.format(bed, genome, scale, prefix)
            os.system(bedgraph_cmd)
            bedgraph_cmd = 'bedtools genomecov -i {} -g {} -bg -scale {} -strand - > {}.neg.bedGraph'.format(bed, genome, scale, prefix)
            os.system(bedgraph_cmd)

        bedsort_cmd = 'bedSort {0}.pos.bedGraph {0}.pos.bedGraph'.format(prefix)
        os.system(bedsort_cmd)
        bedsort_cmd = 'bedSort {0}.neg.bedGraph {0}.neg.bedGraph'.format(prefix)
        os.system(bedsort_cmd)

        bigwig_cmd = 'bedGraphToBigWig {0}.pos.bedGraph {1} {0}.pos.bigwig'.format(prefix, genome)
        os.system(bigwig_cmd)
        bigwig_cmd = 'bedGraphToBigWig {0}.neg.bedGraph {1} {0}.neg.bigwig'.format(prefix, genome)
        os.system(bigwig_cmd)

    bedsort_cmd = 'bedSort {0}.bedGraph {0}.bedGraph'.format(prefix)
    os.system(bedsort_cmd)
    bigwig_cmd = 'bedGraphToBigWig {0}.bedGraph {1} {0}.bigwig'.format(prefix, genome)
    os.system(bigwig_cmd)

    trackline='track name="{}" visibility="full" type=bedGraph color=200,100,0 altColor=0,100,200'.format(trackname)
    for line in fileinput.input(['{}.bedGraph'.format(prefix)], inplace=True):
        if fileinput.isfirstline():
            print(trackline)
        print(line),

if __name__ == '__main__':
    create_wigs()
