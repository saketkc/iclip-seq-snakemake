from ginterval import GIntervalTree
import numpy as np
import click
import pandas

TREE_PRIORITY = ['cds', 'utr5', 'utr3', 'intron']
def determine_region_type(chrom, start, stop, trees):
    for index, tree in enumerate(trees):
        if tree[chrom].search(start, stop):
            return TREE_PRIORITY[index]
    ## None found so must be intergenic
    return 'intergenic'

def determine_region_type_and_length(chrom, start, stop, trees):
    for index, tree in enumerate(trees):
        search_tree = tree[chrom].search(start, stop)
        if search_tree:
            for interval in search_tree:
                length = interval.end-interval.begin
                return '{}_{}'.format(TREE_PRIORITY[index], length)
    ## None found so must be intergenic
    return 'intergenic'

@click.command()
@click.option('--cds', help='CDS regions', required=True)
@click.option('--utr5', help='5 UTR regions', required=True)
@click.option('--utr3', help='3 UTR regions', required=True)
@click.option('--intron', help='Intron regions', required=True)
@click.option('--bed', help='Input BED file', required=True)
@click.option('--length', help='Include length values', is_flag=True)
@click.option('--outbed', help='Output BED file', required=True)

def determine_peak_region(cds, utr5, utr3, intron, bed, length, outbed):
    """
    Parameters
    ----------

    """
    columns = ['chrom', 'start', 'stop', 'name', 'score', 'strand']#, 'pvalue']
    df = pandas.read_table(bed, names=columns, header=None)
    #df = df[np.isfinite(df['pvalue'])]

    utr3_tree = GIntervalTree().load_bed(utr3)
    utr5_tree = GIntervalTree().load_bed(utr5)
    cds_tree = GIntervalTree().load_bed(cds)
    intron_tree = GIntervalTree().load_bed(intron)
    all_trees = [cds_tree, utr5_tree, utr3_tree, intron_tree]

    if not length:
        hits = lambda x: determine_region_type(x.chrom, x.start, x.stop, all_trees)
        df['name'] = df.apply(hits, axis=1)
    else:
        hits = lambda x: determine_region_type_and_length(x.chrom, x.start, x.stop, all_trees, )
        df['name'] = df.apply(hits, axis=1)

    #df['name'] = df.chrom + '__' + df.start.map(str) + '__' + df.stop.map(str) + '__' + df.name
    df.to_csv(outbed, sep='\t', header=False, index=False)

if __name__ == '__main__':
    determine_peak_region()
