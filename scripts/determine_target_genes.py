from ginterval import GIntervalTree
import click
import os
import pandas
from pandas import Series

def target_gene(chrom, start, stop, tree, name_type='name'):
    hits = []
    for interval in tree[chrom].search(start, stop):
        hits.append(interval.data[name_type])
    hits = set(hits)
    hit_names = ','.join(hits)
    return hit_names


@click.command()
@click.option('--bed', help='Input BED file', required=True)
@click.option('--annotation', help='Input annotation file', required=True)
@click.option('--annotationtype', help='gencode/ref_gene/known_gene', default='ref_gene')
@click.option('--outbed', help='Output BED file', required=True)

def determine_target_genes(bed, annotation,
                           annotationtype, outbed):

    columns = ['chrom', 'start', 'stop', 'name', 'score', 'strand']#, 'pvalue']
    df = pandas.read_table(bed, names=columns, header=None)
    annotation_tree = GIntervalTree().load_bed(annotation)#, 'tx', annotationtype)
    hits = lambda x: target_gene(x.chrom, x.start, x.stop, annotation_tree)
    df['name'] = df.apply(hits, axis=1)
    df['name'] = df['name'].replace('', 'XYZNOHITSXYZ')
    df.to_csv(outbed, sep='\t', header=False, index=False)
    all_genes = set(df['name'].str.split(',').apply(Series, 1).stack())
    filename, ext = os.path.splitext(outbed)
    filename = filename+'.genelist.ensembl'
    with open(filename, 'w') as f:
        f.write('\n'.join(all_genes))
    hits = lambda x: target_gene(x.chrom, x.start, x.stop, annotation_tree)#, 'name2')
    df['name'] = df.apply(hits, axis=1)
    df['name'] = df['name'].replace('', 'INTERGENIC')
    df.to_csv(outbed, sep='\t', header=False, index=False)
    all_genes = set(df['name'].str.split(',').apply(Series, 1).stack())
    filename, ext = os.path.splitext(outbed)
    filename = filename+'.genelist.genename'
    with open(filename, 'w') as f:
        f.write('\n'.join(all_genes))


if __name__ == '__main__':
    determine_target_genes()


