from ginterval import GIntervalTree
import click
import pandas
import re

TREE_PRIORITY = ['cds', 'utr5', 'utr3', 'intron', 'miRNA', 'lincRNA']

def determine_region_type(chrom, start, stop, trees, gene_tree, genenames):
    region_type = None
    hits = []
    tree_index = None
    for index, tree in enumerate(trees):
        if tree[chrom].search(start, stop):
            region_type = TREE_PRIORITY[index]
            tree_index = index
            break
    ## None found so must be intergenic
    if not region_type:
        region_type = 'intergenic'
    if region_type!='intergenic':
        for interval in gene_tree[chrom].search(start, stop):

            #gene_name = genenames.loc[interval.data['name']]['name']
            ## Keep ENSEMBL Ids but remove the version info
            gene_name = re.sub(r'\.[0-9]+','', interval.data['name'])
            gene_type = genenames.loc[interval.data['name']]['type']
            hits.append((gene_name, gene_type))
    if len(hits)>1 and region_type!='intergenic':
        ## We don't want to keep antisense genes for annotations
        for index, hit in hits:
            if hit[1] == 'antisense':
                hits.pop(index)

    if len(hits)>1:
        #sys.stderr.write('Encountered hits: {}\n'.format(hits))
        hits = [hits[0]]
    if hits and region_type!='lincRNA':
        return '{}__{}__{}__{}__{}'.format(chrom, start, stop, hits[0][0], region_type)
    elif hits and region_type=='lincRNA':
        region_hits = trees[tree_index][chrom].search(start, stop)
        #assert len(region_hits)==1
        for hit in region_hits:
            region_name = hit.data['name']
            break
        region_name = re.sub(r'\.[0-9]+', '', region_name)
        return '{}__{}__{}__{}'.format(chrom, start, stop, region_name)
    return '{}__{}__{}__{}__{}'.format(chrom, start, stop, 'INTERGENIC' , region_type)

def determine_region_type_and_length(chrom, start, stop, trees, genenames):
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
@click.option('--mirna', help='miRNA annotations', required=True)
@click.option('--lincrna', help='lncRNA annotations', required=True)
@click.option('--gene', help='Input annotation file', required=True)
@click.option('--genename', help='Gene name mappings and type', required=True)
@click.option('--bed', help='Input BED file', required=True)
@click.option('--length', help='Include length values', is_flag=True)
@click.option('--outbed', help='Output BED file', required=True)

def determine_peak_region(cds, utr5, utr3, intron, mirna, lincrna, gene,
                          genename, bed, length, outbed):
    """
    Parameters
    ----------

    """
    columns = ['chrom', 'start', 'stop', 'name', 'score', 'strand']#, 'pvalue']
    df = pandas.read_table(bed, header=None)
    if len(df.columns)==6:
        df.columns = columns
    elif len(df.columns)==3:
        df.columns = columns[:3]
        df['name'] = '.'
        df['score'] = '.'
        df['strand'] = '.'
    else:
        raise RuntimeError('BED file should be 3 or 6 column')

    utr3_tree = GIntervalTree().load_bed(utr3)
    utr5_tree = GIntervalTree().load_bed(utr5)
    cds_tree = GIntervalTree().load_bed(cds)
    intron_tree = GIntervalTree().load_bed(intron)
    gene_tree = GIntervalTree().load_bed(gene)
    mirna_tree = GIntervalTree().load_bed(mirna)
    lincrna_tree = GIntervalTree().load_bed(lincrna)
    all_trees = [cds_tree, utr5_tree, utr3_tree, intron_tree, mirna_tree, lincrna_tree]

    genename = pandas.read_table(genename, names=['id', 'name', 'type'], index_col=0)

    if not length:
        hits = lambda x: determine_region_type(x.chrom, x.start, x.stop, all_trees, gene_tree, genename)
        df['name'] = df.apply(hits, axis=1)
    else:
        hits = lambda x: determine_region_type_and_length(x.chrom, x.start, x.stop, all_trees, gene_tree, genename)
        df['name'] = df.apply(hits, axis=1)

    #df['name'] = df.chrom + '__' + df.start.map(str) + '__' + df.stop.map(str) + '__' + df.namedi
    split_col = df['name'].str.split('__').str
    all_genes = split_col.get(3)
    split_col = df['name'].str.split('__').str
    all_regions = split_col.get(4)
    gene_counts = all_genes.value_counts()
    region_counts = all_regions.value_counts()

    gene_counts_df = pandas.DataFrame.from_dict(dict(gene_counts), orient='index')
    gene_names_df = pandas.DataFrame({'gene_names': list(dict(gene_counts).keys())})
    gene_counts_df.columns=['count']
    gene_counts_df.sort_values(by=['count'], ascending=[0], inplace=True)
    region_counts_df = pandas.DataFrame.from_dict(dict(region_counts), orient='index')
    region_counts_df.columns=['count']
    region_counts_df.sort_values(by=['count'], ascending=[0], inplace=True)
    gene_counts_df.to_csv(outbed+'.gene_counts', sep='\t', header=False)#, index=False)
    gene_names_df.to_csv(outbed+'.gene_names', sep='\t', header=False, index=False)
    region_counts_df.to_csv(outbed+'.region_counts', sep='\t', header=False)#, index=False)

    region_counts_df.to_csv()
    df.to_csv(outbed, sep='\t', header=False, index=False)

if __name__ == '__main__':
    determine_peak_region()
