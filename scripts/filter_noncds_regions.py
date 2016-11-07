#!/usr/bin/env python
"""
Filter noncds or regions with no gene names
for input to GREAT
"""
import click
import pandas


@click.command()
@click.option('--bed', help='Input BED file', required=True)
@click.option('--outbed', help='Output BED file', required=True)
@click.option('--intergenic', is_flag=True, help='filter intergenic regions')
@click.option('--nogenehits', is_flag=True, help='filter XYZNOHITSZYZ')
def filter_noncds_regions(bed, outbed, intergenic, nogenehits):
    columns = ['chrom', 'start', 'stop', 'name', 'score', 'strand']#, 'pvalue']
    df = pandas.read_table(bed, names=columns, header=None)
    if intergenic:
        df_filtered = df[df['name']=='intergenic']
    else:
        df_filtered = df[df['name']=='XYZNOHITSXYZ']
    df_filtered['score'] = 1000
    df_filtered.to_csv(outbed, sep='\t', header=False, index=False)

if __name__ == '__main__':
    filter_noncds_regions()
