#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
import seaborn as sns
import click
import pandas as pd
import matplotlib.pylab as plt
@click.command()
@click.option('--infile', help='Output file', required=True)
@click.option('--outprefix', help='Prefix')

def plot_size_wise(infile, outprefix):
    df = pd.read_table(infile, names=['chrom', 'start', 'stop', 'name', 'score', 'strand'])
    df = df[df['name']!='intergenic']
    #df['domain'] = None
    #df['tsizes'] = None
    #df[['domain', 'tsizes']] = df['name'].split('_')
    df['domain'], df['tsizes'] = zip(*df['name'].str.split('_').tolist())
    df['tsizes'] = df['tsizes'].astype(float)
    df['plengths'] = df['stop']-df['start']
    #df['peaks/size'] = df['plengths']/df['tsizes']
    df['peak_enrichment'] = df['score']
    #df = df[df['peaks/size']<=0.1]
    sns_plot = sns.violinplot(x='domain', y='peak_enrichment', data=df,hue='domain')#, order=['hg19_human_to_mouse', 'hg19_mouse_to_human', 'mm10_mouse_to_human', 'mm10_human_to_mouse'],hue="type")
    sns_plot.set_xticklabels(rotation=90,labels=[])
    sns_plot.set_title('Distribution of enrichment values among domains')

    fig = sns_plot.get_figure()
    fig.savefig('{}.png'.format(outprefix))

if __name__ == '__main__':
    plot_size_wise()

