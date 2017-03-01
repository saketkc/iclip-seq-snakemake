#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
import seaborn as sns
sns.set_style('whitegrid')

import click
import pandas as pd

@click.command()
@click.option('--infile', help='Mouse fiel', required=True)
@click.option('--outfile', help='Output file', required=True)
@click.option('--title', help='Plot title',required= True)

def plot_utr_distribution(infile, outfile, title):
    df = pd.read_table(infile, names=['chrom', 'start', 'stop', 'name', 'score', 'strand'])
    title += ' Total peaks: {}'.format(len(df.index))
    ax = df['name'].value_counts().plot(kind='bar', rot=0, color='black').set_title(title)
    fig = ax.get_figure()
    fig.savefig(outfile)

if __name__ == '__main__':
    plot_utr_distribution()
