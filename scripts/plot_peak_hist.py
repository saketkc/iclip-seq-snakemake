#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
import seaborn as sns
import click
import pandas as pd
import matplotlib.pylab as plt
sns.set_style('whitegrid')

BED_COL_NAMES = ['chrom', 'start', 'stop', 'name', 'signalValue' ,'strand',]

@click.command()
@click.option('--inbed', help='Input BED',required=True)
@click.option('--outprefix', help='Output prefix', required=True)
@click.option('--title', help='Title', required=True)
def plot_hist_peak(inbed, outprefix, title):
    df = pd.read_table(inbed, names=BED_COL_NAMES)
    df['peak_length'] = df['stop']-df['start']
    peak_lengths = df['peak_length'].tolist()
    X = []
    Y = []
    peak_lengths = df['peak_length'].tolist()
    total_peaks = len(peak_lengths)
    for length in set(peak_lengths):
        X.append(length)
        Y.append(peak_lengths.count(length))
    #fig = plt.figure()
    #plt.hist(df['peak_length'].tolist())
    ax = sns.barplot(X,Y)
    plt.legend()
    ax.set_xlabel('Peak length')
    ax.set_ylabel('Count')
    ax.set_title('{} | Total peaks = {}'.format(title, total_peaks))
    fig = ax.get_figure()
    fig.savefig('{}.png'.format(outprefix))

if __name__ == '__main__':
    plot_hist_peak()
