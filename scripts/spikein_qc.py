#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import numpy as np
import scipy
import pandas as pd
from collections import OrderedDict
import matplotlib
from Bio import SeqIO
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import cm
import seaborn as sns

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description='Spikein-based QC of poly(A) tail calling.')
parser.add_argument(
    '-i', metavar='input', type=str, help="Input TSV.", required=True)
parser.add_argument(
    '-o', metavar='output', type=str, help="Output TSV.", required=True)
parser.add_argument(
    '-d', metavar='data_output', type=str, help="Data output TSV.", required=True)
parser.add_argument(
    '-m', metavar='medians', type=str, help="Medians TSV.", required=True)
parser.add_argument(
    '-r', metavar='output', type=str, help="Report PDF.", required=True)
parser.add_argument(
    '-s', metavar='spikein', type=str, help="Spike-in fasta.", required=True)


def _read_spikein_data(fas):
    """ Read in spike-in info."""
    ids = [r.id for r in SeqIO.parse(fas, "fasta")]
    tmp = {}
    for name in ids:
        tmp[name] = int(name.rsplit("_", 1)[1])
    df = pd.DataFrame(OrderedDict([('contig', list(tmp.keys())), ('TrueTail', list(tmp.values()))]))
    return df


def _make_boxplots(df, mdf, pages):
    """ Spike-in boxplots. """
    ax = sns.boxplot(x="TrueTail", y="polya_length", data=df)
    ax.set_title("Spike-in QC")
    for r in mdf.itertuples():
        if not np.isnan(r.polya_length):
            ax.text(r.Index, r.polya_length + 0.5, np.round(r.polya_length, 2), horizontalalignment='center', size='x-small', color='w', weight='semibold')
    plt.tight_layout()
    pages.savefig()
    plt.clf()
    ax = sns.boxplot(x="TrueTail", y="polya_length", data=df, showfliers=False)
    ax.set_title("Spike-in QC (with outliers removed)")
    for r in mdf.itertuples():
        if not np.isnan(r.polya_length):
            ax.text(r.Index, r.polya_length + 0.5, np.round(r.polya_length, 2), horizontalalignment='center', size='x-small', color='w', weight='semibold')
    plt.tight_layout()
    pages.savefig()
    plt.clf()


def _make_regplot(df, pages, title, jitter=0):
    ax = sns.regplot(x="TrueTail", y="polya_length", data=df, fit_reg=True, x_jitter=jitter)
    ax.set_title(title)
    if len(df.polya_length.values) == 0:
        return
    r, p = scipy.stats.pearsonr(df.TrueTail.values, df.polya_length.values)
    rtext = "Pearson r ={:.3f}, p-value={:.3f}".format(r, p)
    ax.text(1.0, 0.98, rtext, horizontalalignment='right', verticalalignment='center', transform=ax.transAxes)
    plt.tight_layout()
    pages.savefig()
    plt.clf()


if __name__ == '__main__':
    args = parser.parse_args()
    pages = PdfPages(args.r)
    sns.set_style("whitegrid")

    tails = pd.read_csv(args.i, sep="\t")

    sp = _read_spikein_data(args.s)

    df = pd.merge(sp, tails, how='left', on='contig')

    tails = tails[~tails.contig.isin(sp.contig.values)]
    tails.to_csv(args.d, sep="\t", index=False)
    del tails

    mdf = df.groupby('TrueTail').median()
    mdf = mdf.reset_index()[['TrueTail', 'polya_length']]
    _make_boxplots(df, mdf, pages)
    dfn = df.dropna()
    del df
    _make_regplot(mdf.copy().dropna(), pages, "Spike-in QC: regression (medians)")
    _make_regplot(dfn, pages, "Spike-in QC: regression", jitter=3)

    pages.close()

    mdf.to_csv(args.m, sep="\t", index=False)
    dfn.to_csv(args.o, sep="\t", index=False)
