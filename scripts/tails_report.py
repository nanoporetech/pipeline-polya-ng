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
    description='Overview report of poly(A) tail lengths.')
parser.add_argument(
    '-i', metavar='input', type=str, help="Input TSV.", required=True)
parser.add_argument(
    '-m', metavar='medians', type=str, help="Medians TSV.", required=True)
parser.add_argument(
    '-r', metavar='report', type=str, help="Report PDF.", required=True)
parser.add_argument('-x', action="store_true", help="Plot per-transcript distributions.", default=False)


def _make_distplot(data, title, label, xlab, ylab, pages):
    """ Make distplot with median. """
    ax = sns.distplot(data, kde=False, hist_kws={"label": label}, norm_hist=False)
    ax.set_title(title)
    ax.set_xlabel(xlab)
    ax.set_xlabel(xlab)
    ax.legend(loc='best')
    pages.savefig()
    plt.clf()


def _make_boxplot(df, med, title, pages):
    """ Make boxplot. """
    ax = sns.boxplot(x="polya_length", data=df, showfliers=False, orient='v')
    ax.set_title(title)
    ax.text(0, med + 0.5, np.round(med, 2), horizontalalignment='center', size='x-small', color='w', weight='semibold')
    pages.savefig()
    plt.clf()


if __name__ == '__main__':
    args = parser.parse_args()
    pages = PdfPages(args.r)
    sns.set_style("whitegrid")

    tails = pd.read_csv(args.i, sep="\t")

    mdf = tails[['contig', 'polya_length']].groupby(['contig']).polya_length.agg(['median', 'count']).reset_index()
    mdf = mdf.rename(columns={"median": "polya_length"})
    mdf.sort_values(by="count", ascending=False, inplace=True)
    mdf.to_csv(args.m, sep="\t", index=False)

    _make_distplot(tails.polya_length.values, title="Global tail length distribution.", label="Median: {:.2f}".format(tails.polya_length.median()), xlab="Tail length", ylab="Count", pages=pages)
    _make_boxplot(tails, med=tails.polya_length.median(), title="Global tail length distibution (no outliers).", pages=pages)

    if args.x:
        for tr in mdf.contig.values:
            med = mdf[mdf.contig == tr].polya_length.values[0]
            _make_distplot(tails[tails.contig == tr].polya_length, title="Tail length distribution: {}".format(tr), label="Median: {:.2f}".format(med), xlab="Tail length", ylab="Count", pages=pages)

            _make_boxplot(tails[tails.contig == tr], med=med, title="Tail length distibution (no outliers): {}".format(tr), pages=pages)

    pages.close()
