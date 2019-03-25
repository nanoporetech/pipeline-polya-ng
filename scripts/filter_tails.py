#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import cm
import seaborn as sns

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description='Filter poly(A) tail lengths estimated by nanopolish.')
parser.add_argument(
    '-i', metavar='input', type=str, help="Input TSV.")
parser.add_argument(
    '-o', metavar='output', type=str, help="Output TSV.")
parser.add_argument(
    '-r', metavar='output', type=str, help="Report PDF.")
parser.add_argument(
    '-t', metavar='reptsv', type=str, help="Report TSV.")


if __name__ == '__main__':
    args = parser.parse_args()
    pages = PdfPages(args.r)

    tails = pd.read_csv(args.i, sep="\t")
    stats = tails.qc_tag.value_counts().to_frame()
    stats.to_csv(args.t, sep="\t", index=True)

    stats.plot(kind='barh', title="Nanopolish QC stats", fontsize=5)
    pages.savefig()
    plt.clf()
    pages.close()
    
    tails = tails[tails.qc_tag == "PASS"]
    tails.to_csv(args.o, sep="\t", index=False)

