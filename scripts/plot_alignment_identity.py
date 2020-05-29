#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import os
import sys
import gzip
import logging
import argparse
import matplotlib
matplotlib.use('Agg')

from matplotlib import pyplot as plt
#matplotlib v3.2.1

LOG = logging.getLogger(__name__)

__version__ = "1.2.0"
__author__ = ("Xingguo Zhang",)
__email__ = "113178210@qq.com"
__all__ = []


def read_tsv(file, sep=None):
    """
    read tsv joined with sep
    :param file: file name
    :param sep: separator
    :return: list
    """
    LOG.info("reading message from %r" % file)

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)


def class_identity_coverage(file, gather=100):

    identity = []
    coverage = []
    clist = []
    slist = []
    ident = []
    perce = []
    n = 0

    for line in read_tsv(file, '\t'):
        n+=1
        ident.append(float(line[1]))
        perce.append(float(line[2]))

        if n<gather:
            continue

        identity.append(sum(ident)/gather)
        imin = min(ident)
        if imin>=95:
            csize = 5
        elif imin>=90:
            csize = 4
        elif imin>=85:
            csize = 3
        elif imin>=80:
            csize = 2
        else:
            csize = 0
        clist.append(csize)

        coverage.append(sum(perce)/gather)
        cmin = min(perce)
        if cmin>=95:
            ssize = 50
        elif cmin>=90:
            ssize = 40
        elif cmin>=85:
            ssize = 30
        else:
            ssize = 20
        slist.append(ssize)
        ident = []
        perce = []
        n = 0

    return identity, coverage, clist, slist


def plot_identity(file, out, gather):

    identity, coverage, clist, slist = class_identity_coverage(file, gather)

    plt.style.use('ggplot')
    fig, ax = plt.subplots(figsize=(10, 6),)
    ax.spines['top'].set_visible(False) #去掉上边框
    ax.spines['bottom'].set_visible(False) #去掉下边框
    ax.spines['left'].set_visible(False) #去掉左边框
    ax.spines['right'].set_visible(False) #去掉右边框

    ax.scatter(coverage, identity, c=clist, s=slist, edgecolor="black", alpha=0.5)
#    ax.set_xlim(95, 100)
#    ax.set_ylim(85, 100)
    font = {'weight': 'bold','size': 12,}
    ax.set_ylabel('Alignment identity %', font)
    ax.set_xlabel('Percentage of read aliged', font)
    plt.xticks()
    plt.savefig("%s.alignment_identity.pdf" % out)
    plt.savefig("%s.alignment_identity.png" % out, dpi=700)


def add_hlep(parser):

    parser.add_argument('-i', '--input', metavar='STR', type=str, required=True,
        help='Input file.')
    parser.add_argument('-g', '--gather', metavar='INT', type=int, default=200,
        help='ISet the number of gather,default=200.')
    parser.add_argument('-o', '--out', metavar='STR',type=str, default="out",
        help='Out name')

    return parser


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
    description='''
name:
    plot_alignment_identity.py  Draw the alignment of the reads and assembly results with identity

attention:
    plot_alignment_identity.py -i identity_coverage.tsv
    plot_alignment_identity.py -i identity_coverage.tsv -o txt
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep(parser).parse_args()

    plot_identity(args.input, args.out, args.gather)


if __name__ == "__main__":

    main()
