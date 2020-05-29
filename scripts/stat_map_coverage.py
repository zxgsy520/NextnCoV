#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import logging
import argparse
import collections


LOG = logging.getLogger(__name__)

__version__ = "1.0.0"
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


def read_stdin(sep=None):

    for line in sys.stdin:
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)


def read_depth(file, depth):

    depth = depth.strip().split(',')
    data = collections.OrderedDict()
    reflen = 0
    mapbase = 0

    for i in depth:
        data[i] = 0

    if file=='':
        fd = read_stdin('\t')
    else:
        fd = read_tsv(file, '\t')

    for line in fd:
        dp = float(line[2])
        reflen += 1
        mapbase += dp
        for i in depth:
            if int(i)>dp:
                break
            data[i] += 1

    return data, reflen, mapbase


def stat_map_coverage(files, clean, depth, out):

    cdata = {}

    for line in read_tsv(clean, '\t'):
        if line[0]=="Sample":
            continue
        cdata[line[0]] = int(line[1].replace(',', ''))

    fo = open(out, 'w')
    string = ""

    for i in depth.strip().split(','):
        string += '{0}x cov rate(%)\t'.format(i)

    fo.write('#Sample\tRefer Length(bp)\t{0}Reads Mapping bases(bp)\tReads mapping rate(%)\n'.format(string))

    for file in files:
        name = file.split('/')[-1]
        string = ""
        data, reflen, mapbase = read_depth(file, depth)

        if '--' in name:
            name = name.split('--')[1].split('.bam')[0]
        else:
            name = name.split('.')[0]

        for i in data:
            string += '{0:.2f}\t'.format(data[i]*100.0/reflen)

        fo.write('{0}\t{1:,}\t{2}{3:,}\t{4:.2f}\n'.format(name, reflen, string, mapbase, mapbase*100.0/cdata[name]))

    fo.close()


def add_help(parser):

    parser.add_argument('-i', '--input', nargs='+', metavar='FILE', type=str, required=True,
        help='Set the input file.')
    parser.add_argument('-c', '--clean', metavar='FILE', type=str, required=True,
        help='Statistics of clean reads after allegations.')
    parser.add_argument('-d', '--depth', metavar='LIST', type=str, default='3,10,50',
        help='Set the statistics coverage depth list, default=3,10,50')
    parser.add_argument('-o', '--out', metavar='STR', type=str, default='out',
        help='Set the prefix of the output file.')

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
    stat_map_coverage.py  Count the number of bases with different coverage depths.
note:
    The format of the input file:
        The first behavior contig id, the second generation behavior base position, and the third behavior coverage depth.
        Bases covering a depth of 0 also require output.

attention:
    stat_map_coverage -i *.depth -c clean.stat_reads.tsv
    stat_map_coverage -i *.depth -c clean.stat_reads.tsv -d 3,10,50 -o name

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_help(parser).parse_args()

    stat_map_coverage(args.input, args.clean, args.depth, args.out)


if __name__ == "__main__":

    main()
