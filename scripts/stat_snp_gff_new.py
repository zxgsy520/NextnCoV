#!/usr/bin/env python
#coding:utf-8

import os
import re
import sys
import logging
import argparse

import collections

LOG = logging.getLogger(__name__)

__version__ = "1.1.0"
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


def read_snp_gff(file):

    data = collections.OrderedDict()
    data['snp'] = [0, 0, 0, 0, 0]
    data['del'] = [0, 0, 0, 0, 0]
    data['ins'] = [0, 0, 0, 0, 0]
    data['complex'] = [0, 0, 0, 0, 0]

    for line in read_tsv(file, '\t'):
        note = line[8].strip().split()
        types = note[0].split('=')[-1]
        data[types][0] +=1

        if "synonymous_variant" in note:
            data[types][1] += 1
        elif "missense_variant" in note:
            data[types][2] += 1
        elif "intergenic_region" in note:
            data[types][3] += 1
        else:
            data[types][4] += 1

    return data


def list2str(rlist):

    dlist = []

    for i in rlist:
        dlist.append(str(i))

    return dlist


def stat_snp_gff(files):

    print('#Sample\tTypes\tNumber\tSynonymous Variant\tMissense Variant\tIntergenic Region\tOther')
    for file in files:
        name = file.split('/')[-1]

        if '--' in name:
            name = name.split('--')[0].split('.')[0]
        else:
            name = name.split('.')[0]

        data = read_snp_gff(file)
        for types in data:
            print("%s\t%s\t%s" % (name, types, '\t'.join(list2str(data[types]))))


def add_help(parser):

    parser.add_argument('-i', '--input', nargs='+', metavar='FILE', type=str, required=True,
        help='Set the input file(snp.gff).')

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
    stat_snp_gff.py  Statistics of various SNP mutation types.

attention:
    stat_snp_gff.py -i snp.gff >stat.snp.tvs
    stat_snp_gff.py -i *.snp.gff >stat.snp.tvs

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_help(parser).parse_args()

    stat_snp_gff(args.input)


if __name__ == "__main__":

    main()
