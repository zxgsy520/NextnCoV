#!/usr/bin/env python
#coding:utf-8

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


def note2dict(note):
    
    r = {}

    for line in note.split(';'):
        if not line:
            continue
        keys, value = line.split('=')
        r[keys] = value

    return r


def read_snp_gff(file):

    dlist = [0, 0, 0, 0, 0]

    for line in read_tsv(file, '\t'):
        note = note2dict(line[8])

        if line[2]!="SNP":
            continue
        dlist[0] +=1

        if  note["VEP"]=="synonymous_variant":
            dlist[1] += 1
        elif note["VEP"]=="missense_variant":
            dlist[2] += 1
        elif  note["VEP"]=="intergenic_variant":
            dlist[3] += 1
        else:
            dlist[4] += 1

    return dlist


def list2str(rlist):

    dlist = []

    for i in rlist:
        dlist.append(str(i))

    return dlist
    

def stat_snp_gff(files):

    print('#Sample\tSNP Number\tSynonymous Variant\tMissense Variant\tIntergenic Region\tOther')
    for file in files:
        name = file.split('/')[-1]

        if '--' in name:
            name = name.split('--')[0].split('.')[0]
        else:
            name = name.split('.')[0]

        dlist = read_snp_gff(file)

        print("%s\t%s" % (name, '\t'.join(list2str(dlist))))


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
