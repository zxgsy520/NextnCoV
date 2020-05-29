#!/usr/bin/env python
#coding:utf-8

import os
import re
import sys
import json
import logging
import argparse

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


def read_snp_gff(file):

    print("##gff-version 3")
    for line in read_tsv(file, '\t'):
        note = line[8].strip().split()

        if "note=snp" not in note:
            continue

        if ("synonymous_variant" in note) or ("missense_variant" in note) or ("intergenic_region" in note):
            print('\t'.join(line))

    return 0


def add_help(parser):

    parser.add_argument('-i', '--input', metavar='FILE', type=str, required=True,
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
    filter_snp_gff.py  filter snp.

attention:
    filter_snp_gff.py -i snp.gff >filter_snp.gff

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_help(parser).parse_args()

    read_snp_gff(args.input)

if __name__ == "__main__":

    main()
