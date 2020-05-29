#!/usr/bin/env python
#coding:utf-8

import os
import re
import sys
import json
import logging
import argparse

import collections

LOG = logging.getLogger(__name__)

__version__ = "1.0.0"
__author__ = ("Xingguo Zhang",)
__email__ = "113178210@qq.com"
__all__ = []


def read_fasta(file):
    '''Read fasta file'''

    if file.endswith(".gz"):
        fp = gzip.open(file)
    elif file.endswith(".fasta") or file.endswith(".fa"):
        fp = open(file)
    else:
        raise Exception("%r file format error" % file)

    seq = ''

    for line in fp:
        line = line.strip()

        if not line:
            continue
        if not seq:
            seq += "%s\n" % line.strip(">").split()[0]
            continue
        if line.startswith(">"):
            line = line.strip(">").split()[0]
            seq = seq.split('\n')

            yield [seq[0], seq[1]]
            seq = ''
            seq += "%s\n" % line
        else:
            seq += line

    seq = seq.split('\n')
    if len(seq)==2:
        yield [seq[0], seq[1]]
    fp.close()


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


def snp_gff2fasta(snpgff, refer):

    data = {}
    name = snpgff.split('/')[-1]

    if '--' in name:
        name = name.split('--')[1].split('.')[0]
    else:
        name = name.split('.')[0]

    for line in read_fasta(refer):
        data[line[0]] = list(line[1].upper())

    for line in read_tsv(snpgff, '\t'):
        if line[0] not in data:
            continue
        position = int(line[3])-1
        note = note2dict(line[8])
        old, new = note['REF'], note['ALT']
        data[line[0]][position] = new

    for i in data:
        print('>%s\n%s' % (name, ''.join(data[i])))


def add_help(parser):

    parser.add_argument('-g', '--gff', metavar='FILE', type=str, required=True,
        help='Input snp file(snp.gff).')
    parser.add_argument('-r', '--refer', metavar='FILE', type=str, required=True,
        help='Input reference sequence file(.fa, .fasta)')

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
    snp_gff2fasta.py  Generate concencus sequence from gff file

attention:
    snp_gff2fasta.py --gff snp.gff -r reference.fasta >concencus.fasta

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_help(parser).parse_args()

    snp_gff2fasta(args.gff, args.refer)


if __name__ == "__main__":

    main()
