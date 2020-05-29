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


def ref_vcf(file):

    data = collections.OrderedDict()

    for line in read_tsv(file, '\t'):
        data[line[1]] = line

    return data


def filter_vcf(vcf, depth, mindepth=5, qscore=20):

    data = ref_vcf(vcf)
    name = vcf.split('/')[-1]
    if '--' in name:
        name = name.split('--')[0].split('.')[0]
    else:
        name = name.split('.')[0]

    fo = open("%s.clean.vcf" % name, "w")
    fo.write("##fileformat=VCFv4.1\n")
    fo.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
    r = []
    for line in read_tsv(depth, '\t'):
        if line[1] not in data:
            continue
        if int(line[2])<mindepth:
            continue
        lvcf = data[line[1]]
        if len(lvcf[3])!=len(lvcf[4]) or len(lvcf[4])>1:
            continue
        if float(lvcf[-1].split(':')[-1])<qscore:
            continue
        fo.write("%s\n" % ("\t".join(lvcf)))
        r.append([lvcf[0], lvcf[1], lvcf[3], lvcf[4]])
    fo.close()

    return r, name


def vcf2fasta(r, refer, name):

    data = {}

    for line in read_fasta(refer):
        data[line[0]] = list(line[1].upper())

    for line in r:
        if line[0] not in data:
            continue
        position = int(line[1])-1
        data[line[0]][position] = line[3]

    for i in data:
        print('>%s\n%s' % (name, ''.join(data[i])))


def add_help(parser):

    parser.add_argument('-v', '--vcf', metavar='FILE', type=str, required=True,
        help='Input snp file(snp.vcf).')
    parser.add_argument('-d', '--depth', metavar='FILE', type=str, required=True,
        help='Input the coverage depth file for each site(depth.xls)')
    parser.add_argument('-r', '--refer', metavar='FILE', type=str, required=True,
        help='Input reference sequence file(.fa, .fasta)')
    parser.add_argument('-md', '--mindepth', metavar='INT', type=int, default=5,
        help='Set the minimum depth of coverage for filtering, default=5.')
    parser.add_argument('-qs', '--qscore', metavar='FLOAT', type=float, default=20.0,
        help='Set the minimum quality score for filtering, default=20.0.')

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
    vcf2fasta.py  Generate concencus sequence from vcf file

attention:
    vcf2fasta.py -v snp.vcf -d depth.xls -r reference.fasta >concencus.fasta

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_help(parser).parse_args()

    r, name = filter_vcf(args.vcf, args.depth, args.mindepth, args.qscore)
    vcf2fasta(r, args.refer, name)


if __name__ == "__main__":

    main()
