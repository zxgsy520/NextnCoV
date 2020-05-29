#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import logging
import argparse

import collections

LOG = logging.getLogger(__name__)

__version__ = "v1.2.0"
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

            yield seq[0], seq[1]
            seq = ''
            seq += "%s\n" % line
        else:
            seq += line

    seq = seq.split('\n')
    if len(seq)==2:
        yield seq[0], seq[1]


def read_tsv(file, sep=None):

    LOG.info("reading message from %r" % file)

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)


def stat_mutation(alt, mutlist):

    alt = alt.split(',')
    mut = []

    mut.append(0)
    for i in range(len(alt)):
        mut.append(0)

    for i in mutlist:
        if i in ['.']:
            continue
        mut[0] += 1
        if i in ['0']:
            continue
        i = int(i)
        mut[i] += 1

    return alt, mut


def read_vcf(file):

    for line in read_tsv(file, '\t'):
        seqid = line[0]
        pos = line[1]
        ref = line[3]
        alt, mut = stat_mutation(line[4], line[9::])

        yield [seqid, pos, ref, alt, mut]


def filter_snp(alt, mut):

    total = mut[0]
    r = [0, 0, 0, 0]

    for i in range(len(alt)):
        if len(alt[i])!=1:
            continue
        if alt[i] not in ['A', 'T', 'C', 'G']:
            continue
        if alt[i]=='A':
            r[0] = mut[i+1]
        elif alt[i]=='T':
            r[1] = mut[i+1]
        elif alt[i]=='C':
            r[2] = mut[i+1]
        else:
            r[3] = mut[i+1]
    return total, r


def stat_group_vcf(fasta, vcf):

    dseq = {}
    for seqid, seq in read_fasta(fasta):
        dseq[seqid] = seq.upper()

    if len(dseq)!=1:
        raise Exception("The program only supports mutation statistics for single reads")

    seq = list(dseq.values())[0]

    print("ID\tPOD\tREF\tNUM\tA\tT\tC\tG\tMUT_RATE(‰)")
    for line in read_vcf(vcf):
        seqid = line[0]
        pos = line[1]
        ref = line[2]

        if seq[int(pos)-1]!=ref:
            LOG.info("%s, real=%s, ref=%s" % (pos, seq[int(pos)-1], ref))
            ref = seq[int(pos)-1]
        total, mut = filter_snp(line[3], line[4])
        print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8:.2f}".format(seqid, pos, ref, total, mut[0], mut[1], mut[2], mut[3], sum(mut)*1000.0/total))

    return 0


def add_help_args(parser):

    parser.add_argument('-fa', '--fasta', metavar='FILE', type=str, required=True,
        help='Input reference gene sequence.')
    parser.add_argument('-v', '--vcf', metavar='FILE', type=str, required=True,
        help='Input aligned vcf file.')
    return parser


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
name:
     stat_group_vcf.py Statistic primer mutation rate.

attention:
     stat_group_vcf.py  --fasta ref.fa --vcf ref.vcf >stat_mutation_rate.tsv
version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    args = add_help_args(parser).parse_args()

    stat_group_vcf(args.fasta, args.vcf)


if __name__ == "__main__":

    main()
