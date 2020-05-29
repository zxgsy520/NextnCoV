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


AMINO_ACID = {"Gly": "G", "Ala": "A", "Val": "V", "Leu":"L", "Ile": "I",
    "Pro": "P", "Phe": "F", "Tyr": "Y", "Trp": "W", "Ser": "S",
    "Thr": "T", "Cys": "C", "Met": "M", "Asn": "N", "Gln": "Q",
    "Asp": "D", "Glu":"E", "Lys": "K", "Arg": "R", "His": "H"}


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
            seq += "%s\n" % line.strip(">")
            continue
        if line.startswith(">"):
            line = line.strip(">")
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



def supper(string):

    string = string[0].upper() + string[1::]

    return string


def read_gff(file):

    data = collections.OrderedDict()

    for line in read_tsv(file, '\t'):
        line = line[8].strip().split()

        if line[8]=='synonymous_variant':
            continue

        name = line[11]
        line = re.findall(r'[0-9]+|[a-z]+', line[10].lower())
        convert = [AMINO_ACID[supper(line[1])], AMINO_ACID[supper(line[3])], int(line[2])]

        if name not in data:
            data[name] = []

        data[name].append(convert)


    return data


def read_protein(file):

    data = collections.OrderedDict()

    for line in read_fasta(file):
        name = line[0].split(' [')[1]
        name = name.strip(']').split('=')[1]
        data[name] = list(line[1])
    return data


def snp2protein(gff, protein, out):

    snp = read_gff(gff)
    data = read_protein(protein)

    for i in snp:
        for line in snp[i]:
            data[i][line[2]-1] = line[1]

    fh = open(out, 'w')
    for i in data:
        fh.write(">%s\n%s\n" % (i, ''.join(data[i])))

    fh.close()


def add_help(parser):

    parser.add_argument('-sg', '--snpgff', metavar='FILE', type=str, required=True,
        help='Input snp annotation result file (snp.gff).')
    parser.add_argument('-p', '--protein', metavar='FILE', type=str, default='/nextomics/Pipeline/NextnCoV/v1.0.0/database/ref_protein.fa',
        help='Input reference protein sequence, default=/nextomics/Pipeline/NextnCoV/v1.0.0/database/ref_protein.fa.')
    parser.add_argument('-o', '--out', metavar='FILE', type=str, default='out.fa',
        help='Output file.')

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
    snp2protein.py  Modify protein sequence by SNP annotation results.

attention:
    snp2protein -sg smp.gff -o out.fa
    snp2protein -sg smp.gff -p ref_protein.fa -o out.fa

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_help(parser).parse_args()

    snp2protein(args.snpgff, args.protein, args.out)


if __name__ == "__main__":

    main()
