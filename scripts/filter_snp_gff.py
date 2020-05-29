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

__version__ = "1.1.0"
__author__ = ("Xingguo Zhang",)
__email__ = "113178210@qq.com"
__all__ = []


AMINO_ACID = {"Gly": "G", "Ala": "A", "Val": "V", "Leu":"L", "Ile": "I",
    "Pro": "P", "Phe": "F", "Tyr": "Y", "Trp": "W", "Ser": "S",
    "Thr": "T", "Cys": "C", "Met": "M", "Asn": "N", "Gln": "Q",
    "Asp": "D", "Glu":"E", "Lys": "K", "Arg": "R", "His": "H"}


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


def note2norm(seqid, note, depth):

    note = note.strip().split()
    ref, alt = note[1].split('=>')
    vep = note[8]

    if vep in ["synonymous_variant", "missense_variant"]:
        pref, psite, palt = re.findall(r'[0-9]+|[a-z]+', note[10].lower())[1::]
        pref = AMINO_ACID[supper(pref)]
        palt = AMINO_ACID[supper(palt)]
        site = re.findall(r'[0-9]+|[A-Z]+', note[9])[0]
        gene = note[11]
        norm = "ID={id};REF={ref};ALT={alt};VEP={vep},{gene}:p.{psite}{pref}>{palt},\
{gene}:c.{site}{ref}>{alt};coverage={depth};".format(
            id=seqid, ref=ref, alt=alt, vep=vep, psite=psite,
            pref=pref, palt=palt, gene=gene, site=site, depth=depth)
    else:
        if vep=="intergenic_region":
            vep = "intergenic_variant"
        norm = "ID={id};REF={ref};ALT={alt};VEP={vep}".format(
            id=seqid, ref=ref, alt=alt, vep=vep)

    return norm


def filter_snp_gff(file, depth):

    data = collections.OrderedDict()

    for line in read_tsv(file, '\t'):
        note = line[8].strip().split()
        snps = note[0].split('=')[1].upper()

        if "note=snp" not in note:
            continue
        if ("synonymous_variant" in note) or ("missense_variant" in note) or ("intergenic_region" in note):
            line[1] = "snpEff:4.3"
            line[2] = snps
            data[line[3]] = line

    print("##gff-version 3")
    for line in read_tsv(depth, '\t'):
        if line[1] not in data:
            continue
        snps = data[line[1]]
        snps[-1] = note2norm(snps[0], snps[-1], line[2])
        print('\t'.join(snps))

    return 0


def add_help(parser):

    parser.add_argument('-i', '--input', metavar='FILE', type=str, required=True,
        help='Set the input file(snp.gff).')
    parser.add_argument('-d', '--depth', metavar='FILE', type=str, required=True,
        help='Input the coverage depth file for each site(depth.xls)')

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
    filter_snp_gff.py -i snp.gff -d depth.xls >filter_snp.gff

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_help(parser).parse_args()

    filter_snp_gff(args.input, args.depth)

if __name__ == "__main__":

    main()
