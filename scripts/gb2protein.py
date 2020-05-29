#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import gzip
import logging
import argparse

from Bio import SeqIO

LOG = logging.getLogger(__name__)

__version__ = "v1.2.0"
__author__ = ("Xingguo Zhang",)
__email__ = "113178210@qq.com"
__all__ = []


def gb2protein(file):

    LOG.info("reading message from %r" % file)

    if file.endswith(".gz"):
        fh = gzip.open(file)
    else:
        fh = open(file)

    if file.endswith('.gbff.gz') or file.endswith('.gb.gz') or file.endswith('.gbff') or  file.endswith('.gb'):
        fh = SeqIO.parse(fh, "genbank")
    else:
        raise Exception("Unknown format!")

    for record in fh:
        for rna_dct in record.features :
            rna_seq = rna_dct.extract(record.seq)

            if rna_dct.type != "CDS" or len(rna_seq)%3!=0:
                continue

            try:
                protein = str(rna_dct.qualifiers['translation'][0])
                if 'locus_tag' not in rna_dct.qualifiers:
                    rna_name = str(rna_dct.qualifiers['gene'][0])
                else:
                    rna_name = str(rna_dct.qualifiers['locus_tag'][0])
            except:
                LOG.info("%r is wrong" %  rna_name)
                continue

            print(">%s\n%s" % (rna_name, protein))


def add_help_args(parser):

    parser.add_argument('genbank', help='')

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


version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    args = add_help_args(parser).parse_args()

    gb2protein(args.genbank)


if __name__ == "__main__":

    main()
