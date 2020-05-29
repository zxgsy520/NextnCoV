#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import logging
from Bio import SeqIO

LOG = logging.getLogger(__name__)

__version__ = "1.1.0"
__author__ = ("Xingguo Zhang, Junpeng Fan",)
__email__ = "1131978210@qq.com, jpfan@whu.edu.cn"
__all__ = []


def process_assembly(fasta, _format, min_length, completeness, topology, moltype, organism, strain, gcode):

    n = 1

    if _format:
        records = sorted(SeqIO.parse(fasta, "fasta"), key=lambda k: len(k), reverse=True)
    else:
        records = SeqIO.parse(fasta, "fasta")

    for record in records:

        if len(record) < min_length:
            break
        record.seq = record.seq.upper()
        record.description = record.description.split(None, 1)
        
        if len(record.description)>=2:
            record.description = record.description[1]
        else:
            record.description = ""

        if _format:
            record.id = _format % n
        if completeness:
            record.description += " [completeness=%s]" % completeness
        if moltype:
            record.description += " [moltype=%s]" % moltype
        if organism:
            record.description += " [organism=%s]" % organism
        if strain:
            record.description += " [strain=%s]" % strain
        if gcode:
            record.description += " [gcode=%s]" % gcode
        if topology:
            record.description += " [topology=%s]" % topology

        print(record.format("fasta"))
        n += 1


def add_args(parser):

    parser.add_argument("fasta", help="")
    parser.add_argument("--format", help="contig name format")
    parser.add_argument("--min_length", type=int, help="")
    parser.add_argument("--completeness", type=str, help="")
    parser.add_argument("--topology", choices=["circular", "linear"], default="circular",
        help="topology, default=circular")
    parser.add_argument("--moltype", choices=["DNA", "ss-RNA"], default="DNA",
        help="moltype, default=DNA.")
    parser.add_argument("--organism", help="organism name")
    parser.add_argument("--strain", help="strain")
    parser.add_argument("--gcode", type=int, help="genetic code")

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

    parser = add_args(parser)
    args = parser.parse_args()
    process_assembly(args.fasta, args.format, args.min_length, args.completeness, args.topology, args.moltype, args.organism, args.strain, args.gcode)


if __name__ == "__main__":
    main()

