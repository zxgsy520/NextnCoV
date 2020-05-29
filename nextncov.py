#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

from nextncov.parser import *
from nextncov.ncovqc import ncovqc
from nextncov.ncovsnp import ncovsnp
from nextncov.ncovann import ncovann
from nextncov.ncovall import ncovall
from nextncov import __version__, __email__, __author__


def add_nextncov_parser(parser):

    subparsers = parser.add_subparsers(
        title='command',
        dest='commands')
    subparsers.required = True

    ncovqc_parser = subparsers.add_parser('ncovqc', help="quality control")
    ncovqc_parser = add_ncovqc_args(ncovqc_parser)
    ncovqc_parser.set_defaults(func=ncovqc)

    ncovsnp_parser = subparsers.add_parser('ncovsnp', help="snp analysis")
    ncovsnp_parser = add_ncovsnp_args(ncovsnp_parser)
    ncovsnp_parser.set_defaults(func=ncovsnp)

    ncovann_parser = subparsers.add_parser('ncovann', help="annotation analysis")
    ncovann_parser = add_ncovann_args(ncovann_parser)
    ncovann_parser.set_defaults(func=ncovann)

    ncovall_parser = subparsers.add_parser('ncovall', help="nCoV analysis")
    ncovall_parser = add_ncovall_args(ncovall_parser)
    ncovall_parser.set_defaults(func=ncovall)

    return parser


def main():

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
Performing coronavirus analysis

version: %s
contact:  %s <%s>\
        """ % (__version__, " ".join(__author__), __email__))

    parser = add_nextncov_parser(parser)
    args = parser.parse_args()

    args.func(args)

    return parser.parse_args()


if __name__ == "__main__":
    main()
