#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import logging
import argparse

LOG = logging.getLogger(__name__)

__version__ = "1.0.0"
__author__ = ("Xingguo Zhang",)
__email__ = "1131978210@qq.com"
__all__ = []


def sed(file, old, new, rmax):

    for line in open(file, 'r'):
        line = line.strip().replace(old, new, rmax)
        print(line)


def add_help_args(parser):

    parser.add_argument("input", help="")
    parser.add_argument("--old", help="")
    parser.add_argument("--new", help="")
    parser.add_argument("--rmax", type=int, default=10, help="")

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

    parser = add_help_args(parser)
    args = parser.parse_args()

    sed(args.input, args.old, args.new, args.rmax)


if __name__ == "__main__":
    main()
