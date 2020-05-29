#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import gzip
import logging
import argparse

from collections import OrderedDict
from os.path import abspath, expanduser

LOG = logging.getLogger(__name__)

__version__ = "1.1.0"
__author__ = ("Xingguo Zhang",)
__email__ = "1131978210@qq.com"
__all__ = []



class GffRecord(object):

    #__slots__ = ["seqid", "source", "type", "start", "end", "score", strand, phase, attrs]
    def __init__(self, seqid, source, type, start, end, score, strand, phase, attrs):
        try:
            assert "\n" not in seqid
            assert "\n" not in source
            assert "\n" not in type
            assert "\n" not in start
            assert "\n" not in end
            assert "\n" not in score
            assert "\n" not in strand
            assert "\n" not in phase
            assert "\n" not in attrs
        except AssertionError:
            raise ValueError("Invalid GFF record data")

        self.seqid = seqid
        self.source = source
        self._type = type
        self.start = int(start)
        self.end = int(end)
        assert self.start <= self.end, "%s %s %s" % (self.seqid, self.start, self.end)

        #self.length = self.end-self.start+1
        self.score = score
        self.strand = strand
        self.phase = phase
        self._attrs = attrs
        self.attributes = self._split_attr(attrs)

    @property
    def length(self):
        return self.end - self.start + 1

    def _split_attr(self, attributes):
        r = OrderedDict()
        contents = attributes.split(";")
        for content in contents:
            if not content:
                continue
            if "=" not in content:
                print("%r is not a good formated attribute: no tag!")
                continue
            tag, value = content.split("=", 1)
            r[tag] = value

        return r

    def to_string(self):
        attr = []
        for key, value in self.attributes.items():
            if key in "ID":
                attr.insert(0, "%s=%s" % (key, value))
            else:
                attr.append("%s=%s" % (key, value))
        r = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (self.seqid, self.source, self._type, self.start, self.end, self.score, self.strand, self.phase, ";".join(attr))
        return r

    @property
    def type(self):
        return self._type

    @type.setter
    def type(self, s):
        self._type = s

    @classmethod
    def from_string(cls,s):
        try:
            assert "\n" not in s
            parts = s.split("\t")
            assert len(parts) == 9
            seqid, source, type, start, end, score, strand, phase, attributes = parts
           # assert strand in "+-"
        except AssertionError:
            raise ValueError("%r not recognized as a valid GFF record" % s)

        return GffRecord(seqid, source, type, start, end, score, strand, phase, attributes)


def open_gff(fn):
    filename = abspath(expanduser(fn))
    if filename.endswith(".gz"):
        ofs = gzip.open(filename, 'r')
    elif filename.endswith(".dexta"):
        ofs = stream_stdout("undexta -vkU -w60 -i", filename)
    else:
        ofs = open(filename)

    for line in ofs.readlines():
        line = line.strip()
        if line.startswith('#'):
            continue
        if len(line) > 1:
            yield GffRecord.from_string(line)

    ofs.close()


def read_gff(file):
    
    pseqid = ""
    r = OrderedDict()

    for record in sorted(open_gff(file), key=lambda d: d.seqid):
        if record.seqid != pseqid:
            pseqid = record.seqid
            if "gene" in r:
                yield r
                r = OrderedDict()
                if record.type in "gene|mRNA":
                    r["gene"] = record
                
            print(">Feature %s" % pseqid)
        
        if record.type in "gene|mRNA":
            if "gene" not in r:
                r["gene"] = record
            elif len(r)>=2:
                yield r
                r = OrderedDict()
                r["gene"] = record
            else:
                continue
        elif record.type in "CDS|rRNA|ncRNA|tRNA":
            if record.type not in r:
                r[record.type] = [record]
                continue
            else:
                r[record.type].append(record)
                continue
        else:
             continue

    if "gene" in r:
        yield r


def obtain_locus_tag(note):
    
    r = ""
    if "Parent" in note:
        r = note["Parent"]
    elif "ID" in note:
        r = note["ID"]
    elif "Name" in note:
        r = note["Name"]
    elif "gene" in note:
        r = note["gene"]
    else:
        raise ValueError("Invalid GFF record data")

    return r


def gff2tbl(file):

    for redict in read_gff(file):
        if redict["gene"].strand == "-":
            strand = "-"
        else:
            strand = "+"
        for gtype in redict:
            if gtype == "gene":
                record = redict[gtype]

                if strand == "+":
                    start, end = record.start, record.end
                else:
                    start, end = record.end, record.start
                _attr = record.attributes
                locus_tag = obtain_locus_tag(_attr)

                if "gene" in _attr:
                    gene = "\n\t\tgene\t%s" % _attr["gene"]
                else:
                    gene = ""

                print("""\
{start}\t{end}\tgene{gene}
\t\tlocus_tag\t{locus_tag}""".format(**locals()))

            else:
                n = 0
                for record in redict[gtype]:
                    n += 1
                    if strand == "+":
                        start, end = record.start, record.end
                    else:
                        start, end = record.end, record.start
                    if n<=1:
                        print("%s\t%s\t%s" % (start, end, record.type))
                    else:
                        print("%s\t%s" % (start, end))

                    if n < len(redict[gtype]):
                        continue
                    _attr = record.attributes
                    locus_tag = obtain_locus_tag(_attr)
                    ctype = record.type

                    if "gene" in _attr:
                        gene = "\t\tgene\t%s\n" % _attr["gene"]
                    else:
                        gene = ""
                    
                    print("""\
{gene}\t\tlocus_tag\t{locus_tag}""".format(**locals()))

                    for k in ["EC_number", "inference", "note", "product"]:
                        if k in record.attributes:
                            for i in record.attributes[k].split(","):
                                print("\t\t%s\t%s" % (k, i.replace("%2C", ",")))


def set_args():
    args = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                   description="""
description:

version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    args.add_argument("file", help="")

    return args.parse_args()


def main():
    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )

    args = set_args()

    gff2tbl(args.file)


if __name__ == "__main__":
    main()

