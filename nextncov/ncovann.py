#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import logging
import argparse

from nextncov.config import *
from nextncov.common import check_path, check_paths, mkdir, get_version
from dagflow import DAG, Task, ParallelTask, do_dag
from nextncov.parser import add_ncovann_args

LOG = logging.getLogger(__name__)

__author__ = ("Xingguo Zhang",)
__email__ = "1131978210qq.com"
__version__ = "v1.0.0"


def create_ncovann_task(genome, name, refgff, job_type, work_dir, out_dir):

    ann_task = Task(
        id="annotate_%s" % name,
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1",
        script="""
export PATH={tbl2asn}:$PATH
{script}/process_assembly.py {genome} --topology linear --moltype ss-RNA --completeness complete --gcode 1  --organism 'Unknow' --strain {name} > {name}.genomic.fasta
{script}/sed.py {refgff} --old MN908947.3 --new {name} > {name}.genomic.gff
{script}/gff2tbl.py {name}.genomic.gff >{name}.genomic.tbl
tbl2asn -i {name}.genomic.fasta -V b -s T
mv {name}.genomic.gbf {name}.genomic.gb
{script}/gb2protein.py {name}.genomic.gb >{name}.protein.fasta
cp {name}.genomic.fasta {name}.genomic.gff {name}.genomic.sqn {name}.genomic.gb {name}.protein.fasta {out_dir}
""".format(tbl2asn=TBL2ASN_BIN,
            script=SCRIPTS,
            genome=genome,
            name=name,
            refgff=refgff,
            out_dir=out_dir)
        )

    return ann_task


def run_ncovann(genomes, refgff, concurrent, refresh, job_type, work_dir, out_dir):

    genomes = check_paths(genomes)
    refgff = check_path(refgff)
    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)
    dag = DAG("ncovann")

    for genome in genomes:
        name = genome.split('/')[-1]
        if '--' in name:
            name = name.split('--')[0].split('.')[0]
        else:
            name = name.split('.')[0]

        ann_task = create_ncovann_task(
            genome=genome,
            name=name,
            refgff=refgff,
            job_type=job_type,
            work_dir=work_dir,
            out_dir=out_dir)

        dag.add_task(ann_task)

    do_dag(dag, concurrent, refresh)

    return 0


def ncovann(args):

    run_ncovann(
        genomes=args.genomes,
        refgff=args.refgff,
        concurrent=args.concurrent,
        refresh=args.refresh,
        job_type=args.job_type,
        work_dir=args.work_dir,
        out_dir=args.out_dir)


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
        description='''

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    parser = add_ncovann_args(parser)
    args = parser.parse_args()

    ncovann(args)


if __name__ == "__main__":

    main()
