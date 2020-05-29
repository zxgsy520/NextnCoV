#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import json
import logging
import argparse

from collections import OrderedDict
from nextncov.config import *
from nextncov.common import check_path, check_paths, mkdir, get_version
from dagflow import DAG, Task, ParallelTask, do_dag
from nextncov.parser import add_ncovqc_args

LOG = logging.getLogger(__name__)

__author__ = ("Xingguo Zhang",)
__email__ = "1131978210qq.com"
__version__ = "v1.0.0"


def stat_reads_task(reads, name, thread, job_type, work_dir, out_dir):

    task = Task(
        id="stat_reads_%s" % name,
        work_dir=work_dir,
        type=job_type,
        option="-pe smp %s" % thread,
        script="""
export PATH={python}:$PATH
python {scripts}/stat_barcode.py --input {reads} --out {name}.stat_reads.tsv >data.js
cp {name}.stat_reads.tsv {out_dir}
""".format(scripts=SCRIPTS,
            python=PYTHON_BIN,
            reads=reads,
            name=name,
            out_dir=out_dir)
    )

    return task, os.path.join(work_dir, "%s.stat_reads.tsv" % name)



def map_ref_tasks(reads, names, reference, thread, job_type, work_dir, out_dir, platform="PromethION"):

    option = OrderedDict()
    option["minimap2"] = {
        "version": get_version(SOFTWARE_VERSION["minimap2"]),
        "option": "%s" % SEQUENCER[platform]["minimap2"]
    }
    option["map_paf2reads"] = {
        "version": get_version(SOFTWARE_VERSION["map_paf2reads"]),
        "option": "default"
    }

    tasks = ParallelTask(
        id="minimap",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp %s" % thread,
        script="""
export PATH={minimap}:{map_paf2reads}:$PATH
minimap2 -t {thread} {x} {reference} {{reads}} |map_paf2reads -r {{reads}} -s {{name}}.identity_coverage.tsv >{{name}}.clean_reads.fa
{scripts}/plot_alignment_identity.py -i {{name}}.identity_coverage.tsv -o {{name}}
cp {{name}}.alignment_identity.png {{name}}.alignment_identity.pdf {out_dir}
gzip -c {{name}}.clean_reads.fa >{out_dir}/{{name}}.clean_reads.fa.gz
""".format(
            minimap=MINIMAP_BIN,
            map_paf2reads=MAP_PAF2READS,
            scripts=SCRIPTS,
            reference=reference,
            x=SEQUENCER[platform]["minimap2"],
            thread=thread,
            out_dir=out_dir
        ),
        reads=reads,
        name=names
    )

    return tasks, os.path.join(work_dir, "*.clean_reads.fa"), option


def run_ncovqc(reads, reference, thread, job_type, concurrent, refresh, work_dir, out_dir):

    reference = check_path(reference)
    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)
    reads = check_paths(reads)
    names = []

    for i in reads:
        name = i.split('/')[-1]
        if '--' in name:
            name = name.split('--')[1].split('.bam')[0]
        else:
            name = name.split('.')[0]
        names.append(name)

    options = {
        "software": OrderedDict(),
        "database": OrderedDict()
    }

    dag = DAG("ncovqc")
    raw_task, raw_stat = stat_reads_task(
        reads=" ".join(reads),
        name="raw",
        thread=thread,
        job_type=job_type,
        work_dir=work_dir,
        out_dir=out_dir)

    map_tasks, clean_reads, option = map_ref_tasks(
        reads=reads,
        names=names,
        reference=reference,
        thread=thread,
        job_type=job_type,
        work_dir=work_dir,
        out_dir=out_dir)
    options["software"] = option

    clean_task, clean_stat = stat_reads_task(
        reads=clean_reads,
        name="clean",
        thread=thread,
        job_type=job_type,
        work_dir=work_dir,
        out_dir=out_dir)
    dag.add_task(raw_task)
    dag.add_task(*map_tasks)
    dag.add_task(clean_task)
    clean_task.set_upstream(*map_tasks)
    do_dag(dag, concurrent, refresh)

    return clean_reads, clean_stat, options


def ncovqc(args):

    clean_reads, clean_stat, options = run_ncovqc(
        reads=args.reads,
        reference=args.reference,
        thread=args.thread,
        job_type=args.job_type,
        concurrent=args.concurrent,
        refresh=args.refresh,
        work_dir=args.work_dir,
        out_dir=args.out_dir)

    with open(os.path.join(args.out_dir, "ncovqc.json"), "w") as fh:
        json.dump(options, fh, indent=2)


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
        description='''
name:
    ncovqc Quality control of data.

attention:
    ncovqc --reads name.reads.fq --reference reference.fa

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    parser = add_ncovqc_args(parser)
    args = parser.parse_args()

    ncovqc(args)


if __name__ == "__main__":

    main()
