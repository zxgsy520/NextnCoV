#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import json
import logging
import argparse

from nextncov.config import *
from nextncov.common import check_path, check_paths, mkdir, get_version, read_files, read_tsv
from dagflow import DAG, Task, ParallelTask, do_dag
from nextncov.ncovqc import run_ncovqc
from nextncov.ncovsnp import run_ncovsnp
from nextncov.ncovann import run_ncovann
from nextncov.parser import add_ncovall_args

LOG = logging.getLogger(__name__)

__author__ = ("Xingguo Zhang",)
__email__ = "1131978210qq.com"
__version__ = "v1.0.0"


def create_report_task(author, reviewer, project, pid, raw, clean, coverage, snp, ncov_json, fig_identity, fig_coverage, job_type, work_dir, out_dir):

    task = Task(
        id="ncov_report",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1",
        script="""
#export PATH={python}:$PATH

{scripts}/report_ncov.py --author {author} --reviewer {reviewer} \\
--project {project} --pid {pid} \\
--raw {raw} \\
--clean {clean} \\
--coverage {coverage} \\
--snp {snp} \\
--ncov_json {ncov_json} \\
--fig_identity {fig_identity} \\
--fig_coverage {fig_coverage} \\
--html {html} \\
--out {out_dir}
""".format(python=PYTHON_BIN,
            scripts=SCRIPTS,
            html=os.path.join(TEMPLATE, "html"),
            author=author,
            reviewer=reviewer,
            project=project,
            pid=pid,
            raw=raw,
            clean=clean,
            coverage=coverage,
            snp=snp,
            ncov_json=ncov_json,
            fig_identity=fig_identity,
            fig_coverage=fig_coverage,
            out_dir=out_dir)
        )

    return task


def run_ncovall(author, reviewer, project, pid, reads, primer, reffa, refgb, refgff, thread, concurrent, refresh, job_type, work_dir, out_dir):

    primer = check_path(primer)
    reffa = check_path(reffa)
    refgb = check_path(refgb)
    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)
    reads = check_paths(reads)

    work_dict = {
        "data": "01_data",
        "snp": "02_snp",
        "ann": "03_annotate",
    }

    for k, v in work_dict.items():
        mkdir(os.path.join(work_dir, v))
        mkdir(os.path.join(out_dir, v))

    clean_reads, clean_stat, options = run_ncovqc(
        reads=reads,
        reference=primer,
        thread=thread,
        job_type=job_type,
        concurrent=concurrent,
        refresh=refresh,
        work_dir=os.path.join(work_dir, work_dict["data"]),
        out_dir=os.path.join(out_dir, work_dict["data"]))

    reads = clean_reads.split('/')
    clean_reads = read_files('/'.join(reads[0:-1]), reads[-1])

    genomes, noption = run_ncovsnp(
        reads=clean_reads,
        reffa=reffa,
        refgb=refgb,
        thread=thread,
        concurrent=concurrent,
        refresh=refresh,
        job_type=job_type,
        work_dir=os.path.join(work_dir, work_dict["snp"]),
        out_dir=os.path.join(out_dir, work_dict["snp"]),
        clean=clean_stat)

    run_ncovann(
        genomes=genomes,
        refgff=refgff,
        concurrent=concurrent,
        refresh=refresh,
        job_type="local",
        work_dir=os.path.join(work_dir, work_dict["ann"]),
        out_dir=os.path.join(out_dir, work_dict["ann"])
        )

    options['software'].update(noption['software'])

    with open(os.path.join(out_dir, "nextncov.json"), "w") as fh:
        json.dump(options, fh, indent=2)

    out_data = os.path.join(out_dir, work_dict["data"])
    out_snp = os.path.join(out_dir, work_dict["snp"])
    best_name = ""
    best_cover = 0

    for line in read_tsv(os.path.join(out_snp, "stat_map_coverage.tsv"), '\t'):
        if float(line[4])>=best_cover:
            best_cover = float(line[4])
            best_name = line[0]

    dag = DAG("report")

    report_task = create_report_task(
        author=author,
        reviewer=reviewer,
        project=project,
        pid=pid,
        raw=os.path.join(out_data, "raw.stat_reads.tsv"),
        clean=os.path.join(out_data, "clean.stat_reads.tsv"),
        coverage=os.path.join(out_snp, "stat_map_coverage.tsv"),
        snp=os.path.join(out_snp, "stat.snp.tsv"),
        ncov_json=os.path.join(out_dir, "nextncov.json"),
        fig_identity=os.path.join(out_data, "%s.alignment_identity.png" % best_name),
        fig_coverage=os.path.join(out_snp, "%s.depth.png" % best_name),
        job_type="local",
        work_dir=work_dir,
        out_dir=out_dir)
    dag.add_task(report_task)
    do_dag(dag, concurrent, refresh)


def ncovall(args):

    run_ncovall(
        author=args.author,
        reviewer=args.reviewer,
        project=args.project,
        pid=args.pid,
        reads=args.reads,
        primer=args.primer,
        reffa=args.reference,
        refgb=args.genbank,
        refgff=args.gff,
        thread=args.thread,
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
name:
    ncovsnp Perform snp analysis.

attention:
    ncovsnp --reads name.reads.fq --reference reference.fa --genbank reference.gb

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    parser = add_ncovsll_args(parser)
    args = parser.parse_args()

    ncovall(args)


if __name__ == "__main__":

    main()
