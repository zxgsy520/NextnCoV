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
from nextncov.parser import add_ncovsnp_args

LOG = logging.getLogger(__name__)

__author__ = ("Xingguo Zhang",)
__email__ = "1131978210qq.com"
__version__ = "v1.0.0"


def stat_mapcover_snp(reads, clean, depths, snps, job_type, work_dir, out_dir):

    task = Task(
        id="stat_cover_snp",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1",
        script="""
export PATH={python}:$PATH

if [ -f "{clean}" ];then
    cp {clean} clean.stat_reads.tsv
else
    python {scripts}/stat_barcode.py --input {reads} --out clean.stat_reads.tsv >data.js
fi

python {scripts}/stat_map_coverage.py --input {depths} --clean clean.stat_reads.tsv --out stat_map_coverage.tsv
python {scripts}/stat_snp_gff.py --input {snps} >stat.snp.tsv
cp stat_map_coverage.tsv stat.snp.tsv {out_dir}
""".format(python=PYTHON_BIN,
            scripts=SCRIPTS,
            clean=clean,
            reads=reads,
            depths=depths,
            snps=snps,
            out_dir=out_dir)
        )

    return task


def create_ncovsnp_tasks(read, name, reffa, refgb, thread, job_type, work_dir, out_dir):

    option = OrderedDict()
    option["samtools"] = {
        "version": get_version(SOFTWARE_VERSION["samtools"]),
        "option": "default"
    }
    option["medaka"] = {
        "version": get_version(SOFTWARE_VERSION["medaka"]),
        "option": "medaka_variant (parameters: default)"
    }
    option["snippy"] = {
        "version": get_version(SOFTWARE_VERSION["snippy"]),
        "option": "--minfrac 0.9"
    }

    snp_task = Task(
        id="medaka_%s" % name,
        work_dir=work_dir,
        type=job_type,
        option="-pe smp %s" % thread,
        script="""
export PATH={samtools}:{minimap}:{medaka}:{python}:$PATH
source {activate} medaka

minimap2 -ax map-ont -t {thread} {reffa} {read} |samtools sort -@ 5 - >{name}.sort.bam
samtools index {name}.sort.bam
samtools depth -a {name}.sort.bam -d 0 --reference {reffa} >{name}.depth.xls
medaka_variant -f {reffa} -i {name}.sort.bam -t {thread} -b 1000 -p {name}
cp medaka_variant/round_1.vcf {name}.raw.vcf
python {script}/vcf2fasta.py --vcf {name}.raw.vcf --depth {name}.depth.xls --refer {reffa} >{name}.raw_concencus.fasta
rm -rf medaka_variant {name}.sort.bam
#cp {name}.concencus.fasta {out_dir}
""".format(minimap=MINIMAP_BIN,
            samtools=SAMTOOLS_BIN,
            medaka=MEDAKA_BIN,
            activate=MEDAKA_ENV,
            script=SCRIPTS,
            python=PYTHON_BIN,
            read=read,
            name=name,
            reffa=reffa,
            thread=thread,
            out_dir=out_dir)
        )

    snippy_task = Task(
        id="snippy_%s" % name,
        work_dir=work_dir,
        type=job_type,
        option="-pe smp %s" % thread,
        script="""
export PATH={snippy}:{python}:$PATH
python {script}/plot_depth_stat.py -i {name}.depth.xls -w 1 -o {name}
snippy --cpus {thread} --outdir {name} --minfrac 0.9 --ref {refgb} -ctgs {name}.raw_concencus.fasta
python {script}/filter_snp_gff.py -i {name}/snps.gff -d {name}.depth.xls >{name}.snps.gff
python {script}/snp_gff2fasta.py --gff {name}.snps.gff --refer {reffa} >{name}.concencus.fasta
cp {name}/snps.vcf {name}.snps.vcf
rm -rf {name}
cp {name}.concencus.fasta {out_dir}
cp {name}.depth.png {name}.depth.pdf {name}.snps.gff {out_dir}
""".format(snippy=SNIPPY_BIN,
            python=PYTHON_BIN,
            script=SCRIPTS,
            refgb=refgb,
            reffa=reffa,
            name=name,
            thread=thread,
            out_dir=out_dir)
        )

    snippy_task.set_upstream(snp_task)

    return snp_task, snippy_task, os.path.join(out_dir, "%s.concencus.fasta" % name), option


def run_ncovsnp(reads, reffa, refgb, thread, concurrent, refresh, job_type, work_dir, out_dir, clean=""):

    reads = check_paths(reads)
    reffa = check_path(reffa)
    refgb = check_path(refgb)
    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)
    options = {
        "software": OrderedDict(),
        "database": OrderedDict()
    }

    dag = DAG("ncovsnp")
    option = OrderedDict()
    depths = os.path.join(work_dir, "*/*.depth.xls")
    snps = os.path.join(work_dir, "*/*.snps.gff")

    stat_map_task = stat_mapcover_snp(
        reads=' '.join(reads),
        clean=clean,
        depths=depths,
        snps=snps,
        job_type="local",
        work_dir=work_dir,
        out_dir=out_dir)
    genomes = []

    for read in reads:
        name = read.split('/')[-1]
        if '--' in name:
            name = name.split('--')[0].split('.')[0]
        else:
            name = name.split('.')[0]

        name_work = mkdir(os.path.join(work_dir, name))
        snp_task, snippy_task, concencus, option = create_ncovsnp_tasks(
            read=read,
            name=name,
            reffa=reffa,
            refgb=refgb,
            thread=thread,
            job_type=job_type,
            work_dir=name_work,
            out_dir=out_dir)
        genomes.append(concencus)
        dag.add_task(snp_task)
        dag.add_task(snippy_task)
        stat_map_task.set_upstream(snp_task)
        stat_map_task.set_upstream(snippy_task)
    options["software"] = option
    dag.add_task(stat_map_task)
    do_dag(dag, concurrent, refresh)

    return genomes, options


def ncovsnp(args):

    genomes, options = run_ncovsnp(
        reads=args.reads,
        reffa=args.reference,
        refgb=args.genbank,
        thread=args.thread,
        concurrent=args.concurrent,
        refresh=args.refresh,
        job_type=args.job_type,
        work_dir=args.work_dir,
        out_dir=args.out_dir)

    with open(os.path.join(args.out_dir, "ncovsnp.json"), "w") as fh:
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
    ncovsnp Perform snp analysis.

attention:
    ncovsnp --reads name.reads.fq --reference reference.fa --genbank reference.gb

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    parser = add_ncovsnp_args(parser)
    args = parser.parse_args()

    ncovsnp(args)


if __name__ == "__main__":

    main()
