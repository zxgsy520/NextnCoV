#!/usr/bin/env python
# -*- coding: utf-8 -*-

from nextncov.config import *

__all__ = ["add_ncovqc_args", "add_ncovsnp_args",  "add_ncovann_args", "add_ncovall_args"]

def add_workflow_args(parser):
    """
    add workflow arguments to parser
    :param parser: argparse object
    :return: parser
    """

    workflow_group = parser.add_argument_group(title="Workflow arguments", )
    workflow_group.add_argument("--concurrent", metavar="INT", type=int, default=10,
        help="Maximum number of jobs concurrent  (default: 10).")
    workflow_group.add_argument("--refresh", metavar="INT", type=int, default=30,
        help="Refresh time of log in seconds  (default: 30).")
    workflow_group.add_argument("--job_type", choices=["sge", "local"], default="local",
        help="Jobs run on [sge, local]  (default: local).")
    workflow_group.add_argument("--work_dir", metavar="DIR", type=str, default=".",
        help="Work directory (default: current directory).")
    workflow_group.add_argument("--out_dir", metavar="DIR", type=str, default=".",
        help="Output directory (default: current directory).")

    return parser


def add_ncovqc_args(parser):

    parser.add_argument("-r", "--reads", metavar="FILE", nargs='+', type=str, required=True,
        help="Input sequencing reads.")
    parser.add_argument("-ref", "--reference", metavar="FILE", type=str, required=True,
        help="Input reference sequence.")
    parser.add_argument("--thread",  metavar="INT", type=int, default=5,
        help="Input the number of threads running in a single program, default=5")
    parser = add_workflow_args(parser)

    return parser


def add_ncovsnp_args(parser):

    parser.add_argument("-r", "--reads", metavar="FILE", nargs='+', type=str, required=True,
        help="Input sequencing reads.")
    parser.add_argument("-ref", "--reference", metavar="FILE", type=str, default="%s/reference.fasta" % DATABASE,
        help="Input reference sequence, default=reference.fasta.")
    parser.add_argument("-gb", "--genbank", metavar="FILE", type=str, default="%s/reference.gb" % DATABASE,
        help="Input the genbank format of the reference sequence, default=reference.gb.")
    parser.add_argument("--thread",  metavar="INT", type=int, default=5,
        help="Input the number of threads running in a single program, default=5")
    parser = add_workflow_args(parser)

    return parser


def add_ncovann_args(parser):
    
    parser.add_argument("-g", "--genomes", metavar="FILE", nargs='+', type=str, required=True,
        help="Input genomes.")
    parser.add_argument("-gff", "--refgff", metavar="FILE", type=str, default="%s/reference.gff3" % DATABASE,
        help="Input reference gff, default=reference.gff3.")
    parser = add_workflow_args(parser)

    return parser


def add_ncovall_args(parser):

    parser.add_argument("-gff", "--gff", metavar="FILE", type=str, default="%s/reference.gff3" % DATABASE,
        help="Input reference gff, default=reference.gff3.")
    parser.add_argument('-a', '--author', metavar='SRT', type=str, default='张兴国',
        help='Report writer，default=张兴国')
    parser.add_argument('--reviewer', metavar='SRT', type=str, default='朱明飞',
        help='Report reviewer，default=朱明飞')
    parser.add_argument('--project', metavar='SRT', type=str, default="NextnCoV",
        help='project name,default=NextnCoV.')
    parser.add_argument('-i', '--pid', metavar='SRT', type=str, default="NextnCoV",
        help='Project id,default=NextnCoV.')
    parser.add_argument("-p", "--primer", metavar="FILE", type=str, default="%s/pm98.fasta" % DATABASE,
        help="Input sequencing primer sequence, default=pm98.fasta.")

    
    parser = add_ncovsnp_args(parser)

    return parser
