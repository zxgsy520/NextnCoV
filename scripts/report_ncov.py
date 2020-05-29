#!/export/software/Bases/Python/v2.7.15/bin/python
# -*- coding: utf-8 -*-

import sys
reload(sys)
sys.setdefaultencoding('utf8')
import json
import argparse
import logging
import os.path
import shutil
import collections

from datetime import datetime
from jinja2 import Template
from docxtpl import DocxTemplate
try:
    from ConfigParser import ConfigParser
except:
    from configparser import ConfigParser

LOG = logging.getLogger(__name__)

__version__ = "1.0.0"
__author__ = ("Xingguo Zhang",)
__email__ = "1131978210@qq.com"
__all__ = []


def mkdir(d):
    """
    from FALCON_KIT
    :param d:
    :return:
    """
    d = os.path.abspath(d)
    if not os.path.isdir(d):
        LOG.debug('mkdir {!r}'.format(d))
        os.makedirs(d)
    else:
        LOG.debug('mkdir {!r}, {!r} exist'.format(d, d))

    return d


def read_config(cfg):
    """
    read config fron ini
    :param cfg:
    :return:
    """
    check_paths(cfg)

    r = {}
    config = ConfigParser()
    config.read(cfg)

    for section in config.sections():
        r[section] = {}

        for option in config.options(section):
            value = config.get(section, option).strip().decode("utf-8")
            r[section][option] = value

    return r


def read_tsv(file):

    r = []

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        r.append(line.split("\t"))

    return r


def remove_comma(string):

    string = string.replace(',', '')

    return string


def read_tsv(file, sep=None):

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)


def read_table_data(raw, clean):

    data = collections.OrderedDict()
    r = []
    rdata = 0
    cdata = 0
    sample_num = 0

    for line in read_tsv(raw, '\t'):
        data[line[0]] = [line[0], line[1], line[2]]
        rdata += int(remove_comma(line[1]))
        sample_num += 1
    #Sample	Raw Bases(bp)	Raw Reads	Clean Bases(bp)	Clean Reads	Mean Length(bp)	Effective(%)  Off Target
    for line in read_tsv(clean, '\t'):
        rlist = data[line[0]]
        cdata += int(remove_comma(line[1]))
        bontarget = "{0:.2f}".format(int(remove_comma(line[1]))*100.0/int(remove_comma(rlist[1])))
        rontarget = int(remove_comma(line[2]))*100.0/int(remove_comma(rlist[2]))
        if rontarget>=100:
            rontarget = "{0:.2f}".format(99.99)
        else:
            rontarget = "{0:.2f}".format(rontarget)
        dlist = rlist + [line[1], line[2], line[3], bontarget, rontarget]

        r.append(dlist)

    return r, rdata, cdata, sample_num


def read_table_coverage(file):

    return read_tsv(file)


def read_table_snp(file):

    return read_tsv(file)


def run_report(author, reviewer, project, pid, sequencer,
    raw, clean, coverage, snp, ncov_json,
    figure_identity, figure_coverage, tpl_html, out_dir):

    out_dir = mkdir(out_dir)
    now = datetime.now()

    data_table, rdata, cdata, sample_num = read_table_data(raw, clean)
    coverage_table = read_table_coverage(coverage)
    snp_table = read_table_snp(snp)

    r = {
        "author": author,
        "reviewer": reviewer,
        "year": now.year,
        "month": now.month,
        "day": now.day,
        "project": project,
        "id": pid,
        "sequencer": sequencer,
        "sample_num": sample_num,
        "raw_bases": "{0:,}".format(rdata),
        "clean_bases": "{0:,}".format(cdata),
        "table_data": data_table,
        "table_coverage": coverage_table,
        "table_snp": snp_table,
        "software": {},
        "database": {}

    }

    #r.update(read_config(cfg)["general"])
    for i in [ncov_json]:
        with open(i) as fh:
            js = json.load(fh)
            for k in js:
                r[k].update(js[k])

    # html_report
    #print(r)
    for i in ["images", "static"]:
        temp = os.path.join(out_dir, i)
        if os.path.exists(temp):
            shutil.rmtree(temp)
        shutil.copytree(os.path.join(tpl_html, i), temp)

    shutil.copy(figure_identity, os.path.join(out_dir, "images/name.identity.png"))
    shutil.copy(figure_coverage, os.path.join(out_dir, "images/coverage.png"))

    for j in ["index.html", "main.html"]:
        tpl = Template(open(os.path.join(tpl_html, j)).read().decode("utf-8"))

        with open(os.path.join(out_dir, j), "w") as fh:
            fh.write(tpl.render(r).encode("utf-8"))
    # html_report
    return r


def add_report_help(parser):

    parser.add_argument('-a', '--author', metavar='SRT', type=str, default='张兴国',
        help='Report writer，default=张兴国')
    parser.add_argument('--reviewer', metavar='SRT', type=str, default='朱明飞',
        help='Report reviewer，default=朱明飞')
    parser.add_argument('-p', '--project', metavar='SRT', type=str, required=True,
        help='project name')
    parser.add_argument('-i', '--pid', metavar='SRT', type=str, required=True,
        help='Project id')
    parser.add_argument("--platform", choices=["PromethION", "GridION", "RSII", "Sequel"], default="PromethION",
        help='Select sequencing platform, default=PromethION.')
    parser.add_argument('-o', '--out', metavar='FILE', type=str, default='./',
        help='Output result directory.')
    table_group = parser.add_argument_group(title='Tables', )
    table_group.add_argument('-r', '--raw', metavar='FILE', type=str, required=True,
        help='Input raw data statistics.')
    table_group.add_argument('-c', '--clean', metavar='FILE', type=str, required=True,
        help='Data statistics after inputting quality control.')
    table_group.add_argument('-cov', '--coverage', metavar='FILE', type=str, required=True,
        help='Input coverage statistics.')
    table_group.add_argument('-s', '--snp', metavar='FILE', type=str, required=True,
        help='Input snp statistics.')
    json_group = parser.add_argument_group(title='Json', )
    json_group.add_argument('-j', '--ncov_json', metavar='FILE', type=str, required=True,
        help='Input the jison file that the software runs.')
    figure_group = parser.add_argument_group(title='Figure', )
    figure_group.add_argument('-fi', '--fig_identity', metavar='FILE', type=str, required=True,
        help='Input alignment identity diagram.')
    figure_group.add_argument('-fc', '--fig_coverage', metavar='FILE', type=str, required=True,
        help='Input coverage depth map.')
    template_group = parser.add_argument_group(title='Template', )
    template_group.add_argument('--html', metavar='FILE', type=str, default='/nextomics/Pipeline/NextnCoV/v1.0.0/template/html',
        help='Input html report template(default=/nextomics/Pipeline/NextnCoV/v1.0.0/template/html).')

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

    parser = add_report_help(parser)
    args = parser.parse_args()

    run_report(args.author, args.reviewer, args.project, args.pid, args.platform,
        args.raw, args.clean, args.coverage, args.snp, args.ncov_json,
        args.fig_identity, args.fig_coverage, args.html,  args.out)


if __name__ == "__main__":
    main()
