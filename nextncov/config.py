import os.path
from collections import OrderedDict

ROOT = "/nextomics/Pipeline/NextnCoV/v1.0.0"
SCRIPTS = os.path.join(ROOT, "scripts")
TEMPLATE = os.path.join(ROOT, "template")
BIN = os.path.join(ROOT, "nextncov")

#PYNTHON_bin
PYTHON_BIN="/nextomics/Software/Base/miniconda3/bin"
#MINIMAP
MINIMAP_BIN = "/software/minimap2"
#SAMTOOLS_BIN
SAMTOOLS_BIN="/export/software/Bases/samtools/v1.4/bin"
#MAP_PAF2READS
MAP_PAF2READS = "/nextomics/Software/self/map_paf2reads"
#SELF_BIN
SELF_BIN = "/nextomics/Software/self/bin"
#MEDAKA_BIN
MEDAKA_BIN = "/nextomics/Software/Base/miniconda3/envs/medaka/bin"
MEDAKA_ENV = "/nextomics/Software/Base/miniconda3/bin/activate"
#MEDAKA_BIN
SNIPPY_BIN = "/nextomics/Software/meta/snippy/v4.6.0/bin"
#TBL2ASN_BIN
TBL2ASN_BIN = "/nextomics/Software/meta/tbl2asn"


SEQUENCER = OrderedDict([
    ("RSII", {"minimap2": "-secondary=no -c -x map-pb"}
    ),
    ("Sequel", {"minimap2": "-secondary=no -c -x map-pb"}
    ),
    ("GridION", {"minimap2": "-secondary=no -c -x map-ont"}
    ),
    ("PromethION", {"minimap2": "-secondary=no -c -x map-ont"}
    ),
    ("illumina", {"minimap2": "-secondary=no -c -x sr"}
    ),
    ("mgi", {"minimap2": "-secondary=no -c -x sr"}
    ),
])

# SOFTWARE VERSION

SOFTWARE_VERSION = {
    "minimap2": {
        "GETVER": "%s/minimap2 --version" % MINIMAP_BIN,
        "REGEXP": "\d+\.\d+",
        "MINVER": "2.17",
    },
    "map_paf2reads": {
        "GETVER": "%s/map_paf2reads -h 2>&1|grep 'version'" % MAP_PAF2READS,
        "REGEXP": "\d+\.\d+\.\d+",
        "MINVER": "1.0.0",
    },
    "samtools": {
        "GETVER": "%s/samtools --version 2>&1|grep 'samtools'" % SAMTOOLS_BIN,
        "REGEXP": "\d+\.\d+",
        "MINVER": "1.4",
    },
    "medaka": {
        "GETVER": "%s/medaka --version" % MEDAKA_BIN,
        "REGEXP": "\d+\.\d+\.\d+",
        "MINVER": "0.10.0"
    },
    "snippy": {
        "GETVER": "%s/snippy --version 2>&1|grep 'snippy'" % SNIPPY_BIN,
        "REGEXP": "\d+\.\d+\.\d+",
        "MINVER": "4.6.0"
    },
}


DATABASE = "/nextomics/Pipeline/NextnCoV/v1.0.0/database/"
