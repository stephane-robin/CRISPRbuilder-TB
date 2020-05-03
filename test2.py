import pathlib

import xmltodict # interprets XML files like JSON files
import subprocess # allows to connect to input/output/error pipes of processes
from os import mkdir, chdir, remove, listdir, rename, getcwd, system, stat
from pickle import load, dump
from xlrd import open_workbook
from openpyxl import load_workbook
import csv
import shutil
from shutil import copyfile
from Bio import pairwise2, SeqIO
from Bio.pairwise2 import format_alignment
from Bio import Entrez # provides code to access NCBI over the Web
from datetime import datetime
import argparse
import csv
from collections import namedtuple
from pathlib import PurePath
from pathlib import Path
import pathlib
import tempfile

with tempfile.TemporaryDirectory() as fp:
    p = str(PurePath(fp, 'nb.txt'))
    #with open(p, 'w') as f_nb:
        #pass
    # TODO check  /tmp/nb.txt
    # TODO system("cat " + p_shuffled + " | grep '>' | wc -l > " + p)
    # TODO ??? why keep nb in a temp file /tmp/nb.txt ?
    with open('essai.txt', 'r') as f_in, open(p, 'w') as f_out:
        lignes = f_in.readlines()
        cpt = 0
        for elt in lignes:
            cpt += elt.count('>')
        f_out.write(str(cpt))

    nb = eval(open(p).read().split('\n')[0])
    print(nb)

"""
with open('essai.txt', 'r') as f_in, open('essai2.txt', 'w') as f_out:
    lignes = f_in.readlines()
    cpt = 0
    for elt in lignes:
        cpt += elt.count('>')
    f_out.write(str(cpt))
"""

