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
import time

item = 'SRR8368696'
rep = str(PurePath('sequences', item))
repitem = str(PurePath('sequences', item, item))
p_shuffled = str(PurePath(rep, item + '_shuffled.fasta'))
P_FASTA = str(PurePath('tmp', 'snp.fasta'))
debut = time.time()*1000
result = subprocess.run(["blastp", "-num_threads",
                                                     "8", "-query", P_FASTA,
                                                     "-evalue",
                                                     "1e-5", "-task", "blastn",
                                                     "-db", repitem, "-outfmt",
                                                     "10 sseq"],
                                                    stdout=subprocess.PIPE)
fin = time.time()*1000
print(fin-debut)

print(type(result))

"""
with open('essai.txt', 'r') as f_in, open('essai2.txt', 'w') as f_out:
    lignes = f_in.readlines()
    cpt = 0
    for elt in lignes:
        cpt += elt.count('>')
    f_out.write(str(cpt))
"""

