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



def cat(p1, p2, p_shuffled):
    with open(p1, 'r') as f1, open(p2, 'r') as f2, open(p_shuffled, 'w') as \
            f_shuffled:
        lignes1 = f1.readlines()
        for elt in lignes1:
            f_shuffled.write(elt)
        lignes2 = f2.readlines()
        for elt in lignes2:
            f_shuffled.write(elt)

p1 = str(PurePath('essai.txt'))
p2 = str(PurePath('essai2.txt'))
p_shuffled = str(PurePath('essai_shuffled.txt'))

cat(p1, p2, p_shuffled)

