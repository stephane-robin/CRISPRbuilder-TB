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
import time




debut = time.time()*1000

fin = time.time()*1000
print('time: ', fin - debut)


"""
def ff():
    demi_longueur = 20
    Lignee_renvoyee = {}

    with open('data/lineage.csv', 'rt') as f:
        csv_reader = csv.reader(f, delimiter=',', quotechar='"')
        next(csv_reader)
        for row in csv_reader:
            if row[1] != '':
                lignee = row[0].strip()
                pos0 = int(row[1].strip())
                pos = pos0 - 1
                source = row[3].strip()[0]
                cible = row[3].strip()[2]
                assert h37Rv[pos] == source
                seq1 = h37Rv[pos - demi_longueur:pos + demi_longueur + 1]
                seq2 = seq1[:demi_longueur] + cible + seq1[demi_longueur + 1:]
                Lignee_renvoyee[pos] = (seq1, seq2, lignee)

    print("We have selected specific reads to compare with different lineages")
    return Lignee_renvoyee

print(ff())
"""


