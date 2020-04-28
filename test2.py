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

with open('data/lineage2.csv', 'r', newline='') as csvin, \
        open('data/lineage3.csv', 'w', newline='') as csvout:
    csv_reader = csv.reader(csvin, delimiter=',', quotechar='"',
                            quoting=csv.QUOTE_MINIMAL)
    csv_writer = csv.writer(csvout, delimiter=',',
                            quotechar='"',
                            quoting=csv.QUOTE_MINIMAL)
    for row in csv_reader:
        if row[0].strip() != '4.1.3':
            csv_writer.writerow(row)
