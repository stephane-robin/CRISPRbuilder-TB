import xmltodict  # interprets XML files like JSON files
import subprocess  # allows to connect to input/output/error pipes of processes
from os import mkdir, chdir, remove, listdir, rename, getcwd, system, \
    stat  # system
from pickle import load, dump
from xlrd import open_workbook
from openpyxl import load_workbook
import csv
import shutil  # TESTER UTILITE D'IMPORTER TOUT LE MODULE
from shutil import copyfile
from Bio import pairwise2, SeqIO
from Bio.pairwise2 import format_alignment
from Bio import Entrez  # provides code to access NCBI over the Web
from datetime import datetime

item = 'SRR8368689'
rep = '../REP/sequences/'+item+'/'
completed = subprocess.run(['makeblastdb', '-in', rep+item+'_shuffled.fasta',
                            '-dbtype', 'nucl','-title', item, '-out', rep+item])
