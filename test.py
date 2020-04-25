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
import csv
from collections import namedtuple




wwb = open_workbook('../data/Coll_62_SNPs_(copie).xlsx')
wws = wwb.sheet_by_index(0)

for row in range(1, wws.nrows):

    chaine = ''
    for i in range(0, 16):
        chaine += str(wws.cell_value(row, i)).replace('.0', '')+', '
        liste_SRA = chaine.strip().split(',')
    with open('essai.csv', 'w', newline='') as file:
        c = csv.writer(file, delimiter=',', quotechar='"',
                       quoting=csv.QUOTE_MINIMAL)
        c.writerow(liste_SRA)

"""
wb = load_workbook(filename='../data/Coll_62_SNPs_(copie).xlsx',
                    read_only=True)
fiche = wb['Feuil1']

with open('../data/Coll_62_SNPs_(copie).xlsx', 'rb') as f:
    chaine_SRA = f.read()
    print(chaine_SRA)
    liste_SRA = chaine_SRA.strip().split()

for row in fiche.iter_rows(min_row=2):
    if row[1].value != None:


        with open('essai.csv', 'w', newline='') as file:
            c = csv.writer(file, delimiter=',', quotechar='"',
                                    quoting=csv.QUOTE_MINIMAL)
            c.writerow(row)
"""


