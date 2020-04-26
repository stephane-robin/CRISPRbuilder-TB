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

chaine_csv = input('Please apply the following format when adding a new '
                       'line to the lineage.csv database:\n'
                       '1. use comas between the different fields\n'
                       '2. don\'t use apostrophe or quote marks around the elemts '
                       'of the field\n'
                       '3. fill up the different fields in this order:\n'
                       'lineage, Position, Gene coord., Allele change, Codon number,'
                       'Codon change, Amino acid change, Locus Id, Gene name, Gene '
                       'type, Type of mutation, 5\' gene, 3\' gene, Strand, '
                       'Sublineage surname, Essential, Origin\n')
chaine_csv = chaine_csv.strip()
liste_csv = chaine_csv.split(',')
liste_csv = [u.strip() for u in liste_csv]
if len(liste_csv) > 17:
    print('The line you wrote doesn\'t match the number of fields in '
              'lineage.csv. Please proceed again.')
else:
    if len(liste_csv) < 17:
        liste_tmp = []
        for i in range(17 - len(liste_csv)):
            liste_tmp.append(' ')
        liste_csv.extend(liste_tmp)

    with open('snps.csv', 'a', newline='') as f:
        c = csv.writer(f, delimiter=',', quotechar=' ',
                           quoting=csv.QUOTE_MINIMAL)
        c.writerow(liste_csv)




