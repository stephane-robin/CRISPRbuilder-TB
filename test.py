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
from pathlib import PurePath


dico_afr = {}
item = 'SRR8368689'
dico_afr[item] = {}
repitem = 'REP/sequences/' + item + '/' + item
p_shuffled = str(PurePath('REP', 'sequences', item, item
                              + '_shuffled.fasta'))


print("We're creating a database for Blast")
completed = subprocess.run(['makeblastdb', '-in', p_shuffled,
                                        '-dbtype', 'nucl', '-title', item,
                                        '-out', repitem])
assert completed.returncode == 0


print(f"The spoligotypes are being blasted")
dico_afr[item]['spoligo'] = ''
dico_afr[item]['spoligo_new'] = ''

p_spoligo_old = str(PurePath('data',
                                         'spoligo_old.fasta'))
p_spoligo_new = str(PurePath('data',
                                         'spoligo_new.fasta'))
p_old_blast = str(PurePath('/tmp/' + item +
                                       "_old.blast"))
p_new_blast = str(PurePath('/tmp/' + item +
                                       "_new.blast"))

completed = subprocess.run("blastn -num_threads 12 -query " +
                                       p_spoligo_old + " -evalue 1e-6 -task "
                                       "blastn -db " + repitem + " -outfmt '10 "
                                       "qseqid sseqid sstart send qlen length "
                                       "score evalue' -out " + p_old_blast,
                                       shell=True)
assert completed.returncode == 0

completed = subprocess.run("blastn -num_threads 12 -query " +
                                       p_spoligo_new + " -evalue 1e-6 -task "
                                       "blastn -db " + repitem + " -outfmt '10 "
                                       "qseqid sseqid sstart send qlen length "
                                       "score evalue' -out " + p_new_blast,
                                       shell=True)
assert completed.returncode == 0

#print("We're writing the spoligotypes obtained in the csv file")

for pos, spol in enumerate(['old', 'new']):
    p_blast = str(PurePath('/tmp/' + item + '_' +
                                       spol + '.blast'))
    p_fasta = str(PurePath('data', 'spoligo_' +
                                       spol + '.fasta'))

    with open(p_blast) as f:
        matches = f.read()
        nb = open(p_fasta).read().count('>')
        for k in range(1, nb + 1):
            #if 'espaceur'+spol.capitalize()+str(k) in matches:
            #if matches.count('espaceur'+spol.capitalize()+str(k)+',')/dico_afr[item]['couverture']>0.05:
            if matches.count('espaceur' + spol.capitalize() + str(k)
                                         + ',') >= 5:
                dico_afr[item]['spoligo' + ['', '_new'][pos]] \
                                    += '\u25A0'
            else:
                dico_afr[item]['spoligo' + ['', '_new'][pos]] \
                                    += '\u25A1'

    dico_afr[item]['spoligo' + ['', '_new'][pos] + '_nb'] = [
                        matches.count('espaceur' + spol.capitalize() + str(k) +
                                      ',') for k in range(1, nb + 1)]


print("     " + dico_afr[item]['spoligo'])
print("     " + str(dico_afr[item]['spoligo_nb']))
print("     " + dico_afr[item]['spoligo_new'])
print("     " + str(dico_afr[item]['spoligo_new_nb']))




